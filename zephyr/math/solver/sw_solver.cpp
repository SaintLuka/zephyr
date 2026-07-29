#include <zephyr/math/solver/sw_solver.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/math/funcs.h>

namespace zephyr::math {

using namespace geom;
using namespace swe;

using utils::threads;
using utils::mpi;

SwSolver::SwSolver(double bed)
    : SwSolver(ConstBed::create(bed)) {

}

SwSolver::SwSolver(IBed::Ref bed) {
    m_bed = bed;
    if (!bed) {
        std::cerr << "Nullptr bed, set const -1 level\n";
        m_bed = ConstBed::create(-1.0);
    }
    m_nf = HLL::create();
    m_CFL = 0.5;
    m_limiter = Limiter("MC");
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

SwSolver::Parts SwSolver::add_types(EuMesh& mesh) {
    part.init = mesh.add<PState>("init");
    part.d_dx = mesh.add<PState>("d_dx");
    part.d_dy = mesh.add<PState>("d_dy");
    part.half = mesh.add<PState>("half");
    part.next = mesh.add<PState>("next");
    part.bed  = mesh.add<double>("bed");
    part.slope= mesh.add<Vector2d>("slope");

    mesh.for_each([this](EuCell& cell) {
        auto[bed, slope] = m_bed->get(cell.center());
        cell[part.bed] = bed;
        cell[part.slope] = slope;
    });
    return part;
}

void SwSolver::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void SwSolver::set_accuracy(int acc) {
    if (acc < 1 || acc > 2) {
        std::cerr << "SmSolver warning: available accuracy 1 or 2\n";
    }
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
}

void SwSolver::set_method(Fluxes method) {
    m_nf = NumFlux::create(method);
}

void SwSolver::set_limiter(const std::string& lim) {
    m_limiter = Limiter(lim);
}

double SwSolver::CFL() const {
    return m_CFL;
}

double SwSolver::dt() const {
    return m_dt;
}

void SwSolver::set_max_dt(double dt) {
    m_max_dt = dt;
}

PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
    if (flag != Boundary::WALL) {
        return zc;
    }

    PState zn(zc);
    Vector2d Vn = normal.head<2>() * zc.velocity.dot(normal.head<2>());
    zn.velocity = zc.velocity - 2.0 * Vn;

    return zn;
}

void SwSolver::update(EuMesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        mesh.sync(part.init);
        fluxes(mesh);
    }
    else {
        mesh.sync(part.init);
        compute_grad(mesh);

        mesh.sync(part.d_dx, part.d_dy);
        fluxes_stage1(mesh);

        mesh.sync(part.half);
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void SwSolver::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        double h = cell[part.init].depth;
        if (h < swe::min_depth) {
            // Почему?
            return std::numeric_limits<double>::max();
        }
        double c = std::sqrt(swe::g * h);
        Vector2d v = cell[part.init].velocity;
        return std::min(cell.hx(), cell.hy()) / (v.norm() + c);

        // Условие строго из Куликовского
        //return 1.0 / ((std::abs(v.x()) + c) / cell.hx() + (std::abs(v.y()) + c) / cell.hy());
    });

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void SwSolver::compute_grad(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        auto grad = gradient::LSM<PState>(cell, part.init, boundary_value);
        grad = gradient::limiting<PState>(cell, m_limiter, grad, part.init, boundary_value);

        cell[part.d_dx] = grad.x;
        cell[part.d_dy] = grad.y;
    });
}

void SwSolver::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];
        double bed_c = cell[part.bed];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto face: cell.faces()) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            double bed_n = bed_c;
            if (!face.is_boundary()) {
                z_n = face.neib(part.init);
                bed_n = face.neib(part.bed);
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Уровень дна на грани
            double bed_f = std::max(bed_c, bed_n);

            // гидростатическая реконструкция Audusse & Bristeau
            zm.depth = std::max(0.0, z_c.surf(bed_c) - bed_f);
            zp.depth = std::max(0.0, z_n.surf(bed_n) - bed_f);

            // Численный поток на грани
            Flux loc_flux = m_nf->flux(zm, zp);

            // Поправка Audusse & Bristeau (well-balanced схема)
            loc_flux.momentum.x() += 0.5 * g * (std::pow(z_c.depth, 2) - std::pow(zm.depth, 2));
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Новое значение примитивных переменных
        cell[part.next] = PState(q_c);
    });
}

void SwSolver::fluxes_stage1(EuMesh &mesh) const {
    /*
    mesh.for_each([this](EuCell &cell) {
        // Ячейка в изоляции
        if (cell[part.wait] > 0) {
            cell[part.half] = cell[part.init];
            return;
        }

        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = neib[part.init];
            }
            else {
                z_n = boundary_value(z_c, normal, face.flag());
            }
            Vector3d neib_c = (face.flag() == Boundary::INNER ?
                face.neib_center() : face.symm_point(cell_c));

            auto face_extra = FaceExtra::Direct(
                    z_c, cell[part.d_dx], cell[part.d_dy], cell[part.d_dz],
                    z_n, neib[part.d_dx], neib[part.d_dy], neib[part.d_dz],
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_c);

            // Восстанавливаем после интерполяции
            zm.energy = m_eos->energy_rP(zm.density, zm.pressure);

            // При некорректной интерполяции
            if (zm.is_bad(*m_eos)) { zm = z_c; }

            // Переводим в локальную систему координат
            zm.to_local(normal);

            // Численный поток на грани
            Flux loc_flux(zm);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area(m_axial);
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * m_dt / cell.volume(m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = cell.volume() / cell.volume(m_axial);
            q_c.momentum.y() += 0.5 * coeff * z_c.pressure * m_dt;
        }

        // Значение примитивных переменных на полушаге
        cell[part.half] = PState(q_c, *m_eos);
        if (cell[part.half].is_bad(*m_eos)) {
            cell[part.half] = z_c;
        }
    });
    */
}

void SwSolver::fluxes_stage2(EuMesh &mesh) const {
    /*
    mesh.for_each([this](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор на полуслое
        PState z_c = cell[part.init];

        // Примитивный вектор на полуслое
        PState z_ch = cell[part.half];

        // Консервативный вектор в ячейке на прошлом шаге
        QState q_c(z_c);

        // Переменная для потока (суммирование по промежуточным)
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto  normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа (на предыдущем и на полушаге)
            PState z_n, z_nh;
            if (!face.is_boundary()) {
                z_n  = neib[part.init];
                z_nh = neib[part.half];
            }
            else {
                z_n  = boundary_value(z_c,  normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            PState zm, zp;

            // Одна из ячеек в изоляции
            bool robust = cell[part.wait] > 0 || neib[part.wait] > 0;
            if (robust) {
                zm = z_ch;
                zp = z_nh;
            }
            else {
                Vector3d neib_c = (face.flag() == Boundary::INNER ?
                    face.neib_center() : face.symm_point(cell_c));

                // Параметры интерполяции с предыдущего (!) слоя
                auto face_extra = FaceExtra::Direct(
                        z_c, cell[part.d_dx], cell[part.d_dy], cell[part.d_dz],
                        z_n, neib[part.d_dx], neib[part.d_dy], neib[part.d_dz],
                        cell_c, neib_c, face_c);

                // Интерполяция на грань со стороны ячейки
                zm = face_extra.m(z_ch);

                // Восстанавливаем после интерполяции
                zm.energy = m_eos->energy_rP(zm.density, zm.pressure);

                // При некорректной интерполяции
                if (zm.is_bad(*m_eos)) { zm = z_ch; }

                // Интерполяция на грань со стороны соседа
                if (!face.is_boundary()) {
                    zp = face_extra.p(z_nh);

                    // Восстанавливаем после интерполяции
                    zp.energy = m_eos->energy_rP(zp.density, zp.pressure);

                    // При некорректной интерполяции
                    if (zp.is_bad(*m_eos)) { zp = z_nh; }
                }
                else {
                    zp = boundary_value(zm, normal, face.flag());
                }
            }

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            auto loc_flux = robust? HLL::calc_flux(zm, zp, *m_eos) : m_nf->flux(zm, zp, *m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area(m_axial);
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume(m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = cell.volume() / cell.volume(m_axial);
            q_c.momentum.y() += coeff * z_ch.pressure * m_dt;
        }

        // Значение примитивных переменных на новом слое
        cell[part.next] = PState(q_c, *m_eos);
        if (cell[part.next].is_bad(*m_eos)) {
            // Сажаем ячейку в изоляцию
            cell[part.next] = z_c;
            cell[part.wait] = 3;
        }
        else {
            // Уменьшаем срок изоляции
            if (cell[part.wait] > 0) {
                --cell[part.wait];
            }
        }
    });
    */
}

void SwSolver::swap(EuMesh &mesh) const {
    mesh.swap(part.init, part.next);
}

Distributor SwSolver::distributor(const std::string& type) const {
    return Distributor::empty();
    /*
    if (type != "const" && type != "slope") {
        throw std::runtime_error("SwSolver error: unknown m_distributor type '" + type + "'");
    }
    
    using mesh::Children;
    
    Distributor distr;

    // Консервативное суммирование
    distr.merge = [this](const Children &children, EuCell &parent) {
        QState q_p;
        for (auto child: children) {
            QState q_ch(child[part.init]);
            q_p.arr() += q_ch.arr() * child.volume();
        }
        q_p.arr() /= parent.volume();
        PState z_p(q_p, *m_eos);
        parent[part.init] = z_p;
        parent[part.wait] = 0;
    };

    // Снос копированием
    auto split_const = [this](const EuCell &parent, Children &children) {
        const PState& z_p = parent[part.init];
        for (auto child: children) {
            child[part.init] = z_p;
        }
    };
    
    // Снос по градиентам
    auto split_slope = [this](const EuCell &parent, Children &children) {
        const PState& z_p  = parent[part.init];
        const PState& d_dx = parent[part.d_dx];
        const PState& d_dy = parent[part.d_dy];
        const PState& d_dz = parent[part.d_dz];

        auto P = m_eos->pressure_re(z_p.density, z_p.energy, {.deriv = true});

        Vector3d grad_e = {
                (d_dx.pressure - P.dR * d_dx.density) / P.dE,
                (d_dy.pressure - P.dR * d_dy.density) / P.dE,
                (d_dz.pressure - P.dR * d_dz.density) / P.dE
        };

        bool bad_grad = false;
        for (auto child: children) {
            Vector3d dr = child.center() - parent.center();

            PState z_ch = parent[part.init];

            z_ch.density = z_p.density +
                           d_dx.density * dr.x() +
                           d_dy.density * dr.y() +
                           d_dz.density * dr.z();

            z_ch.velocity = z_p.velocity +
                            (z_p.density / z_ch.density) * (
                                    d_dx.velocity * dr.x() +
                                    d_dy.velocity * dr.y() +
                                    d_dz.velocity * dr.z()
                            );

            z_ch.energy = z_p.energy - 0.5 * (z_p.velocity - z_ch.velocity).squaredNorm() +
                          (z_p.density / z_ch.density) * grad_e.dot(dr);

            z_ch.pressure = m_eos->pressure_re(z_ch.density, z_ch.energy);

            if (z_ch.is_bad(*m_eos)) {
                bad_grad = true;
                break;
            }

            child[part.init] = z_ch;
            child[part.wait] = 0;
        }

        // Не удалось сделать интерполяцию в одну из дочерних ячеек,
        // выполняем простой перенос
        if (bad_grad) {
            for (auto child: children) {
                child[part.init] = z_p;
                child[part.wait] = 0;
            }
        }
    };
    
    if (type == "const") {
        // Снос копированием
        distr.split = split_const;
    }
    else {
        // Снос по градиентам
        distr.split = split_slope;
    }

    return distr;
    */
}

#if 0
void SwSolver::set_flags(EuMesh &mesh) const {
    if (!mesh.adaptive()) { return; }

    compute_grad(mesh);

    // Пороги (относительные) на разбиение
    const double xi_dens = 0.05;
    const double xi_pres = 0.05;

    for (auto cell: mesh) {
        //cell.set_flag(1); continue;
        cell.set_flag(-1);

        double dens = cell[part.init].density;
        double pres = cell[part.init].pressure;

        double dens_split = xi_dens * std::abs(dens);
        double pres_split = xi_pres * std::abs(pres);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            double dens_n = face.neib(part.init).density;
            double pres_n = face.neib(part.init).pressure;

            // Большой перепад плотностей или давлений
            if (std::abs(dens_n - dens) > dens_split ||
                std::abs(pres_n - pres) > pres_split) {
                cell.set_flag(1);
                break;
            }

            // Пороги минимум в два раза меньше
            if (std::abs(dens_n - dens) > 0.4 * dens_split ||
                std::abs(pres_n - pres) > 0.4 * pres_split) {
                cell.set_flag(0);
            }
        }
    }
}
#endif

void SwSolver::set_flags(EuMesh &mesh) const {
#if 0
    if (!mesh.adaptive()) { return; }

    mesh.sync(part.init);
    compute_grad(mesh);


    // Пороги (относительные) на разбиение
    const double xi_dens = 0.05;
    const double xi_pres = 0.05;

    for (auto cell: mesh) {
        cell.set_flag(-1);

        // ---------------------- SLOPE CRITERION ---------------------------

        double dens = cell[part.init].density;
        double pres = cell[part.init].pressure;

        double dens_split = xi_dens * std::abs(dens);
        double pres_split = xi_pres * std::abs(pres);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            double dens_n = face.neib(part.init).density;
            double pres_n = face.neib(part.init).pressure;

            // Большой перепад плотностей или давлений
            if (std::abs(dens_n - dens) > dens_split ||
                std::abs(pres_n - pres) > pres_split) {
                cell.set_flag(1);
                break;
            }

            // Пороги минимум в два раза меньше
            if (std::abs(dens_n - dens) > 0.4 * dens_split ||
                std::abs(pres_n - pres) > 0.4 * pres_split) {
                cell.set_flag(0);
            }
        }
        /*
        // ---------------------- CHI CRITERION ---------------------------
        const auto& zc = cell[part.init];
        const auto& dzcx = cell[part.d_dx];
        const auto& dzcy = cell[part.d_dy];
        const auto& dzcz = cell[part.d_dz];

        Matrix3d dens_A = Matrix3d::Zero();
        Matrix3d dens_B = Matrix3d::Zero();
        Matrix3d pres_A = Matrix3d::Zero();
        Matrix3d pres_B = Matrix3d::Zero();

        double full_area = 0.0;
        for(auto& face: cell.faces()) {
            full_area += face.area();
        }

        const double eps = 0.001;
        Vector3d cell_c = cell.center();

        for(auto& face: cell.faces()) {
            auto neib = face.neib();
            Vector3d normal = face.normal();
            Vector3d neig_c = neib.center();
            Vector3d face_c = face.center();

            const auto& zn = neib[part.init];
            const auto& dznx = neib[part.d_dx];
            const auto& dzny = neib[part.d_dy];
            const auto& dznz = neib[part.d_dz];

            double S  = face.area();
            auto   Sn = normal * S;

            // Значения на гранях
            Vector3d drc = face_c - cell_c;
            PState zf = zc.arr() + dzcx.arr() * drc.x() + dzcy.arr() * drc.y() + dzcz.arr() * drc.z();

            //linear interpolation for derivatives at edge
            double t = (face_c - cell_c).dot(normal);
            t       /= (neig_c - cell_c).dot(normal);

            std::array<PState, 3> dzf = {
                dzcx.arr() + t * (dznx.arr() - dzcx.arr()),
                dzcy.arr() + t * (dzny.arr() - dzcy.arr()),
                dzcz.arr() + t * (dznz.arr() - dzcz.arr())
            };

            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    dens_A(i, j) += 0.5 * (dzf[i].density * Sn[j] + dzf[j].density * Sn[i]);
                    dens_B(i, j) += 0.5 * (fabs(dzf[i].density) + fabs(dzf[j].density)) * S;
                    dens_B(i, j) += eps * fabs(zf.density) * S / full_area;

                    pres_A(i, j) += 0.5 * (dzf[i].pressure * Sn[j] + dzf[j].pressure * Sn[i]);
                    pres_B(i, j) += 0.5 * (fabs(dzf[i].pressure) + fabs(dzf[j].pressure)) * S;
                    pres_B(i, j) += eps * fabs(zf.pressure) * S / full_area;
                }
            }
        }

        double dens_norm_a = dens_A.squaredNorm();
        double dens_norm_b = dens_B.squaredNorm() + 1.e-10;
        double pres_norm_a = pres_A.squaredNorm();
        double pres_norm_b = pres_B.squaredNorm() + 1.e-10;

        double dens_chi = std::sqrt( dens_norm_a / dens_norm_b );
        double pres_chi = std::sqrt( pres_norm_a / pres_norm_b );

        const double chi_p = 0.15; // Верхний порог
        const double chi_m = 0.10; // Нижний порог

        if (dens_chi > chi_p || pres_chi > chi_p) {
            cell.set_flag(1);
        }
        else if (dens_chi < chi_m && pres_chi < chi_m) {
            cell.set_flag(-1);
        }
        else {
            cell.set_flag(0);
        }
        */
    }
#endif
}

} // namespace zephyr::math