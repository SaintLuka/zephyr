#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/math/calc/weno.h>

namespace zephyr::math {

using namespace geom;
using namespace smf;

using utils::threads;
using utils::mpi;


SmFluid::SmFluid(Eos::Ptr eos) : m_eos(eos) {
    m_nf = HLLC::create();
    m_CFL = 0.5;
    m_axial = false;
    m_limiter = Limiter("MC");
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

SmFluid::Parts SmFluid::add_types(EuMesh& mesh) {
    part.init = mesh.add<PState>("init");
    part.d_dx = mesh.add<PState>("d_dx");
    part.d_dy = mesh.add<PState>("d_dy");
    part.d_dz = mesh.add<PState>("d_dz");
    part.half = mesh.add<PState>("half");
    part.next = mesh.add<PState>("next");
    return part;
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void SmFluid::set_accuracy(int acc) {
    // m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
    m_acc = acc;
}

void SmFluid::set_axial(bool axial) {
    m_axial = axial;
}

void SmFluid::set_method(Fluxes method) {
    m_nf = NumFlux::create(method);
}

void SmFluid::set_limiter(const std::string& lim) {
    m_limiter = Limiter(lim);
}

double SmFluid::CFL() const {
    return m_CFL;
}

double SmFluid::dt() const {
    return m_dt;
}

void SmFluid::set_max_dt(double dt) {
    m_max_dt = dt;
}

PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
    if (flag != Boundary::WALL) {
        return zc;
    }

    PState zn(zc);
    Vector3d Vn = normal * zc.velocity.dot(normal);
    zn.velocity = zc.velocity - 2.0 * Vn;

    return zn;
}

void SmFluid::update(EuMesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        mesh.sync(part.init);
        fluxes(mesh);
    }
    else if (m_acc == 5) {
        mesh.sync(part.init);
        fluxes_weno(mesh);
    }
    else {
        mesh.sync(part.init);
        compute_grad(mesh);

        mesh.sync(part.d_dx, part.d_dy, part.d_dz);
        fluxes_stage1(mesh);

        mesh.sync(part.half);
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void SmFluid::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        //double c = m_eos->sound_speed_rP(cell(part.density), cell(part.pressure));
        //return cell.incircle_diameter() / (cell(part.velocity).norm() + c);
        double c = m_eos->sound_speed_rP(cell[part.init].density, cell[part.init].pressure);
        return cell.incircle_diameter() / (cell[part.init].velocity.norm() + c);
    });

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void SmFluid::compute_grad(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        auto grad = gradient::LSM<PState>(cell, part.init, boundary_value);
        grad = gradient::limiting<PState>(cell, m_limiter, grad, part.init, boundary_value);

        cell[part.d_dx] = grad.x;
        cell[part.d_dy] = grad.y;
        cell[part.d_dz] = grad.z;
    });
}

void SmFluid::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = face.neib(part.init);
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Численный поток на грани
            Flux loc_flux = m_nf->flux(zm, zp, *m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area(m_axial);
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume(m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = cell.volume() / cell.volume(m_axial);
            q_c.momentum.y() += coeff * z_c.pressure * m_dt;
        }

        // Новое значение примитивных переменных
        cell[part.next] = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_weno(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;

        // Cосед слева от рассматриваемой ячейки (I_{i-1})
        auto face_left = cell.face(mesh::Side2D::L);
        auto normal_left = face_left.normal();
        auto cell_left = face_left.neib();

        // Примитивные вектора соседей слева (включая I_{i-2} и I_{i-3})
        PState neib_left;
        PState neib_one_left;
        PState neib_two_left;

        // Расчет примитивных векторов слева от рассматриваемой ячейки
        if (!face_left.is_boundary()) {
            neib_left = cell_left[part.init];

            // Сосед через одного слева от рассматриваемой ячейки (I_{i-2})
            auto face_one_left = cell_left.face(mesh::Side2D::L);
            auto normal_one_left = face_one_left.normal();
            auto cell_one_left = face_one_left.neib();

            if (!face_one_left.is_boundary()) {
                neib_one_left = cell_one_left[part.init];

                // Сосед через два слева от рассматриваемой ячейки (I_{i-3})
                auto face_two_left = cell_one_left.face(mesh::Side2D::L);
                auto normal_two_left = face_two_left.normal();
                auto cell_two_left = face_two_left.neib();

                if (!face_two_left.is_boundary()) {
                    neib_two_left = cell_two_left[part.init];
                } else {
                    neib_two_left = boundary_value(cell_one_left[part.init], normal_two_left, face_two_left.flag());
                }
            } else {
                neib_one_left = boundary_value(cell_left[part.init], normal_one_left, face_one_left.flag());
            }
        } else {
            neib_left = boundary_value(z_c, normal_left, face_left.flag());
            neib_one_left = boundary_value(z_c, normal_left, face_left.flag());
        }

        // Cосед справа от рассматриваемой ячейки (I_{i+1})
        auto face_right = cell.face(mesh::Side2D::R);
        auto normal_right = face_right.normal();
        auto cell_right = face_right.neib();

        // Примитивные вектора соседей справа (включая I_{i+2} и I_{i+3})
        PState neib_right;
        PState neib_one_right;
        PState neib_two_right;

        // Расчет примитивных векторов справа от рассматриваемой ячейки
        if (!face_right.is_boundary()) {
            neib_right = cell_right[part.init];

            // Сосед через одного справа от рассматриваемой ячейки (I_{i+2})
            auto face_one_right = cell_right.face(mesh::Side2D::R);
            auto normal_one_right = face_one_right.normal();
            auto cell_one_right = face_one_right.neib();

            if (!face_one_right.is_boundary()) {
                neib_one_right = cell_one_right[part.init];
            } else {
                neib_one_right = boundary_value(cell_right[part.init], normal_one_right, face_one_right.flag());
            }
        } else {
            neib_right = boundary_value(z_c, normal_right, face_right.flag());
            neib_one_right = boundary_value(z_c, normal_right, face_right.flag());
            neib_two_right = boundary_value(z_c, normal_right, face_right.flag());
        }

        // Преобразование примитивных векторов в консервативные
        QState neib_left_c(neib_left);
        QState neib_one_left_c(neib_one_left);
        QState neib_two_left_c(neib_two_left);

        QState neib_right_c(neib_right);
        QState neib_one_right_c(neib_one_right);
        QState neib_two_right_c(neib_two_right);

        // Консервативные векторы u_{i+1/2}^{-} и u_{i+1/2}^{+}
        QState q_iplus12_minus;
        QState q_iminus12_plus;

        QState q_iminus12_minus;
        QState q_iplus12_plus;

        // Расчет u_{i+1/2} и u_{i-1/2}
        for (int i = 0; i < q_c.arr().size(); i++) {
            // WENO5 для u_{i+1/2}^{-} и u_{i-1/2}^{+}
            WENO5 stencil_i{neib_one_left_c.arr()[i], neib_left_c.arr()[i],
                                q_c.arr()[i], neib_right_c.arr()[i], neib_one_right_c.arr()[i]};
            q_iplus12_minus.arr()[i] += stencil_i.p();
            q_iminus12_plus.arr()[i] += stencil_i.m();

            // WENO5 для u_{i-1/2}^{-} и u_{i+1/2}^{+}
            WENO5 stencil_i_minus1{neib_two_left_c.arr()[i], neib_one_left_c.arr()[i],
                                neib_left_c.arr()[i], q_c.arr()[i], neib_right_c.arr()[i]};
            q_iminus12_minus.arr()[i] += stencil_i_minus1.p();
            q_iplus12_plus.arr()[i] += stencil_i_minus1.m();
        }

        // Преобразование консервативных векторов в примитивные
        PState z_iplus12_minus{q_iplus12_minus, *m_eos};
        PState z_iminus12_plus{q_iminus12_plus, *m_eos};

        PState z_iminus12_minus{q_iminus12_minus, *m_eos};
        PState z_iplus12_plus{q_iplus12_plus, *m_eos};

        // Численные потоки f_{i+1/2} и f_{i-1/2}
        Flux f_iplus12 = m_nf->flux(z_iplus12_minus.in_local(normal_right), z_iplus12_plus.in_local(normal_right), *m_eos);
        Flux f_iminus12 = m_nf->flux(z_iminus12_minus.in_local(normal_left), z_iminus12_plus.in_local(normal_left), *m_eos);

        f_iplus12.to_global(normal_right);
        f_iminus12.to_global(normal_left);

        // Суммируем потоки
        flux.arr() += f_iplus12.arr() * face_right.area(m_axial);
        flux.arr() += f_iminus12.arr() * face_left.area(m_axial);

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume(m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = cell.volume() / cell.volume(m_axial);
            q_c.momentum.y() += coeff * z_c.pressure * m_dt;
        }

        // Новое значение примитивных переменных
        cell(part.next) = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_stage1(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
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
            Vector3d neib_c = (face.flag() == Boundary::ORDINARY ?
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
    });
}

void SmFluid::fluxes_stage2(EuMesh &mesh) const {
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
            Vector3d neib_c = (face.flag() == Boundary::ORDINARY ?
                face.neib_center() : face.symm_point(cell_c));

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell[part.d_dx], cell[part.d_dy], cell[part.d_dz],
                    z_n, neib[part.d_dx], neib[part.d_dy], neib[part.d_dz],
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_ch);

            // Восстанавливаем после интерполяции
            zm.energy = m_eos->energy_rP(zm.density, zm.pressure);

            // При некорректной интерполяции
            if (zm.is_bad(*m_eos)) { zm = z_ch; }

            // Интерполяция на грань со стороны соседа
            PState zp;
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

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            auto loc_flux = m_nf->flux(zm, zp, *m_eos);
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
    });
}

void SmFluid::swap(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        cell[part.init] = cell[part.next];
    });
}

Distributor SmFluid::distributor(const std::string& type) const {
    if (type != "const" && type != "slope") {
        throw std::runtime_error("SmFluid error: unknown m_distributor type '" + type + "'");
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
    };

    // Снос копированием
    auto split_const = [this](const EuCell &parent, Children &children) {
        PState z_p = parent[part.init];
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
        }

        // Не удалось сделать интерполяцию в одну из дочерних ячеек,
        // выполняем простой перенос
        if (bad_grad) {
            for (auto child: children) {
                child[part.init] = z_p;
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
}

#if 0
void SmFluid::set_flags(EuMesh &mesh) const {
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

void SmFluid::set_flags(EuMesh &mesh) const {
    if (!mesh.adaptive()) { return; }

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
    }
}

} // namespace zephyr::math