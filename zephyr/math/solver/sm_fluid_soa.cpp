#include <zephyr/mesh/primitives/bfaces.h>
#include <zephyr/math/solver/sm_fluid_soa.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using namespace geom;
using namespace smf;

using mesh::AmrStorage;
using zephyr::utils::threads;
using zephyr::utils::mpi;

SmFluidSoA::State SmFluidSoA::datatype() {
    return {};
}

SmFluidSoA::SmFluidSoA(Eos::Ptr eos) : m_eos(eos) {
    m_nf = HLLC::create();
    m_CFL = 0.5;
    m_axial = false;
    m_limiter = Limiter("MC");
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

void SmFluidSoA::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void SmFluidSoA::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
}

void SmFluidSoA::set_axial(bool axial) {
    m_axial = axial;
}

void SmFluidSoA::set_method(Fluxes method) {
    m_nf = NumFlux::create(method);
}

void SmFluidSoA::set_limiter(const std::string& lim) {
    m_limiter = Limiter(lim);
}

double SmFluidSoA::CFL() const {
    return m_CFL;
}

double SmFluidSoA::dt() const {
    return m_dt;
}

void SmFluidSoA::set_max_dt(double dt) {
    m_max_dt = dt;
}

#define get_current_sm ([this](Cell& cell) -> auto { return curr[cell.index()]; })

namespace {
PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
    if (flag != Boundary::WALL) {
        return zc;
    }

    PState zn(zc);
    Vector3d Vn = normal * zc.velocity.dot(normal);
    zn.velocity = zc.velocity - 2.0 * Vn;

    return zn;
}
}

void SmFluidSoA::update(Mesh &eu_mesh, SoaMesh& mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        eu_mesh.sync();
        fluxes(eu_mesh, mesh);
    }
    else {
        eu_mesh.sync();
        compute_grad(mesh);

        fluxes_stage1(mesh);

        eu_mesh.sync();
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void SmFluidSoA::compute_dt(SoaMesh& smesh) {
    double dt = smesh.min(
        [this, &smesh](size_t ic) -> double {
            double c = m_eos->sound_speed_rP(curr[ic].density, curr[ic].pressure);
            return smesh.cell_diameter(ic) / (curr[ic].velocity.norm() + c);
        }, std::numeric_limits<double>::infinity());

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void SmFluidSoA::compute_grad(SoaMesh& mesh)  {
    mesh.for_each([this, &mesh](size_t ic) {
        auto grad = gradient::LSM<smf::PState>(mesh, ic, curr.data(), boundary_value);
        grad = gradient::limiting<smf::PState>(mesh, ic, m_limiter, grad, curr.data(), boundary_value);

        d_dx[ic] = grad.x;
        d_dy[ic] = grad.y;
        d_dz[ic] = grad.z;
    });
}

void SmFluidSoA::fluxes(Mesh &mesh, SoaMesh& smesh) {
    mesh.for_each([this](Cell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = curr[cell.index()];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto &normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = curr[face.neib().index()];
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            if (face.center().y() < 1.0e-10) {
                continue;
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
        next[cell.index()] = PState(q_c, *m_eos);
    });
}

void SmFluidSoA::fluxes_stage1(SoaMesh& mesh)  {
    mesh.for_each([this, &mesh](size_t ic) {
        // Центр ячейки
        auto& cell_c = mesh.cell.center[ic];

        // Примитивный вектор в ячейке
        PState z_c = curr[ic];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (size_t k: mesh.cell_faces(ic)) {
            if (mesh.is_undefined_face(k)) {
                continue;
            }

            // Внешняя нормаль и центр грани
            auto &normal = mesh.face.normal[k];
            auto &face_c = mesh.face.center[k];

            // Возвращает саму ячейку, если соседа не существует
            size_t jc = mesh.face.adjacent[k].index;

            // Примитивный вектор соседа
            PState z_n;
            Vector3d neib_c;
            if (!mesh.is_boundary_face(k)) {
                neib_c = mesh.cell_center(jc);
                z_n    = curr[jc];
            }
            else {
                neib_c = mesh.face_symm_point(k, cell_c);
                z_n = boundary_value(z_c, normal, mesh.face.boundary[k]);
            }

            auto face_extra = FaceExtra::Direct(
                    z_c, d_dx[ic], d_dy[ic], d_dz[ic],
                    z_n, d_dx[jc], d_dy[jc], d_dz[jc],
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
            flux.arr() += loc_flux.arr() * mesh.face_area(k, m_axial);
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * m_dt / mesh.cell_volume(ic, m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = mesh.cell_volume(ic) / mesh.cell_volume(ic, m_axial);
            q_c.momentum.y() += 0.5 * coeff * z_c.pressure * m_dt;
        }

        // Значение примитивных переменных на полушаге
        half[ic] = PState(q_c, *m_eos);
    });
}

void SmFluidSoA::fluxes_stage2(SoaMesh& mesh)  {
    mesh.for_each([this, &mesh](int ic) {
        // Центр ячейки
        auto& cell_c = mesh.cell_center(ic);

        // Примитивный вектор на полуслое
        PState z_c = curr[ic];

        // Примитивный вектор на полуслое
        PState z_ch = half[ic];

        // Консервативный вектор в ячейке на прошлом шаге
        QState q_c(z_c);

        // Переменная для потока (суммирование по промежуточным)
        Flux flux;
        for (size_t k: mesh.cell_faces(ic)) {
            if (mesh.is_undefined_face(k)) {
                continue;
            }

            // Внешняя нормаль и центр грани
            auto& normal = mesh.face.normal[k];
            auto& face_c = mesh.face.center[k];

            // Возвращает саму ячейку, если соседа не существует
            int jc = mesh.face.adjacent[k].index;

            // Примитивный вектор соседа (на предыдущем и на полуслое)
            PState z_n, z_nh;
            Vector3d neib_c;
            if (!mesh.is_boundary_face(k)) {
                neib_c = mesh.cell_center(jc);
                z_n    = curr[jc];
                z_nh   = half[jc];
            }
            else {
                neib_c = mesh.face_symm_point(k, cell_c);
                z_n    = boundary_value(z_c, normal, mesh.face.boundary[k]);
                z_nh   = boundary_value(z_ch, normal, mesh.face.boundary[k]);
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, d_dx[ic], d_dy[ic], d_dz[ic],
                    z_n, d_dx[jc], d_dy[jc], d_dz[jc],
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_ch);

            // Восстанавливаем после интерполяции
            zm.energy = m_eos->energy_rP(zm.density, zm.pressure);

            // При некорректной интерполяции
            if (zm.is_bad(*m_eos)) { zm = z_ch; }

            // Интерполяция на грань со стороны соседа
            PState zp;
            if (!mesh.is_boundary_face(k)) {
                zp = face_extra.p(z_nh);

                // Восстанавливаем после интерполяции
                zp.energy = m_eos->energy_rP(zp.density, zp.pressure);

                // При некорректной интерполяции
                if (zp.is_bad(*m_eos)) { zp = z_nh; }
            }
            else {
                zp = boundary_value(zm, normal, mesh.face.boundary[k]);
            }

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            auto loc_flux = m_nf->flux(zm, zp, *m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * mesh.face_area(m_axial);
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / mesh.cell_volume(ic, m_axial)) * flux.arr();

        if (m_axial) {
            double coeff = mesh.cell_volume(ic) / mesh.cell_volume(ic, m_axial);
            q_c.momentum.y() += coeff * z_ch.pressure * m_dt;
        }

        // Значение примитивных переменных на новом слое
        next[ic] = PState(q_c, *m_eos);
    });
}

void SmFluidSoA::swap(SoaMesh& smesh) {
    std::swap(curr, next);
}

Distributor SmFluidSoA::distributor(const std::string& type) const {
    return Distributor::empty();
    //if (type != "const" && type != "slope") {
    //    throw std::runtime_error("SmFluidSoA error: unknown distributor type '" + type + "'");
    //}
    /*
    using mesh::Children;
    
    Distributor distr;

    // Консервативное суммирование
    distr.merge = [this](Children &children, AmrStorage::Item &parent) {
        QState q_p;
        for (auto &child: children) {
            QState q_ch(child(U).get_state());
            q_p.arr() += q_ch.arr() * child.get_volume();
        }
        q_p.arr() /= parent.get_volume();
        PState z_p(q_p, *m_eos);
        parent(U).set_state(z_p);
    };

    // Снос копированием
    auto split_const = [this](AmrStorage::Item &parent, Children &children) {
        PState z_p = parent(U).get_state();
        for (auto &child: children) {
            child(U).set_state(z_p);
        }
    };
    
    // Снос по градиентам
    auto split_slope = [this](AmrStorage::Item &parent, Children &children) {
        PState z_p = parent(U).get_state();
        PState& d_dx = parent(U).d_dx;
        PState& d_dy = parent(U).d_dy;
        PState& d_dz = parent(U).d_dz;

        auto P = m_eos->pressure_re(z_p.density, z_p.energy, {.deriv = true});

        Vector3d grad_e = {
                (d_dx.pressure - P.dR * d_dx.density) / P.dE,
                (d_dy.pressure - P.dR * d_dy.density) / P.dE,
                (d_dz.pressure - P.dR * d_dz.density) / P.dE
        };

        bool bad_grad = false;
        for (auto &child: children) {
            Vector3d dr = child.center - parent.center;

            PState z_ch = parent(U).get_state();

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

            child(U).set_state(z_ch);
        }

        // Не удалось сделать интерполяцию в одну из дочерних ячеек,
        // выполняем простой перенос
        if (bad_grad) {
            for (auto &child: children) {
                child(U).set_state(z_p);
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

void SmFluidSoA::set_flags(Mesh &mesh) {
    if (!mesh.is_adaptive()) {
        return;
    }

    throw std::runtime_error("Not yet #18");

    /*

    
    compute_grad(mesh, get_current_sm);

    // Пороги (относительные) на разбиение
    const double xi_dens = 0.05;
    const double xi_pres = 0.05;

    for (auto cell: mesh) {
        //cell.set_flag(1); continue;
        cell.set_flag(-1);

        double dens = cell(U).density;
        double pres = cell(U).pressure;

        double dens_split = xi_dens * std::abs(dens);
        double pres_split = xi_pres * std::abs(pres);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            auto& zn = face.neib(U);

            // Большой перепад плотностей или давлений
            if (std::abs(zn.density  - dens) > dens_split ||
                std::abs(zn.pressure - pres) > pres_split) {
                cell.set_flag(1);
                break;
            }

            // Пороги минимум в два раза меньше
            if (std::abs(zn.density  - dens) > 0.4 * dens_split ||
                std::abs(zn.pressure - pres) > 0.4 * pres_split) {
                cell.set_flag(0);
            }
        }
    }
     */
}

} // namespace zephyr::math