#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>

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

SmFluid::State SmFluid::add_types(EuMesh& mesh) {
    /*
    data.density  = mesh.add<double>("density");
    data.velocity = mesh.add<Vector3d>("velocity");
    data.pressure = mesh.add<double>("pressure");
    data.energy   = mesh.add<double>("energy");
    */

    data.init = mesh.add<PState>("init");
    data.half = mesh.add<PState>("half");
    data.next = mesh.add<PState>("next");
    data.d_dx = mesh.add<PState>("d_dx");
    data.d_dy = mesh.add<PState>("d_dy");
    data.d_dz = mesh.add<PState>("d_dz");

    return data;
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void SmFluid::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
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
        mesh.sync(data.init);
        fluxes(mesh);
    }
    else {
        mesh.sync(data.init);

        compute_grad(mesh);

        mesh.sync(data.d_dx, data.d_dy, data.d_dz);
        fluxes_stage1(mesh);

        mesh.sync(data.half);
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void SmFluid::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        //double c = m_eos->sound_speed_rP(cell(data.density), cell(data.pressure));
        //return cell.incircle_diameter() / (cell(data.velocity).norm() + c);
        double c = m_eos->sound_speed_rP(cell(data.init).density, cell(data.init).pressure);
        return cell.incircle_diameter() / (cell(data.init).velocity.norm() + c);
    }, 1.0e300);

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void SmFluid::compute_grad(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        auto grad = gradient::LSM<PState>(cell, data.init, boundary_value);
        grad = gradient::limiting<PState>(cell, m_limiter, grad, data.init, boundary_value);

        cell(data.d_dx) = grad.x;
        cell(data.d_dy) = grad.y;
        cell(data.d_dz) = grad.z;
    });
}

void SmFluid::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell(data.init);

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
                z_n = face.neib(data.init);
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
        cell(data.next) = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_stage1(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        PState z_c = cell(data.init);

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
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib(data.init);
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
            }

            auto face_extra = FaceExtra::Direct(
                    z_c, cell(data.d_dx), cell(data.d_dy), cell(data.d_dz),
                    z_n, neib(data.d_dx), neib(data.d_dy), neib(data.d_dz),
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
        cell(data.half) = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_stage2(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор на полуслое
        PState z_c = cell(data.init);

        // Примитивный вектор на полуслое
        PState z_ch = cell(data.half);

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
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n  = neib(data.init);
                z_nh = neib(data.half);
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n  = boundary_value(z_c,  normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell(data.d_dx), cell(data.d_dy), cell(data.d_dz),
                    z_n, neib(data.d_dx), neib(data.d_dy), neib(data.d_dz),
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
        cell(data.next) = PState(q_c, *m_eos);
    });
}

void SmFluid::swap(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        cell(data.init) = cell(data.next);
    });
}

Distributor SmFluid::distributor(const std::string& type) const {
    if (type != "const" && type != "slope") {
        throw std::runtime_error("SmFluid error: unknown distributor type '" + type + "'");
    }
    
    using mesh::Children;
    
    Distributor distr;

    // Консервативное суммирование
    distr.merge = [this](Children &children, EuCell &parent) {
        QState q_p;
        for (auto child: children) {
            QState q_ch(child(data.init));
            q_p.arr() += q_ch.arr() * child.volume();
        }
        q_p.arr() /= parent.volume();
        PState z_p(q_p, *m_eos);
        parent(data.init) = z_p;
    };

    // Снос копированием
    auto split_const = [this](EuCell &parent, Children &children) {
        PState z_p = parent(data.init);
        for (auto child: children) {
            child(data.init) = z_p;
        }
    };
    
    // Снос по градиентам
    auto split_slope = [this](EuCell &parent, Children &children) {
        PState& z_p  = parent(data.init);
        PState& d_dx = parent(data.d_dx);
        PState& d_dy = parent(data.d_dy);
        PState& d_dz = parent(data.d_dz);

        auto P = m_eos->pressure_re(z_p.density, z_p.energy, {.deriv = true});

        Vector3d grad_e = {
                (d_dx.pressure - P.dR * d_dx.density) / P.dE,
                (d_dy.pressure - P.dR * d_dy.density) / P.dE,
                (d_dz.pressure - P.dR * d_dz.density) / P.dE
        };

        bool bad_grad = false;
        for (auto child: children) {
            Vector3d dr = child.center() - parent.center();

            PState z_ch = parent(data.init);

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

            child(data.init) = z_ch;
        }

        // Не удалось сделать интерполяцию в одну из дочерних ячеек,
        // выполняем простой перенос
        if (bad_grad) {
            for (auto child: children) {
                child(data.init) = z_p;
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

void SmFluid::set_flags(EuMesh &mesh) const {
    if (!mesh.adaptive()) {
        return;
    }

    compute_grad(mesh);

    // Пороги (относительные) на разбиение
    const double xi_dens = 0.05;
    const double xi_pres = 0.05;

    for (auto cell: mesh) {
        //cell.set_flag(1); continue;
        cell.set_flag(-1);

        double dens = cell(data.init).density;
        double pres = cell(data.init).pressure;

        double dens_split = xi_dens * std::abs(dens);
        double pres_split = xi_pres * std::abs(pres);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            double dens_n = face.neib(data.init).density;
            double pres_n = face.neib(data.init).pressure;

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

} // namespace zephyr::math