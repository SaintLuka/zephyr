#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using namespace geom;
using namespace smf;

using mesh::AmrStorage;
using zephyr::utils::threads;

static const SmFluid::State U = SmFluid::datatype();

SmFluid::State SmFluid::datatype() {
    return {};
}

SmFluid::SmFluid(Eos::Ptr eos) : m_eos(eos) {
    m_nf = HLLC::create();
    m_CFL = 0.5;
    m_limiter = Limiter("MC");
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void SmFluid::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
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

PState get_current_sm(Cell &cell) {
    return cell(U).get_state();
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

void SmFluid::update(Mesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        fluxes(mesh);
    }
    else {
        compute_grad(mesh, get_current_sm);

        fluxes_stage1(mesh);
        
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void SmFluid::compute_dt(Mesh &mesh) {
    double dt = mesh.min([this](Cell &cell) -> double {
        double dt = std::numeric_limits<double>::infinity();

        // скорость звука
        double c = m_eos->sound_speed_rP(cell(U).density, cell(U).pressure);

        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell(U).velocity.dot(face.normal());

            // Максимальное по модулю СЗ
            double lambda = std::abs(vn) + c;

            // Условие КФЛ
            dt = std::min(dt, cell.volume() / (face.area() * lambda));
        }

        return dt;
    });

    m_dt = std::min(m_CFL * dt, m_max_dt);
}

void SmFluid::compute_grad(Mesh &mesh, const std::function<smf::PState(Cell &)> &get_state)  {
    mesh.for_each([this, &get_state](Cell &cell) {
        //auto grad = gradient::gauss<smf::PState>(cell, get_state, boundary_value);
        auto grad = gradient::LSM<smf::PState>(cell, get_state, boundary_value);
        grad = gradient::limiting<smf::PState>(cell, m_limiter, grad, get_state, boundary_value);

        cell(U).d_dx = grad.x;
        cell(U).d_dy = grad.y;
        cell(U).d_dz = grad.z;
    });
}

void SmFluid::fluxes(Mesh &mesh) {
    mesh.for_each([this](Cell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

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
                z_n = face.neib(U).get_state();
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
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Новое значение примитивных переменных
        cell(U).next = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_stage1(Mesh &mesh)  {
    mesh.for_each([this](Cell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto &normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа
            PState z_n;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib(U).get_state();
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
            }

            auto face_extra = FaceExtra::Direct(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
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
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).half = PState(q_c, *m_eos);
    });
}

void SmFluid::fluxes_stage2(Mesh &mesh)  {
    mesh.for_each([this](Cell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор на полуслое
        PState z_c = cell(U).get_state();

        // Примитивный вектор на полуслое
        PState z_ch = cell(U).half;

        // Консервативный вектор в ячейке на прошлом шаге
        QState q_c(z_c);

        // Переменная для потока (суммирование по промежуточным)
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto &normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа (на предыдущем и на полуслое)
            PState z_n, z_nh;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib(U).get_state();
                z_nh = neib(U).half;
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
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
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на новом слое
        cell(U).next = PState(q_c, *m_eos);
    });
}

void SmFluid::swap(Mesh &mesh) {
    mesh.for_each([](Cell &cell) {
        cell(U).set_state(cell(U).next);
    });
}

Distributor SmFluid::distributor(const std::string& type) const {
    if (type != "const" && type != "slope") {
        throw std::runtime_error("SmFluid error: unknown distributor type '" + type + "'");
    }
    
    using mesh::Children;
    
    Distributor distr;

    // Консервативное суммирование
    distr.merge = [this](Children &children, AmrStorage::Item &parent) {
        QState q_p;
        for (auto &child: children) {
            QState q_ch(child(U).get_state());
            q_p.arr() += q_ch.arr() * child.volume();
        }
        q_p.arr() /= parent.volume();
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
}

void SmFluid::set_flags(Mesh &mesh) {
    if (!mesh.is_adaptive()) {
        return;
    }

    
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
}

} // namespace zephyr::math