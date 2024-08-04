#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>

#include <iomanip>

namespace zephyr::math {

using namespace geom;
using namespace smf;

using mesh::AmrStorage;
using zephyr::utils::threads;

static const SmFluid::State U = SmFluid::datatype();

SmFluid::State SmFluid::datatype() {
    return {};
}

smf::PState get_current_sm(Cell &cell) {
    return cell(U).get_state();
}

SmFluid::SmFluid(const phys::Eos &eos) : m_eos(eos) {
    m_nf = HLLC::create();
    m_CFL = 0.9;
    m_dt = std::numeric_limits<double>::max();
}

void SmFluid::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);    // 1 или 2
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

    m_step += 1;
    m_time += m_dt;
}

void SmFluid::compute_dt(Mesh &mesh) {
    double dt = mesh.min([this](Cell &cell) -> double {
        double dt = std::numeric_limits<double>::infinity();

        // скорость звука
        double c = m_eos.sound_speed_rp(cell(U).density, cell(U).pressure);

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

    m_dt = m_CFL * dt;
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
            Flux loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Новое значение примитивных переменных
        cell(U).next = PState(q_c, m_eos);
    });
}

void SmFluid::compute_grad(Mesh &mesh, const std::function<smf::PState(Cell &)> &get_state)  {
    mesh.for_each([&get_state](Cell &cell) -> void {
        // auto grad = compute_grad_gauss<smf::PState>(cell, get_state);
        auto grad = gradient::LSM<smf::PState>(cell, get_state, boundary_value);

        //auto lim_grad = gradient::limiting<smf::PState>(cell, grad, get_state, boundary_value);

        cell(U).d_dx = grad[0];
        cell(U).d_dy = grad[1];
        cell(U).d_dz = grad[2];
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

            auto face_extra = FaceExtra::ATvL(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_c);

            // Интерполяция на грань со стороны соседа
            PState zp;
            if (!face.is_boundary()) {
                zp = face_extra.p(z_n);
            }
            else {
                zp = boundary_value(zm, normal, face.flag());
            }

            // Восстанавливаем после интерполяции
            zm.energy = m_eos.energy_rp(zm.density, zm.pressure);
            zp.energy = m_eos.energy_rp(zp.density, zp.pressure);

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            Flux loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).half = PState(q_c, m_eos);
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
            auto face_extra = FaceExtra::ATvL(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_ch);

            // Интерполяция на грань со стороны соседа
            PState zp;
            if (!face.is_boundary()) {
                zp = face_extra.p(z_nh);
            }
            else {
                zp = boundary_value(zm, normal, face.flag());
            }

            // Восстанавливаем после интерполяции
            zm.energy = m_eos.energy_rp(zm.density, zm.pressure);
            zp.energy = m_eos.energy_rp(zp.density, zp.pressure);

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            auto loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).next = PState(q_c, m_eos);
    });
}

void SmFluid::swap(Mesh &mesh) {
    mesh.for_each([](Cell &cell) {
        cell(U).set_state(cell(U).next);
    });
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

double SmFluid::dt() const {
    return m_dt;
}

double SmFluid::CFL() const {
    return m_CFL;
}

double SmFluid::get_time() const {
    return m_time;
}

size_t SmFluid::get_step() const {
    return m_step;
}

std::string SmFluid::get_flux_name() const {
    return m_nf->get_name();
}

void SmFluid::set_acc(int acc) {
    m_acc = acc;
}

void SmFluid::set_method(Fluxes method) {
    m_nf = NumFlux::create(method);
}

Distributor SmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](AmrStorage::Item &parent, mesh::Children &children) {

        for (auto &child: children) {

            Vector3d dr = child.center - parent.center;

            // PState child_state = parent(U).get_state().arr() +
            //                      parent(U).d_dx.arr() * dr.x() +
            //                      parent(U).d_dy.arr() * dr.y() +
            //                      parent(U).d_dz.arr() * dr.z();
            // child_state.energy = m_eos.energy_rp(child_state.density, child_state.pressure);
            // child(U).set_state(child_state);

            PState child_state;
            child_state.density = parent(U).density +
                                    parent(U).d_dx.density * dr.x() +
                                    parent(U).d_dy.density * dr.y() +
                                    parent(U).d_dz.density * dr.z(); 
            
            child_state.velocity = parent(U).velocity +
                                   (parent(U).density / child_state.density) * (
                                        parent(U).d_dx.velocity * dr.x() +
                                        parent(U).d_dy.velocity * dr.y() +
                                        parent(U).d_dz.velocity * dr.z()
                                    );

            zephyr::phys::dRdE rhoe = m_eos.pressure_re(parent(U).density, parent(U).energy, {.deriv = true});

            child_state.energy = parent(U).energy -
                                    0.5 * (parent(U).velocity - child_state.velocity).squaredNorm() +
                                 (parent(U).density / child_state.density) * (
                                            (parent(U).d_dx.pressure - rhoe.dR * parent(U).d_dx.density) / rhoe.dE * dr.x() + 
                                            (parent(U).d_dy.pressure - rhoe.dR * parent(U).d_dy.density) / rhoe.dE * dr.y() +
                                            (parent(U).d_dz.pressure - rhoe.dR * parent(U).d_dz.density) / rhoe.dE * dr.z()
                                        );
            
            child_state.pressure = m_eos.pressure_re(child_state.density, child_state.energy);

            child(U).set_state(child_state);
        }

        // PState pc(0, {0,0,0}, 0, 0);
        // QState qc(pc);

        // for (auto &child: children) {
        //     QState qs(child(U).get_state());
        //     qc.arr() += qs.arr() * child.volume();
        // }

        // std::cout << std::setprecision(18) << QState(parent(U).get_state()).arr() * parent.volume()
        //           << "\n" << std::setprecision(18) << qc << "\n" << std::endl; 
    };

    distr.merge = [this](mesh::Children &children, AmrStorage::Item &parent) {
        PState sum;
        for (auto &child: children) {
            sum.arr() += child(U).get_state().arr() * child.volume();
        }
        parent(U).set_state(sum.arr() / parent.volume());
        parent(U).energy = m_eos.energy_rp(parent(U).density, parent(U).pressure);
    };

    return distr;
}

void SmFluid::set_flags(Mesh &mesh) {
    if (!mesh.is_adaptive()) {
        return;
    }

    compute_grad(mesh, get_current_sm);

    for (auto cell: mesh) {
        
        if (cell.flag() == 2)
            continue;

        double p = cell(U).pressure;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            PState t = face.neib(U).get_state() - cell(U).get_state();
            if (abs(t.density) > 0.2 * abs(cell(U).density) || // 0.1
                abs(t.pressure) > 0.5 * abs(cell(U).pressure))    // 0.3
                {
                    need_split = true;
                    break;
                }
        }
        if (need_split) {
            cell.set_flag(1);
        } else {
            cell.set_flag(-1);
        }
    }
}

};