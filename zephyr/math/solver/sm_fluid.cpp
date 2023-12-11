#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/compute_grad.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace smf;
using namespace zephyr::phys;

static const SmFluid::State U = SmFluid::datatype();

[[nodiscard]] SmFluid::State SmFluid::datatype() {
    return {};
}

smf::PState get_smf_pstate(EuCell &cell) {
    return cell(U).get_state();
}

smf::PState get_smf_half(EuCell &cell) {
    return cell(U).half;
}

SmFluid::SmFluid(const phys::Eos &eos, Fluxes flux = Fluxes::HLLC2) : m_eos(eos) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.9;
    m_dt = std::numeric_limits<double>::max();
}

[[nodiscard]] double SmFluid::CFL() const {
    return m_CFL;
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

[[nodiscard]] double SmFluid::dt() const {
    return m_dt;
}

void SmFluid::init_cells(EuMesh &mesh, const phys::ClassicTest &test) {
    // Заполняем начальные данные
    for (auto cell: mesh) {
        cell(U).rho = test.density(cell.center());
        cell(U).v = test.velocity(cell.center());
        cell(U).p = test.pressure(cell.center());
        cell(U).e = m_eos.energy_rp(cell(U).rho, cell(U).p);
    }
}

[[nodiscard]] std::string SmFluid::get_flux_name() const {
    return m_nf->get_name();
}

Vector3d SmFluid::velocity(const Vector3d &c) const {
    return Vector3d::UnitX();
};

void SmFluid::compute_dt(EuMesh &mesh) {
    m_dt = std::numeric_limits<double>::max();
    for (auto cell: mesh) {
        // скорость звука
        double c = m_eos.sound_speed_rp(cell(U).rho, cell(U).p);
        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell(U).v.dot(face.normal());

            // Максимальное по модулю СЗ
            double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

            // Условие КФЛ
            m_dt = std::min(m_dt, cell.volume() / face.area() / lambda);
        }
    }
    m_dt *= m_CFL;
};

void SmFluid::set_accuracy(int acc) {
    m_acc = acc;    // 1 или 2
};

void SmFluid::fluxes_stage1(EuMesh &mesh) {
    for (auto cell: mesh) {
        // Консервативный вектор в ячейке
        QState qc(cell(U).get_state());
        if (m_acc == 1) {
            qc.vec() -= m_dt / cell.volume() * calc_flux(cell).vec();
        }
        if (m_acc == 2) {
            qc.vec() -= 0.5 * m_dt / cell.volume() * calc_flux_extra(cell, false).vec();
        }
        cell(U).half = PState(qc, m_eos);
    }
}

void SmFluid::fluxes_stage2(EuMesh &mesh) {
    if (m_acc == 1) {
        for (auto cell: mesh) 
            cell(U).set_state(cell(U).half);
    }
    if (m_acc == 2) {
        for (auto cell: mesh) {
            // Консервативный вектор в ячейке
            QState qc(cell(U).get_state());
            // Расчет потока f2
            qc.vec() -= m_dt / cell.volume() * calc_flux_extra(cell, true).vec();
            cell(U).next = PState(qc, m_eos);
        }
        for (auto cell: mesh) 
            cell(U).set_state(cell(U).next);
    }
}

void SmFluid::compute_grad(EuMesh &mesh,const std::function<smf::PState(zephyr::mesh::EuCell &)> &to_state) const {
    for (auto cell: mesh) {
        auto grad = math::compute_grad_eu<smf::PState>(cell, to_state);
        cell(U).d_dx = grad[0];
        cell(U).d_dy = grad[1];
        cell(U).d_dz = grad[2];
    }
}



// для первого порядка точности
Flux SmFluid::calc_flux(EuCell &cell) {
    // Примитивный вектор в ячейке
    PState zc = cell(U).get_state();

    // Консервативный вектор в ячейке
    QState qc(zc);

    // Переменная для потока
    Flux flux;
    for (auto &face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState zn(zc);

        //Единмтаенное различие
        if (face.is_boundary()) {
            // Граничные условия типа "стенка"
            // double vn = zc.velocity.dot(face.normal());
            // zn.velocity -= 2.0 * vn * face.normal();
        } else {
            zn = face.neib()(U).get_state();
        }

        // Значение на грани со стороны ячейки
        PState zm = zc.in_local(normal);

        // Значение на грани со стороны соседа
        PState zp = zn.in_local(normal);

        // Численный поток на грани
        auto loc_flux = m_nf->flux(zm, zp, m_eos);
        loc_flux.to_global(normal);

        // Суммируем поток
        flux.vec() += loc_flux.vec() * face.area();
    }

    return flux;
};

/// берем только слой half]
Flux SmFluid::calc_flux_extra(EuCell &cell, bool from_begin) {
    // Примитивный вектор в ячейке
    PState zc = cell(U).half;

    // Переменная для потока
    Flux flux;
    for (auto &face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState zn(zc);

        if (face.is_boundary()) {
            // Граничные условия типа "стенка"
            // double vn = zc.velocity.dot(face.normal());
            // zn.velocity -= 2.0 * vn * face.normal();
        } else {
            if(from_begin)
                zn = face.neib()(U).get_state();
            else
                zn = face.neib()(U).half;
        }

        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();
        Vector3d neib_c = 2 * face_c - cell_c;

        // Второй порядок
        auto fe = FaceExtra::ATvL(
                zc, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                zn, face.neib()(U).d_dx, face.neib()(U).d_dy, face.neib()(U).d_dz,
                cell_c, neib_c, face_c);

        PState zm = fe.m(zc).in_local(normal); // -
        PState zp = fe.p(zn).in_local(normal); // + 

        zm.energy = m_eos.energy_rp(zm.density, zm.pressure);
        zp.energy = m_eos.energy_rp(zp.density, zp.pressure);

        // Численный поток на грани
        auto loc_flux = m_nf->flux(zm, zp, m_eos);
        loc_flux.to_global(normal);

        // Суммируем поток
        flux.vec() += loc_flux.vec() * face.area();
    }

    return flux;
};

void SmFluid::update(EuMesh &mesh) {
    /// @brief Выбор шага интегрирования
    compute_dt(mesh);

    /// @brief Считает градиенты для переменных с основного слоя (только для 2ого порядка)
    if (m_acc == 2) {
        compute_grad(mesh, get_smf_pstate);
    }
    fluxes_stage1(mesh);
    if (m_acc == 2) {
        compute_grad(mesh, get_smf_half);
    }
    fluxes_stage2(mesh);
    /// 
    m_time += m_dt;
    m_step += 1;
};

Distributor SmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](AmrStorage::Item &parent, mesh::Children &children) {
        for (auto &child: children) {
            Vector3d dr = child.center - parent.center;
            PState child_state = parent(U).get_state().vec() +
                                 parent(U).d_dx.vec() * dr.x() +
                                 parent(U).d_dy.vec() * dr.y() +
                                 parent(U).d_dz.vec() * dr.z();
            child_state.energy = m_eos.energy_rp(child_state.density, child_state.pressure);
            child(U).set_state(child_state);
        }
    };

    distr.merge = [this](mesh::Children &children, AmrStorage::Item &parent) {
        PState sum;
        for (auto &child: children) {
            sum.vec() += child(U).get_state().vec() * child.volume();
        }
        parent(U).set_state(sum.vec() / parent.volume());
        parent(U).e = m_eos.energy_rp(parent(U).rho, parent(U).p);
    };

    return distr;
}

[[nodiscard]] double SmFluid::get_time() const {
    return m_time;
};

[[nodiscard]] double SmFluid::get_step() const {
    return m_step;
}
}
