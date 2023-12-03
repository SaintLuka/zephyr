#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/compute_grad.h>

namespace zephyr::math {

using mesh::Storage;
using namespace geom;
using namespace smf;
using namespace zephyr::phys;

static const SmFluid::State U = SmFluid::datatype();

[[nodiscard]] SmFluid::State SmFluid::datatype() {
    return {};
}

smf::PState get_smf_pstate(ICell &cell) {
    return cell(U).get_state();
}

smf::PState get_smf_half(ICell &cell) {
    return cell(U).half;
}

SmFluid::SmFluid(const phys::Eos &eos, Fluxes flux = Fluxes::HLLC2) : m_eos(eos) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.9;
    m_dt = std::numeric_limits<double>::max();
}

double SmFluid::compute_dt(EuCell &cell) {
    double max_area = 0.0;
    for (auto &face: cell.faces()) {
        max_area = std::max(max_area, face.area());
    }
    double dx = cell.volume() / max_area;
    return m_CFL * dx / velocity(cell.center()).norm();
}

void SmFluid::compute_grad(EuCell &cell, int stage) {

}

void SmFluid::fluxes(EuCell &cell, int stage) {

};

void SmFluid::update(EuMesh& mesh, IdealGas &eos) {
    NumFlux::Ptr nf = HLL::create(); // ?

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

void SmFluid::fluxes_stage2(Mesh &mesh) {
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

void SmFluid::compute_grad(Mesh &mesh,const std::function<smf::PState(zephyr::mesh::ICell &)> &to_state) const {
    for (auto cell: mesh) {
        auto grad = math::compute_grad<smf::PState>(cell, to_state);
        cell(U).d_dx = grad[0];
        cell(U).d_dy = grad[1];
        cell(U).d_dz = grad[2];
    }
}

// для первого порядка точности
Flux SmFluid::calc_flux(ICell &cell) {
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
Flux SmFluid::calc_flux_extra(ICell &cell, bool from_begin) {
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

void SmFluid::update(Mesh &mesh) {
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

[[nodiscard]] double SmFluid::get_time() const {
    return m_time;
};

[[nodiscard]] double SmFluid::get_step() const {
    return m_step;
}
}

