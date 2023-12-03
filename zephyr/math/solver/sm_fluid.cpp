#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

namespace zephyr { namespace math {

using namespace zephyr::phys;
using namespace smf;

static const SmFluid::State U = SmFluid::datatype();

SmFluid::State SmFluid::datatype() {
    return {};
}

SmFluid::SmFluid() {
    
};

Vector3d SmFluid::velocity(const Vector3d& c) const {
    return Vector3d::UnitX();
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
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state1();

        // Консервативный вектор в ячейке
        QState qc(zc);

        // Переменная для потока
        Flux flux;
        for (auto& face: cell.faces()) {
            // Внешняя нормаль
            auto &normal = face.normal();

            // Примитивный вектор соседа
            PState zn(zc);

            if (face.is_boundary()) {
                const double d = 0.5;

                // Граничные условия типа "стенка"
                double vn = zc.velocity.dot(face.normal());
                zn.velocity -= 2.0 * vn * face.normal();
            }
            else {
                zn = face.neib()(U).get_state1();
            }

            // Значение на грани со стороны ячейки
            PState zm = zc.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = zn.in_local(normal);

            // Численный поток на грани
            auto loc_flux = nf->flux(zm, zp, eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        // Новое значение в ячейке (консервативные переменные)
        QState Qc = qc.vec() - compute_dt(cell) * flux.vec() / cell.volume();

        // Новое значение примитивных переменных
        PState Zc(Qc, eos);

        cell(U).set_state2(Zc);
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).swap();
    }
};

void SmFluid::solution_step() {

};

double SmFluid::dt() const {
    return m_dt;
};
}
}
