#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

namespace zephyr { namespace math {

using mesh::Storage;
using namespace geom;
using namespace smf;

static const SmFluid::State U = SmFluid::datatype();

[[nodiscard]] SmFluid::State SmFluid::datatype() {
    return {};
};

SmFluid::SmFluid(IdealGas &eos) {
    m_nf = HLL::create();
    m_CFL = 0.5;
    m_dt = std::numeric_limits<double>::max(); 
    m_limiter = Limiter("minmod");
    m_eos = eos;
};

[[nodiscard]] double SmFluid::CFL() const {
    return m_CFL;
};

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
};

[[nodiscard]] double SmFluid::dt() const {
    return m_dt;
};

Vector3d SmFluid::velocity(const Vector3d& c) const {
    return Vector3d::UnitX();
};

double SmFluid::compute_dt(Mesh &mesh) {
    m_dt = std::numeric_limits<double>::max();
    for (auto cell: mesh) {
        // скорость звука
        double c = m_eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell(U).v1.dot(face.normal());

            // Максимальное по модулю СЗ
            double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

            // Условие КФЛ
            m_dt = std::min(m_dt, cell.volume() / face.area() / lambda);
        }
    }
    m_dt *= m_CFL;

    return m_dt;
};

// 1ый порядок
void SmFluid::fluxes(Mesh &mesh) {
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
            auto loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        // Новое значение в ячейке (консервативные переменные)
        QState Qc = qc.vec() - m_dt * flux.vec() / cell.volume();

        // Новое значение примитивных переменных
        PState Zc(Qc, m_eos);

        cell(U).set_state2(Zc);
    }
};

// Рунге-Кутта 2ого порядка
void SmFluid::fluxes2(Mesh &mesh, int stage) {
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
            auto loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        if (stage==1) {
            // Новое значение в ячейке (консервативные переменные)
            QState Qc = qc.vec() - 0.5 * m_dt * flux.vec() / cell.volume();

            // Новое значение примитивных переменных
            PState Zc(Qc, m_eos);

            cell(U).set_state2(Zc);
        }

        if (stage == 2) {
            // Новое значение в ячейке (консервативные переменные)
            QState Qc = qc.vec() - m_dt * flux.vec() / cell.volume();

            // Новое значение примитивных переменных
            PState Zc(Qc, m_eos);

            cell(U).set_state2(Zc);
        }
    }    

    if (stage == 1) {
        fluxes2(mesh, 2);
    }
}

void SmFluid::update(Mesh& mesh) {
    m_dt = compute_dt(mesh);
    fluxes2(mesh, 1);
    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).swap();
    }
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
}

