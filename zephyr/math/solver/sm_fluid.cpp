#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>


namespace zephyr { namespace math {

using mesh::Storage;
using namespace geom;
using namespace smf;
using namespace zephyr::phys;

static const SmFluid::State U = SmFluid::datatype();

[[nodiscard]] SmFluid::State SmFluid::datatype() {
    return {};
}

SmFluid::SmFluid(const phys::Eos &eos, Fluxes flux = Fluxes::HLLC2) : m_eos(eos) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.5;
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

void SmFluid::init_cells(Mesh &mesh, const phys::ClassicTest &test) {
    // Заполняем начальные данные
    for (auto cell: mesh) {
        cell(U).rho1 = test.density(cell.center());
        cell(U).v1 = test.velocity(cell.center());
        cell(U).p1 = test.pressure(cell.center());
        cell(U).e1 = m_eos.energy_rp(cell(U).rho1, cell(U).p1);
    }
}

[[nodiscard]] std::string SmFluid::get_flux_name() const {
    return m_nf->get_name();
}
 


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
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P1
        Flux P1 = calc_flux(cell, 0);
        // Новое значение в ячейке (консервативные переменные)
        QState Qc = qc.vec() - m_dt * P1.vec() / cell.volume();
        // Новое значение примитивных переменных
        PState Zc(Qc, m_eos);
        //
        cell(U).set_state(Zc,2);
    }
};

// Рунге-Кутта 2ого порядка
void SmFluid::fluxes2(Mesh &mesh) {
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P1
        Flux P1 = calc_flux(cell, 0);
        //
        QState Qch = qc.vec() - 0.5 * m_dt / cell.volume() * P1.vec();
        //
        PState Zch(Qch, m_eos);
        //
        cell(U).set_state(Zch, 1);
    }
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P1
        Flux P2 = calc_flux(cell, 1);
        //
        QState Qc = qc.vec() - m_dt / cell.volume() * P2.vec();
        //
        PState Zc(Qc, m_eos);
        //
        cell(U).set_state(Zc, 2);
    }
}

//leapfrog
void SmFluid::fluxes3(Mesh &mesh) {
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P1
        cell(U).P1 = calc_flux(cell, 0);
        //
        QState Qc1 = qc0.vec() - m_dt / cell.volume() * cell(U).P1.vec();
        //
        PState Zc1(Qc1, m_eos);
        //
        cell(U).set_state(Zc1, 1);
    }
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P1
        cell(U).P2 = calc_flux(cell, 1);
        //
        QState Qc2 = qc0.vec() - 0.5 * m_dt / cell.volume() * cell(U).P1.vec() - m_dt / cell.volume() * cell(U).P2.vec();
        //
        PState Zc2(Qc2, m_eos);
        //
        cell(U).set_state(Zc2, 2);
    }    
}

// Рунге-Кутта 4ого порядка
void SmFluid::fluxes4(Mesh &mesh) {
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P1
        cell(U).P1 = calc_flux(cell, 0);
        //
        QState Qc1 = qc0.vec() - m_dt / cell.volume() * cell(U).P1.vec();
        //
        PState Zc1(Qc1, m_eos);
        //
        cell(U).set_state(Zc1, 1);
    }
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P2
        cell(U).P2 = calc_flux(cell, 1);
        //
        QState Qc2 = qc0.vec() - 0.5 * m_dt / cell.volume() * cell(U).P1.vec() - m_dt / cell.volume() * cell(U).P2.vec();
        //
        PState Zc2(Qc2, m_eos);
        //
        cell(U).set_state(Zc2, 2);
    } 
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc2 = cell(U).get_state(2);
        // Консервативный вектор в ячейке
        QState qc2(zc2);
        // Примитивный вектор в ячейке
        PState zc1 = cell(U).get_state(1);
        // Консервативный вектор в ячейке
        QState qc1(zc1);
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P3
        cell(U).P3 = calc_flux(cell, 2);
        //
        QState Qc3 = 2 * qc2.vec() - 2 * qc1.vec() + qc0.vec() - \
                        2 * m_dt / cell.volume() * cell(U).P2.vec() + \
                        0.5 * m_dt / cell.volume() * cell(U).P1.vec();
        //
        PState Zc3(Qc3, m_eos);
        //
        cell(U).set_state(Zc3, 3);
    }
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P3
        cell(U).P4 = calc_flux(cell, 3);
        //
        QState Qc = qc.vec() - m_dt / 6 / cell.volume() * \
                    (0.5 * cell(U).P1.vec() + 2 * cell(U).P2.vec() + cell(U).P3.vec() + cell(U).P4.vec());
        //
        PState Zc(Qc, m_eos);
        //
        cell(U).set_state(Zc, 2);
    }   
}

Flux SmFluid::calc_flux(ICell &cell, int stage) {
    // Примитивный вектор в ячейке
    PState zc = cell(U).get_state(stage);

    // Консервативный вектор в ячейке
    QState qc(zc);
    
    // Переменная для потока
    Flux flux;
    for (auto& face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState zn(zc);

        //Единмтаенное различие
        if (face.is_boundary()) {
            // Граничные условия типа "стенка"
            double vn = zc.velocity.dot(face.normal());
            zn.velocity -= 2.0 * vn * face.normal();
        } 
        else {
            zn = face.neib()(U).get_state(stage);
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


void SmFluid::update(Mesh& mesh) {
    for (auto cell : mesh) {
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
