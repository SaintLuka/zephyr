#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/geom/face.h>


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
    m_limiter = Limiter("minmod");
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

void SmFluid::compute_dt(Mesh &mesh) {
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
};

// 1ый порядок
void SmFluid::fluxes(Mesh &mesh) {
    for (auto cell: mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P1
        cell(U).P1 = calc_flux(cell, 0);
        // Новое значение в ячейке (консервативные переменные)
        QState Qc = qc.vec() - m_dt * cell(U).P1.vec() / cell.volume();
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
        cell(U).P1 = calc_flux_extra(cell, 0);
        //
        QState Qc1 = qc.vec() - 0.5 * m_dt / cell.volume() * cell(U).P1.vec();
        //
        PState Zc1(Qc1, m_eos);
        //
        cell(U).set_state(Zc1, 1);
    }
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc(zc);
        // Расчет потока P2
        cell(U).P2 = calc_flux_extra(cell, 1);
        //
        QState Qc = qc.vec() - m_dt / cell.volume() * cell(U).P2.vec();
        //
        PState Zc(Qc, m_eos);
        //
        cell(U).set_state(Zc, 2);
    }
}

// leapfrog
void SmFluid::fluxes3(Mesh &mesh) {
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc0 = cell(U).get_state(0);
        // Консервативный вектор в ячейке
        QState qc0(zc0);
        // Расчет потока P1
        cell(U).P1 = calc_flux_extra(cell, 0);
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
        cell(U).P2 = calc_flux_extra(cell, 1);
        //
        QState Qc2 = qc0.vec() - 0.5 * m_dt / cell.volume() * (cell(U).P1.vec() + cell(U).P2.vec());
        //
        PState Zc2(Qc2, m_eos);
        //
        cell(U).set_state(Zc2, 2);
    }    
}

void SmFluid::compute_grad(ICell &cell, int stage) {
    PState zx;
    PState zy;
    PState zz;

    PState zc = cell(U).get_state(stage);

    for (auto &face: cell.faces()) {
        auto neib = face.neib();

        PState zn = neib(U).get_state(stage);

        Vector3d S = 0.5 * face.normal() * face.area();  

        zx.vec() += (zc.vec() + zn.vec()) * S.x();
        zy.vec() += (zc.vec() + zn.vec()) * S.y();
        zz.vec() += (zc.vec() + zn.vec()) * S.z(); 
    }

    zx.vec() /= cell.volume();
    zy.vec() /= cell.volume();
    zy.vec() /= cell.volume();

    cell(U).set_zx(zx);
    cell(U).set_zy(zy);
    cell(U).set_zz(zz);
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

Flux SmFluid::calc_flux_extra(ICell &cell, int stage) {
    // Примитивный вектор в ячейке
    PState zc = cell(U).get_state(stage);
    
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
            zn = face.neib()(U).get_state(stage);
        }


        // только для 0 stage - а
        compute_grad(cell, stage);
        compute_grad(*face.neib(), stage);

        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();
        
        // Пока проблемки с периодичностью
        Vector3d neib_c = 2.0 * face_c - cell_c;

        // Второй порядок
        auto fe = FaceExtra::ATvL(
                    zc,        cell(U).zx(),        cell(U).zy(),        cell(U).zz(),
                    zn, face.neib()(U).zx(), face.neib()(U).zy(), face.neib()(U).zz(),
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

void SmFluid::update(Mesh& mesh, int flux_type) {
    compute_dt(mesh);
    
    switch(flux_type) {
        case 1:
            fluxes(mesh);
            break;
        case 2:
            fluxes2(mesh);
            break;
        case 3:
            fluxes3(mesh);
            break;
    }
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
