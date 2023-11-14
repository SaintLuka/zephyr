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

void SmFluid::init_cells(Mesh &mesh, const phys::ClassicTest &test) {
    // Заполняем начальные данные
    for (auto cell: mesh) {
        cell(U).rho1 = test.density(cell.center());
        cell(U).v1 = test.velocity(cell.center());
        cell(U).p1 = test.pressure(cell.center());
        cell(U).e1 = m_eos.energy_rp(cell(U).rho1, cell(U).p1);
        cell(U).half = cell(U).get_state();
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

void SmFluid::set_accuracy(int acc) { 
    m_acc = acc;    // 1 или 2
};

void SmFluid::fluxes_stage1(Mesh& mesh) {
    for (auto cell : mesh) {
        // Примитивный вектор в ячейке
        PState zc = cell(U).get_state();
        // Консервативный вектор в ячейке
        QState qc(zc);
        //
        if (m_acc == 1) {
            cell(U).f1 = calc_flux(cell);
            qc.vec() -= m_dt / cell.volume() * cell(U).f1.vec();
        }
        if (m_acc == 2) {
            cell(U).f1 = calc_flux_extra(cell);
            qc.vec() -= 0.5 * m_dt / cell.volume() * cell(U).f1.vec();
        }
        PState Zc(qc, m_eos);
        cell(U).half = Zc;
    }    
}

void SmFluid::fluxes_stage2(Mesh& mesh) {
    if (m_acc == 1) {
        for (auto cell : mesh) {
            cell(U).set_state(cell(U).half);
        }
    }
    if (m_acc == 2) {
        for (auto cell : mesh) {
            // Примитивный вектор в ячейке
            PState zc = cell(U).get_state();
            // Консервативный вектор в ячейке
            QState qc(zc);
            // Расчет потока f2
            cell(U).f2 = calc_flux_extra(cell);
            //
            qc.vec() -= m_dt / cell.volume() * cell(U).f2.vec();
            //
            PState Zc(qc, m_eos);
            //
            cell(U).next = Zc;
        }
        for (auto cell : mesh) {
            cell(U).set_state(cell(U).next);
        }
    }
}

void SmFluid::compute_grad(Mesh& mesh) {
    for (auto cell : mesh) { 
        
        PState zx;
        PState zy;
        PState zz;

        PState zc = cell(U).half;

        for (auto &face: cell.faces()) {
            auto neib = face.neib();

            PState zn = neib(U).half;

            Vector3d S = 0.5 * face.normal() * face.area();  

            zx.vec() += (zc.vec() + zn.vec()) * S.x();
            zy.vec() += (zc.vec() + zn.vec()) * S.y();
            zz.vec() += (zc.vec() + zn.vec()) * S.z(); 
        }

        zx.vec() /= cell.volume();
        zy.vec() /= cell.volume();
        zz.vec() /= cell.volume();

        cell(U).d_dx = zx;
        cell(U).d_dy = zy;
        cell(U).d_dz = zz;
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
    for (auto& face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState zn(zc);

        //Единмтаенное различие
        if (face.is_boundary()) {
            // Граничные условия типа "стенка"
            // double vn = zc.velocity.dot(face.normal());
            // zn.velocity -= 2.0 * vn * face.normal();
        } 
        else {
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

/// берем только слой half
Flux SmFluid::calc_flux_extra(ICell &cell) { 
    // Примитивный вектор в ячейке
    PState zc = cell(U).half;
    
    // Переменная для потока
    Flux flux;
    for (auto& face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState zn(zc);

        if (face.is_boundary()) {
            // Граничные условия типа "стенка"
            // double vn = zc.velocity.dot(face.normal());
            // zn.velocity -= 2.0 * vn * face.normal();
        } 
        else {
            zn = face.neib()(U).half;
        }

        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();
        
        // Пока проблемки с периодичностью
        Vector3d neib_c = face.neib().center();

        // Второй порядок
        auto fe = FaceExtra::ATvL(
                    zc,        cell(U).d_dx,        cell(U).d_dy,        cell(U).d_dz,
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

void SmFluid::update(Mesh& mesh) {
    /// @brief Выбор шага интегрирования
    compute_dt(mesh);

    /// @brief Считает градиенты для переменных с основного слоя (только для 2ого порядка)
    if (m_acc == 2) {
        compute_grad(mesh);
    }
    /// @brief первая итерация
    fluxes_stage1(mesh);
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
}
