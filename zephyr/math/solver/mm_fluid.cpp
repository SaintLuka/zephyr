#include <zephyr/math/solver/mm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

using mesh::Storage;
using namespace geom;
using namespace mmf;

static const MmFluid::State U = MmFluid::datatype();

[[nodiscard]] MmFluid::State MmFluid::datatype() {
    return {};
}

MmFluid::State to_state(const zephyr::math::mmf::PState &state) {
    return MmFluid::State{state.density, state.density,
                          state.velocity, state.velocity,
                          state.pressure, state.pressure,
                          state.energy, state.energy,
                          state.temperature, state.temperature,
                          state.mass_frac, state.mass_frac};
}

MmFluid::State to_state(ICell &cell) {
    return cell(U);
}


MmFluid::MmFluid(const phys::Materials &mixture, Fluxes flux = Fluxes::GODUNOV) : mixture(mixture) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.5;
    m_dt = std::numeric_limits<double>::max();
}

[[nodiscard]] double MmFluid::CFL() const {
    return m_CFL;
}

void MmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

[[nodiscard]] double MmFluid::dt() const {
    return m_dt;
}

double MmFluid::compute_dt(Mesh &mesh) {
    m_dt = std::numeric_limits<double>::max();
    for (auto cell: mesh) {
        // скорость звука
        double c = mixture.sound_speed_rp(cell(U).rho1, cell(U).p1, cell(U).mass_frac1);
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
}

void MmFluid::fluxes(Mesh &mesh) {
    // Расчет по некоторой схеме
    for (auto cell: mesh) {
        // Примитивный вектор в ячейке
        PState p_self(cell(U).rho1, cell(U).v1, cell(U).p1, cell(U).e1, cell(U).t1, cell(U).mass_frac1);

        // Консервативный вектор в ячейке
        QState q_self(p_self);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto &normal = face.normal();

            // Примитивный вектор соседа
            PState p_neib(p_self);

            if (!face.is_boundary()) {
                auto neib = face.neib();
                p_neib.density = neib(U).rho1;
                p_neib.velocity = neib(U).v1;
                p_neib.pressure = neib(U).p1;
                p_neib.energy = neib(U).e1;
                p_neib.temperature = neib(U).t1;
                p_neib.mass_frac = neib(U).mass_frac1;
            }

            // Значение на грани со стороны ячейки
            PState zm = p_self.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = p_neib.in_local(normal);

            // Численный поток на грани
            auto loc_flux = m_nf->mm_flux(zm, zp, mixture);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        // Новое значение в ячейке (консервативные переменные)
        QState q_self2 = q_self.vec() - m_dt * flux.vec() / cell.volume();

        // Новое значение примитивных переменных
        PState p_self2(q_self2, mixture, p_self.pressure, p_self.energy);

//        std::cout << "step: " << m_step << '\n';
//        std::cout << "cell_idx: " << cell.b_idx() << '\n';
//        std::cout << "p_self: " << p_self.mass_frac << '\n';
//        std::cout << "flux: " << flux.mass_frac << '\n';
//        std::cout << "q_self2: " << q_self2.mass_frac << '\n';
//        std::cout << "p_self2: " << p_self2.mass_frac << "\n\n";

        cell(U).rho2 = p_self2.density;
        cell(U).v2 = p_self2.velocity;
        cell(U).p2 = p_self2.pressure;
        cell(U).e2 = p_self2.energy;
        cell(U).t2 = p_self2.temperature;
        cell(U).mass_frac2 = p_self2.mass_frac;

        if (cell(U).is_bad()) {
            std::cerr << "calc new cell " << cell.b_idx() << " state from step " << m_step << " to " << m_step + 1
                      << "\n";
            std::cerr << cell(U);
            throw std::runtime_error("bad cell");
        }
    }
}

void MmFluid::update(Mesh &mesh) {
    for (auto cell: mesh) {
        if (cell(U).is_bad()) {
            std::cerr << "update cell " << cell.b_idx() << " from step " << m_step << " to " << m_step + 1 << "\n";
            std::cerr << cell(U);
            throw std::runtime_error("bad cell");
        }
        std::swap(cell(U).rho1, cell(U).rho2);
        std::swap(cell(U).v1, cell(U).v2);
        std::swap(cell(U).p1, cell(U).p2);
        std::swap(cell(U).e1, cell(U).e2);
        std::swap(cell(U).t1, cell(U).t2);
        std::swap(cell(U).mass_frac1, cell(U).mass_frac2);
    }
    m_time += m_dt;
    m_step += 1;
}

void MmFluid::set_num_flux(Fluxes flux) {
    m_nf = NumFlux::create(flux);
}

[[nodiscard]] double MmFluid::get_time() const {
    return m_time;
}

[[nodiscard]] size_t MmFluid::get_step() const {
    return m_step;
}

[[nodiscard]] std::string MmFluid::get_flux_name() const {
    return m_nf->get_name();
}

}