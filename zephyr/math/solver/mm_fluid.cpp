#include <zephyr/math/solver/mm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/compute_grad.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace mmf;
using zephyr::utils::threads;

static const MmFluid::State U = MmFluid::datatype();

[[nodiscard]] MmFluid::State MmFluid::datatype() {
    return {};
}

/// Используется при дебаге (в калькуляторе можно набрать to_state(cell) и посмотреть что лежит внутри
MmFluid::State to_state(Cell &cell) {
    return cell(U);
}

mmf::PState get_current(Cell &cell) {
    return cell(U).get_pstate();
}

mmf::PState get_half(Cell &cell) {
    return cell(U).half;
}

MmFluid::MmFluid(const phys::Materials &mixture, Fluxes flux, double g) : mixture(mixture), g(g) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.5;
    m_dt = std::numeric_limits<double>::max();
}

void MmFluid::update(Mesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        fluxes(mesh);
    } else {
        compute_grad(mesh, get_current);
        fluxes_stage1(mesh);
        compute_grad(mesh, get_half);
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

double MmFluid::compute_dt(Mesh &mesh) {
    m_dt = std::numeric_limits<double>::max();
    for (auto cell: mesh) {
        // скорость звука
        double c = mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac);
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

    return m_dt;
}

void MmFluid::fluxes(Mesh &mesh) const {
    mesh.for_each([this](Cell &cell) {
        // Примитивный вектор в ячейке
        PState p_self = cell(U).get_pstate();

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
                p_neib = face.neib()(U).get_pstate();
            } else if (face.flag() == Boundary::WALL) {
                Vector3d Vn = normal * p_self.velocity.dot(normal);
                p_neib.velocity = p_self.velocity - 2 * Vn; // Vt - Vn = p_self.velocity - Vn - Vn
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
        q_self2.momentum.y() -= m_dt * g * p_self.density;

        // Новое значение примитивных переменных
        PState p_self2(q_self2, mixture, p_self.pressure, p_self.energy);

//        std::cout << "step: " << m_step << '\n';
//        std::cout << "cell_idx: " << cell.b_idx() << '\n';
//        std::cout << "p_self: " << p_self.mass_frac << '\n';
//        std::cout << "flux: " << flux.mass_frac << '\n';
//        std::cout << "q_self2: " << q_self2.mass_frac << '\n';
//        std::cout << "p_self2: " << p_self2.mass_frac << "\n\n";

        cell(U).next = p_self2;

        if (cell(U).is_bad()) {
            std::cerr << "calc new cell with idx: " << cell.b_idx() << " state from step " << m_step << " to " << m_step + 1
                      << "\n";
            std::cerr << "x: " << cell.center().x() << ", y: " << cell.center().y() << "\n";
            std::cerr << cell(U);
            throw std::runtime_error("bad cell");
        }
    });
}

void MmFluid::compute_grad(Mesh &mesh, const std::function<mmf::PState(Cell &)> &to_state) const {
    mesh.for_each([&to_state](Cell &cell) -> void {
        auto grad = math::compute_grad<mmf::PState>(cell, to_state);
        cell(U).d_dx = grad[0];
        cell(U).d_dy = grad[1];
        cell(U).d_dz = grad[2];
    });
}

void MmFluid::fluxes_stage1(Mesh &mesh) const {
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= 0.5 * m_dt / cell.volume() * calc_flux_extra(cell, true).vec();
        qc.momentum.y() -= m_dt * g * cell(U).rho;
        cell(U).half = PState(qc, mixture);
    });
}

void MmFluid::fluxes_stage2(Mesh &mesh) const {
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= m_dt / cell.volume() * calc_flux_extra(cell, false).vec();
        qc.momentum.y() -= m_dt * g * cell(U).half.density;
        cell(U).next = PState(qc, mixture);
    });
}

mmf::Flux MmFluid::calc_flux_extra(Cell &cell, bool from_begin) const {
    // Примитивный вектор в ячейке
    PState p_self = cell(U).half;
    if (from_begin)
        p_self = cell(U).get_pstate();

    // Переменная для потока
    Flux flux;
    for (auto &face: cell.faces()) {
        // Внешняя нормаль
        auto &normal = face.normal();

        // Примитивный вектор соседа
        PState p_neib = p_self;
        if (!face.is_boundary()) {
            p_neib = face.neib()(U).half;
            if (from_begin)
                p_neib = face.neib()(U).get_pstate();
        } else if (face.flag() == Boundary::WALL) {
            Vector3d Vn = normal * p_self.velocity.dot(normal);
            p_neib.velocity = p_self.velocity - 2 * Vn; // Vt - Vn = p_self.velocity - Vn - Vn
        }

        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();
        Vector3d neib_c = 2 * face_c - cell_c;

        auto face_extra = FaceExtra::ATvL(
                p_self, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                p_neib, face.neib()(U).d_dx, face.neib()(U).d_dy, face.neib()(U).d_dz,
                cell_c, neib_c, face_c);

        // рассчитываем расстояние на грани слева и справа
        PState p_minus = face_extra.m(p_self).in_local(normal);
        PState p_plus = face_extra.p(p_neib).in_local(normal);

        // исправляем возможные нефизичные значения массовых долей
        p_minus.mass_frac.fix();
        p_plus.mass_frac.fix();

        // пересчитываем энергию и температуру
        p_minus.energy = mixture.energy_rp(p_minus.density, p_minus.pressure, p_minus.mass_frac);
        p_plus.energy = mixture.energy_rp(p_plus.density, p_plus.pressure, p_plus.mass_frac);

        p_minus.temperature = mixture.temperature_rp(p_minus.density, p_minus.pressure, p_minus.mass_frac);
        p_plus.temperature = mixture.temperature_rp(p_plus.density, p_plus.pressure, p_plus.mass_frac);

        // Численный поток на грани
        auto loc_flux = m_nf->mm_flux(p_minus, p_plus, mixture);
        loc_flux.to_global(normal);

        // Суммируем поток
        flux.vec() += loc_flux.vec() * face.area();
    }

    return flux;
}

void MmFluid::swap(Mesh &mesh) {
    size_t step = m_step;
    mesh.for_each([step](Cell &cell) -> void {
        if (cell(U).is_bad()) {
            std::cerr << "update cell " << cell.b_idx() << " from step " << step << " to " << step + 1 << "\n";
            std::cerr << cell(U);
            throw std::runtime_error("bad cell");
        }
        cell(U).set_state(cell(U).next);
        cell(U).half = cell(U).next;
    });
    m_time += m_dt;
    m_step += 1;
}

void MmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

[[nodiscard]] double MmFluid::dt() const {
    return m_dt;
}

[[nodiscard]] double MmFluid::CFL() const {
    return m_CFL;
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

void MmFluid::set_acc(int acc) {
    m_acc = acc;
}

Distributor MmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](AmrStorage::Item &parent, mesh::Children &children) {
        for (auto &child: children) {
            Vector3d dr = child.center - parent.center;
            PState child_state = parent(U).get_pstate().vec() +
                                 parent(U).d_dx.vec() * dr.x() +
                                 parent(U).d_dy.vec() * dr.y() +
                                 parent(U).d_dz.vec() * dr.z();
            child_state.mass_frac = parent(U).mass_frac;
            child_state.energy = mixture.energy_rp(child_state.density, child_state.pressure, child_state.mass_frac);
            child_state.temperature = mixture.temperature_rp(child_state.density, child_state.pressure, child_state.mass_frac);
            child(U).set_state(child_state);
        }
    };

    distr.merge = [this](mesh::Children &children, AmrStorage::Item &parent) {
        PState sum;
        for (auto &child: children) {
            sum.vec() += child(U).get_pstate().vec() * child.volume();
        }
        parent(U).set_state(sum.vec() / parent.volume());
        parent(U).mass_frac.fix();
        parent(U).e = mixture.energy_rp(parent(U).rho, parent(U).p, parent(U).mass_frac);
        parent(U).t = mixture.temperature_rp(parent(U).rho, parent(U).p, parent(U).mass_frac);
    };

    return distr;
}

void MmFluid::set_flags(Mesh &mesh) {
    compute_grad(mesh, get_current);

    for (auto cell: mesh) {
        double p = cell(U).p;
        Fractions mass_frac = cell(U).mass_frac;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            if (abs(face.neib()(U).p - p) > 0.1 * p) {
                need_split = true;
                break;
            }

            // проверяем большое различие в долях веществ
            Fractions neib_mass_frac = face.neib()(U).mass_frac;
            for (int i = 0; i < Fractions::max_size; i++) {
                if (abs(mass_frac[i] - neib_mass_frac[i]) > 0.1) {
                    need_split = true;
                    break;
                }
            }
            if (need_split)
                break;
        }
        if (need_split) {
            cell.set_flag(1);
        } else {
            cell.set_flag(-1);
        }
    }
}

}