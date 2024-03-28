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
    m_CFL = 0.4;
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

    // Рассчитываем плотности компонент и скорости звука
    compute_components_chars(mesh);
}

void MmFluid::compute_components_chars(Mesh &mesh) {
    mesh.for_each([this](Cell &cell) {
        for (int i = 0; i < mixture.size(); i++) {
            if (!cell(U).mass_frac.has(i)) {
                cell(U).densities[i] = 0;
                cell(U).speeds[i] = 0;
                continue;
            }
            cell(U).densities[i] = 1.0 / mixture[i].volume_pt(cell(U).p, cell(U).t);;
            cell(U).speeds[i] = mixture[i].sound_speed_rp(cell(U).densities[i], cell(U).p);
        }
    });
}

double MmFluid::compute_dt(Mesh &mesh) {
    m_dt = mesh.min([this](Cell &cell) -> double {
        double min_dt = std::numeric_limits<double>::max();
        // скорость звука
        double c = mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac, {.T0=cell(U).t});
        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell(U).v.dot(face.normal());

            // Максимальное по модулю СЗ
            double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

            // Условие КФЛ
            min_dt = std::min(min_dt, cell.volume() / face.area() / lambda);
        }

        return min_dt;
    });
    m_dt *= m_CFL;

    return m_dt;
}

void MmFluid::check_state(const PState &next, const QState &qc, const PState &old, const std::string &func_name, size_t cell_idx) {
    if (!next.is_bad())
        return;
    out_mu.lock();
    std::cerr << "Bad cell: " << cell_idx << "\n";
    std::cerr << "Failed to calc PState from QState in " + func_name + "\n";
    std::cerr << "QState: " << qc << "\n";
    std::cerr << "PState: " << next << "\n";
    std::cerr << "Previous PState: " << old << "\n";
    out_mu.unlock();
    exit(1); // чтобы расчёт моментально выключился
    throw std::runtime_error("bad cell");
}

void MmFluid::fluxes(Mesh &mesh) {
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
                p_neib = face.neib(U).get_pstate();
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
            if (loc_flux.is_bad()) {
                std::cerr << "Step: " << m_step << '\n';
                std::cerr << "Bad flux for " << cell.b_idx() << " and " << face.neib().b_idx() << '\n';
                exit(1);
            }
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        // Новое значение в ячейке (консервативные переменные)
        QState q_self2 = q_self.vec() - m_dt * flux.vec() / cell.volume();
        q_self2.momentum.y() -= m_dt * g * p_self.density;

        // Новое значение примитивных переменных
        PState p_self2(q_self2, mixture, p_self.pressure, p_self.temperature);
        check_state(p_self2, q_self2, p_self, "fluxes", cell.b_idx());

//        std::cout << "step: " << m_step << '\n';
//        std::cout << "cell_idx: " << cell.b_idx() << '\n';
//        std::cout << "p_self: " << p_self.mass_frac << '\n';
//        std::cout << "flux: " << flux.mass_frac << '\n';
//        std::cout << "q_self2: " << q_self2.mass_frac << '\n';
//        std::cout << "p_self2: " << p_self2.mass_frac << "\n\n";

        cell(U).next = p_self2;

        if (cell(U).is_bad()) {
            out_mu.lock();
            std::cerr << "calc new cell with idx: " << cell.b_idx() << " state from step " << m_step << " to " << m_step + 1
                      << "\n";
            std::cerr << "x: " << cell.center().x() << ", y: " << cell.center().y() << "\n";
            std::cerr << cell(U);
            out_mu.unlock();
            throw std::runtime_error("bad cell");
        }
    });
}

PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
    if (flag != Boundary::WALL)
        return zc;

    PState zn(zc);
    Vector3d Vn = normal * zc.velocity.dot(normal);
    zn.velocity = zc.velocity - 2 * Vn;

    return zn;
}

void MmFluid::compute_grad(Mesh &mesh, const std::function<mmf::PState(Cell &)> &to_state) {
    mesh.for_each([&to_state](Cell &cell) -> void {
        auto grad = math::compute_gradient_LSM<mmf::PState>(cell, to_state, boundary_value);
        auto grad_lim = math::gradient_limiting<mmf::PState>(cell, grad, to_state, boundary_value);
        cell(U).d_dx = grad_lim[0];
        cell(U).d_dy = grad_lim[1];
        cell(U).d_dz = grad_lim[2];
    });
}

void MmFluid::fluxes_stage1(Mesh &mesh) {
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= 0.5 * m_dt / cell.volume() * calc_flux_extra(cell, true).vec();
        qc.momentum.y() -= m_dt * g * cell(U).rho;
        cell(U).half = PState(qc, mixture, cell(U).p, cell(U).t);
        check_state(cell(U).half, qc, cell(U).get_pstate(), "fluxes_stage1", cell.b_idx());
    });
}

void MmFluid::fluxes_stage2(Mesh &mesh) {
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= m_dt / cell.volume() * calc_flux_extra(cell, false).vec();
        qc.momentum.y() -= m_dt * g * cell(U).half.density;
        cell(U).next = PState(qc, mixture, cell(U).half.pressure, cell(U).half.temperature);
        check_state(cell(U).next, qc, cell(U).half, "fluxes_stage2", cell.b_idx());
    });
}

mmf::Flux MmFluid::calc_flux_extra(Cell &cell, bool from_begin) {
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

        /*
        Vector3d dr = face.center() - cell.center();
        PState p_minus = p_self.vec() +
                         cell(U).d_dx.vec() * dr.x() +
                         cell(U).d_dy.vec() * dr.y() +
                         cell(U).d_dz.vec() * dr.z();
        p_minus.to_local(normal);

        dr = face.center() - neib_c;
        PState p_plus = p_minus;
        if (!face.is_boundary())
            p_plus = p_neib.vec() +
                     face.neib()(U).d_dx.vec() * dr.x() +
                     face.neib()(U).d_dy.vec() * dr.y() +
                     face.neib()(U).d_dz.vec() * dr.z();
        else if (face.flag() == Boundary::WALL) {
            Vector3d Vn = normal * p_minus.velocity.dot(normal);
            p_plus.velocity = p_minus.velocity - 2 * Vn;
        }
        p_plus.to_local(normal);
        */

        // исправляем возможные нефизичные значения массовых долей
        p_minus.mass_frac.fix();
        p_plus.mass_frac.fix();

        if (face.flag() == Boundary::WALL) {
            Vector3d Vn = normal * p_minus.velocity.dot(normal);
            p_plus.velocity = p_minus.velocity - 2 * Vn; // Vt - Vn = p_self.velocity - Vn - Vn
        }

        // пересчитываем энергию и температуру
//        p_minus.sync_temperature_energy_rp(mixture);
//        p_plus.sync_temperature_energy_rp(mixture);

        auto is_bad_state = [](const PState &state) -> bool {
            return std::isinf(state.density) || std::isnan(state.density) ||
                   std::isinf(state.velocity.x()) || std::isnan(state.velocity.x()) ||
                   std::isinf(state.velocity.y()) || std::isnan(state.velocity.y()) ||
                   std::isinf(state.velocity.z()) || std::isnan(state.velocity.z()) ||
                   std::isinf(state.pressure) || std::isnan(state.pressure) ||
                   state.mass_frac.empty();
        };

        if (is_bad_state(p_minus) || is_bad_state(p_plus)) {
            out_mu.lock();
            std::cerr << "Can't compute PState on face\n";
            std::cerr << "p_minus: " << p_minus << '\n' << "p_plus: " << p_plus << '\n';
            std::cerr << "p_self: " << p_self << '\n' << "p_neib: " << p_neib << '\n';
            if (is_bad_state(p_minus)) {
                std::cerr << "dx: " << cell(U).d_dx << '\n' << "dy: " << cell(U).d_dy << '\n' << "dz: " << cell(U).d_dz << '\n';
                std::cerr << "dr: {" << (face_c - cell_c).x() << ", " << (face_c - cell_c).y() << ", " << (face_c - cell_c).z() << "}\n";
            } else {
                std::cerr << "dx: " << face.neib()(U).d_dx << '\n' << "dy: " << face.neib()(U).d_dy << '\n' << "dz: " << face.neib()(U).d_dz << '\n';
                std::cerr << "dr: {" << (face_c - neib_c).x() << ", " << (face_c - neib_c).y() << ", " << (face_c - neib_c).z() << "}\n";
            }
            out_mu.unlock();
            exit(1);
        }

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
//            child_state.mass_frac = parent(U).mass_frac;
            child_state.mass_frac.fix();
            child_state.sync_temperature_energy_rp(mixture, {.T0 = parent(U).t});
            child(U).set_state(child_state);
            if (m_step > 0 && child(U).is_bad1()) {
                std::cerr << "Failed to calc child PState in split\n";
                std::cerr << "Parent PState: " << parent(U).get_pstate() << '\n';
                std::cerr << "Child center: {" << child.center.x() << ", " << child.center.y() << ", " << child.center.z() << "}\n";
                std::cerr << "Parent center: {" << parent.center.x() << ", " << parent.center.y() << ", " << parent.center.z() << "}\n";
                std::cerr << "Parent grad x: " << parent(U).d_dx << '\n';
                std::cerr << "Parent grad y: " << parent(U).d_dy << '\n';
                std::cerr << "Parent grad z: " << parent(U).d_dz << '\n';
                std::cerr << "Child PState: " << child(U).get_pstate() << '\n';
                exit(1);
                throw std::runtime_error("bad cell");
            }
        }
    };

    distr.merge = [this](mesh::Children &children, AmrStorage::Item &parent) {
        QState sum;
        double mean_p = 0.0, mean_t = 0.0;
        for (auto &child: children) {
            sum.vec() += QState(child(U).get_pstate()).vec() * child.volume();
            mean_p += child(U).p * child.volume();
            mean_t += child(U).t * child.volume();
        }
        sum.vec() /= parent.volume();
        mean_p /= parent.volume();
        mean_t /= parent.volume();
        for (auto &b: sum.mass_frac.m_data)
            if (b < 0)
                b = 0;
        PState state(sum, mixture, mean_p, mean_t);
        parent(U).set_state(state);
//        auto [t, e] = mixture.temperature_energy_rp(parent(U).rho, parent(U).p, parent(U).mass_frac, {.T0 = mean_t});
//        parent(U).e = e;
//        parent(U).t = t;
        if (parent(U).is_bad1()) {
            std::cerr << "Failed to calc parent PState in merge\n";
            std::cerr << "QState: " << sum << '\n';
            std::cerr << "PState: " << parent(U).get_pstate() << "\n";
            exit(1);
            throw std::runtime_error("bad cell");
        }
    };

    return distr;
}

void MmFluid::set_flags(Mesh &mesh) {
    compute_grad(mesh, get_current);

    mesh.for_each([this](Cell &cell)-> void {
        double p = cell(U).p;
        double rho = cell(U).rho;
        Fractions mass_frac = cell(U).mass_frac;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            if (abs(face.neib()(U).p - p) > 0.3 * abs(p)) {
                need_split = true;
                break;
            }

            if (abs(face.neib()(U).rho - rho) > 0.1 * rho) {
                need_split = true;
                break;
            }

            // проверяем большое различие в долях веществ
            Fractions neib_mass_frac = face.neib()(U).mass_frac;
            for (int i = 0; i < mixture.size(); i++) {
                if (abs(mass_frac[i] - neib_mass_frac[i]) > 0.02) {
                    need_split = true;
                    break;
                }
            }
        }
        if (need_split) {
            cell.set_flag(1);
        } else {
            cell.set_flag(-1);
        }
    });
}

}