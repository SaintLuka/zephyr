#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/compute_grad.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace smf;
using zephyr::utils::threads;

static const SmFluid::State U = SmFluid::datatype();

[[nodiscard]] SmFluid::State SmFluid::datatype() {
    return {};
}

smf::PState get_current_sm(Cell &cell) {
    return cell(U).get_pstate();
}

smf::PState get_half_sm(Cell &cell) {
    return cell(U).half;
}

SmFluid::SmFluid(const phys::Eos &eos, Fluxes flux) : m_eos(eos) {
    m_nf = NumFlux::create(flux);
    m_CFL = 0.2;
    m_dt = std::numeric_limits<double>::max();
}

void SmFluid::set_accuracy(int acc) {
    m_acc = acc;    // 1 или 2
};

void SmFluid::update(Mesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        fluxes(mesh);
    } else {
        
        mesh.for_each([this](Cell &cell) -> void {
            assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
            assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
            assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
        });

        compute_grad(mesh, get_current_sm);

        mesh.for_each([this](Cell &cell) -> void {
            assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
            assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
            assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
        });

        fluxes_stage1(mesh);

        mesh.for_each([this](Cell &cell) -> void {
            assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
            assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
            assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
        });

        compute_grad(mesh, get_half_sm);
        
        mesh.for_each([this](Cell &cell) -> void {
            assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
            assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
            assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
        });
        
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);

    check_asserts(mesh);
}

double SmFluid::compute_dt(Mesh &mesh) {
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

    return m_dt;
}

void SmFluid::fluxes(Mesh &mesh) {
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
            auto loc_flux = m_nf->flux(zm, zp, m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.vec() += loc_flux.vec() * face.area();
        }

        // Новое значение в ячейке (консервативные переменные)
        QState q_self2 = q_self.vec() - m_dt * flux.vec() / cell.volume();

        // Новое значение примитивных переменных
        PState p_self2(q_self2, m_eos);

        cell(U).next = p_self2;
    });
}

void SmFluid::compute_grad(Mesh &mesh, const std::function<smf::PState(Cell &)> &to_state)  {
    mesh.for_each([&to_state](Cell &cell) -> void {
        auto grad = math::compute_grad_gauss<smf::PState>(cell, to_state);
        cell(U).d_dx = grad[0];
        cell(U).d_dy = grad[1];
        cell(U).d_dz = grad[2];
    });
}

void SmFluid::fluxes_stage1(Mesh &mesh)  {
    check_asserts(mesh);
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= 0.5 * m_dt / cell.volume() * calc_flux_extra(cell, true).vec();
        cell(U).half = PState(qc, m_eos);
        assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
    });
}

void SmFluid::fluxes_stage2(Mesh &mesh)  {
    check_asserts(mesh);
    mesh.for_each([this](Cell &cell) -> void {
        QState qc(cell(U).get_pstate());
        qc.vec() -= m_dt / cell.volume() * calc_flux_extra(cell, false).vec();
        cell(U).next = PState(qc, m_eos);
        assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
    });
}

smf::Flux SmFluid::calc_flux_extra(const Cell &cell, bool from_begin)  {
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

        // auto face_extra = FaceExtra::Triad(
        //         p_self, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
        //         p_neib, face.neib()(U).d_dx, face.neib()(U).d_dy, face.neib()(U).d_dz,
        //         cell_c, neib_c, face_c);

        auto face_extra = FaceExtra::ATvL(
                p_self, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                p_neib, face.neib()(U).d_dx, face.neib()(U).d_dy, face.neib()(U).d_dz,
                cell_c, neib_c, face_c);

        // рассчитываем расстояние на грани слева и справа
        PState p_minus = face_extra.m(p_self).in_local(normal);
        PState p_plus = face_extra.p(p_neib).in_local(normal);

        // std::cout << p_neib << std::endl;
        // std::cout << p_minus << std::endl;
        // std::cout << p_plus << std::endl;     

        // пересчитываем энергию и температуру
        p_minus.energy = m_eos.energy_rp(p_minus.density, p_minus.pressure);
        p_plus.energy = m_eos.energy_rp(p_plus.density, p_plus.pressure);

        assert(p_minus.pressure >= 0 && p_minus.energy >= 0);
        assert(p_plus.pressure >= 0 && p_plus.energy >= 0);

        // if (face.flag() == Boundary::WALL) {
        //     Vector3d Vn = normal * p_minus.velocity.dot(normal);
        //     p_plus.velocity = p_minus.velocity - 2 * Vn; // Vt - Vn = p_self.velocity - Vn - Vn
        // }

        // Численный поток на грани
        auto loc_flux = m_nf->flux(p_minus, p_plus, m_eos);
        loc_flux.to_global(normal);

        // Суммируем поток
        flux.vec() += loc_flux.vec() * face.area();
    }

    return flux;
}

void SmFluid::swap(Mesh &mesh) {
    size_t step = m_step;
    mesh.for_each([step](Cell &cell) -> void {
        assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
        assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
        assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
        cell(U).set_state(cell(U).next);
        cell(U).half = cell(U).next;
    });
    m_time += m_dt;
    m_step += 1;
}

void SmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

[[nodiscard]] double SmFluid::dt() const {
    return m_dt;
}

[[nodiscard]] double SmFluid::CFL() const {
    return m_CFL;
}

[[nodiscard]] double SmFluid::get_time() const {
    return m_time;
}

[[nodiscard]] size_t SmFluid::get_step() const {
    return m_step;
}

[[nodiscard]] std::string SmFluid::get_flux_name() const {
    return m_nf->get_name();
}

void SmFluid::set_acc(int acc) {
    m_acc = acc;
}

Distributor SmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](AmrStorage::Item &parent, mesh::Children &children) {
        for (auto &child: children) {
            Vector3d dr = child.center - parent.center;
            PState child_state = parent(U).get_pstate().vec() +
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
            sum.vec() += child(U).get_pstate().vec() * child.volume();
        }
        parent(U).set_state(sum.vec() / parent.volume());
        parent(U).e = m_eos.energy_rp(parent(U).rho, parent(U).p);
    };

    return distr;
}

void SmFluid::check_asserts(Mesh &mesh, std::string msg) {
    std::cout << msg << std::endl; 
    mesh.for_each([this](Cell &cell) -> void {
        assert(cell(U).get_pstate().pressure >= 0 && cell(U).get_pstate().energy >= 0);
        assert(cell(U).next.pressure >= 0 && cell(U).next.energy >= 0);
        assert(cell(U).half.pressure >= 0 && cell(U).half.energy >= 0);
    });
    std::cout << "OK" << std::endl;
}

void SmFluid::set_flags(Mesh &mesh) {

    check_asserts(mesh);

    compute_grad(mesh, get_current_sm);

    for (auto cell: mesh) {
        double p = cell(U).p;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            if (abs(face.neib()(U).p - p) >= 0.1 * p) {
                need_split = true;
                break;
            }
        }
        if (need_split) {
            cell.set_flag(1);
        } else {
            cell.set_flag(-1);
        }
    }

    check_asserts(mesh);
}

}