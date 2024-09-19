#include <zephyr/math/solver/mm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace mmf;
using zephyr::utils::threads;

static const MmFluid::State U = MmFluid::datatype();

MmFluid::State MmFluid::datatype() {
    return {};
}

MmFluid::MmFluid(const phys::Materials &mixture)
    : mixture(mixture) {
    m_nf = HLLC::create();
    m_CFL = 0.7;
    m_g = 0.0;
    m_dt = std::numeric_limits<double>::max();
}

void MmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void MmFluid::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
}

void MmFluid::set_method(Fluxes method) {
    m_nf = NumFlux::create(method);
}

void MmFluid::set_gravity(double g) {
    m_g = g;
}

double MmFluid::CFL() const {
    return m_CFL;
}

double MmFluid::dt() const {
    return m_dt;
}

mmf::PState get_current_mm(Cell &cell) {
    return cell(U).get_state();
}

PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
    if (flag != Boundary::WALL) {
        return zc;
    }

    PState zn(zc);
    Vector3d Vn = normal * zc.velocity.dot(normal);
    zn.velocity = zc.velocity - 2.0 * Vn;

    return zn;
}

void MmFluid::update(Mesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        fluxes(mesh);
    }
    else {
        compute_grad(mesh, get_current_mm);

        fluxes_stage1(mesh);

        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void MmFluid::compute_dt(Mesh &mesh) {
    double dt = mesh.min([this](Cell &cell) -> double {
        double dt = std::numeric_limits<double>::infinity();

        // скорость звука
        double c = mixture.sound_speed_rp(
                cell(U).density, cell(U).pressure, cell(U).mass_frac,
                {.T0=cell(U).temperature, .alpha=&cell(U).vol_frac});

        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell(U).velocity.dot(face.normal());

            // Максимальное по модулю СЗ
            double lambda = std::abs(vn) + c;

            // Условие КФЛ
            dt = std::min(dt, cell.volume() / (face.area() * lambda));
        }

        return dt;
    });

    m_dt = m_CFL * dt;
}

void MmFluid::fluxes(Mesh &mesh) {
    mesh.for_each([this](Cell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto &normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = face.neib(U).get_state();
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Численный поток на грани
            Flux loc_flux = m_nf->flux(zm, zp, mixture);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Новое значение примитивных переменных
        cell(U).next = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());
    });
}

void MmFluid::compute_grad(Mesh &mesh, const std::function<mmf::PState(Cell &)> &get_state)  {
    mesh.for_each([&get_state](Cell &cell) -> void {
        // auto grad = compute_grad_gauss<mmf::PState>(cell, get_state);
        auto grad = gradient::LSM_orig<mmf::PState>(cell, get_state, boundary_value);

        //auto lim_grad = gradient::limiting<mmf::PState>(cell, grad, get_state, boundary_value);

        cell(U).d_dx = grad.x;
        cell(U).d_dy = grad.y;
        cell(U).d_dz = grad.z;
    });
}

void MmFluid::fluxes_stage1(Mesh &mesh)  {
    mesh.for_each([this](Cell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto &normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа
            PState z_n;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib(U).get_state();
            } else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
            }

            auto face_extra = FaceExtra::ATvL(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_c);

            // Восстанавливаем согласованность состояния
            zm.interpolation_update(mixture);

            // Переводим в локальную систему координат
            zm.to_local(normal);

            // Численный поток на грани
            Flux loc_flux(zm);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).half = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());
    });
}

void MmFluid::fluxes_stage2(Mesh &mesh)  {
    mesh.for_each([this](Cell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор на полуслое
        PState z_c = cell(U).get_state();

        // Примитивный вектор на полуслое
        PState z_ch = cell(U).half;

        // Консервативный вектор в ячейке на прошлом шаге
        QState q_c(z_c);

        // Переменная для потока (суммирование по промежуточным)
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto &normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа (на предыдущем и на полуслое)
            PState z_n, z_nh;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib(U).get_state();
                z_nh = neib(U).half;
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::ATvL(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState zm = face_extra.m(z_ch);
            zm.interpolation_update(mixture);

            // Интерполяция на грань со стороны соседа
            PState zp;
            if (!face.is_boundary()) {
                zp = face_extra.p(z_nh);
                zp.interpolation_update(mixture);
            }
            else {
                zp = boundary_value(zm, normal, face.flag());
            }

            // Переводим в локальную систему координат
            zm.to_local(normal);
            zp.to_local(normal);

            // Численный поток на грани
            auto loc_flux = m_nf->flux(zm, zp, mixture);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).next = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());
    });
}

void MmFluid::swap(Mesh &mesh) {
    mesh.for_each([](Cell &cell) {
        cell(U).set_state(cell(U).next);
    });
}

Distributor MmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](AmrStorage::Item &parent, mesh::Children &children) {
        for (auto &child: children) {
            Vector3d dr = child.center - parent.center;
//            PState child_state = parent(U).get_state().arr() +
//                                 parent(U).d_dx.arr() * dr.x() +
//                                 parent(U).d_dy.arr() * dr.y() +
//                                 parent(U).d_dz.arr() * dr.z();
//            child_state.mass_frac.fix();
//            child_state.sync_temperature_energy_rp(mixture, {.T0 = parent(U).t});
//            child(U).set_state(child_state);
            PState shift = parent(U).d_dx.arr() * dr.x() +
                           parent(U).d_dy.arr() * dr.y() +
                           parent(U).d_dz.arr() * dr.z();
            child(U).density = parent(U).density + shift.density;
            child(U).velocity = parent(U).velocity + parent(U).density / child(U).density * shift.velocity;
            child(U).mass_frac = parent(U).mass_frac.arr() + parent(U).density / child(U).density * shift.mass_frac.arr();
            child(U).mass_frac.normalize();
            child(U).energy = parent(U).energy + (parent(U).velocity - child(U).velocity).squaredNorm() / 2 + parent(U).density / child(U).density * shift.energy;
            child(U).pressure = mixture.pressure_re(child(U).density, child(U).energy, child(U).mass_frac, {.P0 = parent(U).pressure, .T0 = parent(U).temperature});
            child(U).temperature = mixture.temperature_rp(child(U).density, child(U).pressure, child(U).mass_frac, {.T0 = parent(U).temperature});
            if (child(U).is_bad1()) {
                std::cerr << "Failed to calc child PState in split\n";
                std::cerr << "Parent PState: " << parent(U).get_state() << '\n';
                std::cerr << "Child center: {" << child.center.x() << ", " << child.center.y() << ", " << child.center.z() << "}\n";
                std::cerr << "Parent center: {" << parent.center.x() << ", " << parent.center.y() << ", " << parent.center.z() << "}\n";
                std::cerr << "Parent grad x: " << parent(U).d_dx << '\n';
                std::cerr << "Parent grad y: " << parent(U).d_dy << '\n';
                std::cerr << "Parent grad z: " << parent(U).d_dz << '\n';
                std::cerr << "Child PState: " << child(U).get_state() << '\n';
                exit(1);
                throw std::runtime_error("bad cell");
            }
        }
    };

    distr.merge = [this](mesh::Children &children, AmrStorage::Item &parent) {
        QState sum;
        double mean_p = 0.0, mean_t = 0.0;
        for (auto &child: children) {
            sum.arr() += QState(child(U).get_state()).arr() * child.volume();
            mean_p += child(U).pressure * child.volume();
            mean_t += child(U).temperature * child.volume();
        }
        sum.arr() /= parent.volume();
        mean_p /= parent.volume();
        mean_t /= parent.volume();
        for (auto &b: sum.comp_mass.m_data)
            if (b < 0)
                b = 0;
        PState state(sum, mixture, mean_p, mean_t, Fractions::NaN());
        parent(U).set_state(state);
//        auto [t, e] = mixture.temperature_energy_rp(parent(U).density, parent(U).p, parent(U).mass_frac, {.T0 = mean_t});
//        parent(U).e = e;
//        parent(U).t = t;
        if (parent(U).is_bad1()) {
            std::cerr << "Failed to calc parent PState in merge\n";
            std::cerr << "QState: " << sum << '\n';
            std::cerr << "PState: " << parent(U).get_state() << "\n";
            exit(1);
            throw std::runtime_error("bad cell");
        }
    };

    return distr;
}

void MmFluid::set_flags(Mesh &mesh) {
    compute_grad(mesh, get_current_mm);

    mesh.for_each([this](Cell &cell) -> void {
        double p = cell(U).pressure;
        double rho = cell(U).density;
        Fractions mass_frac = cell(U).mass_frac;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            if (std::abs(face.neib()(U).pressure - p) > 0.3 * abs(p)) {
                need_split = true;
                break;
            }

            if (abs(face.neib()(U).density - rho) > 0.1 * rho) {
                need_split = true;
                break;
            }

            // проверяем большое различие в долях веществ
            Fractions neib_mass_frac = face.neib()(U).mass_frac;
            for (int i = 0; i < mixture.size(); i++) {
                if (abs(mass_frac[i] - neib_mass_frac[i]) > 0.01) {
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