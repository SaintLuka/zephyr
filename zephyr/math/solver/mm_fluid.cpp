#include <zephyr/math/funcs.h>
#include <zephyr/math/solver/mm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/geom/sections.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace mmf;
using zephyr::utils::threads;

static const MmFluid::State U = MmFluid::datatype();

double MmFluid::State::vol_frac(int idx) const {
    if (mass_frac.has(idx)) {
        return mass_frac[idx] * density / densities[idx];
    } else {
        return 0.0;
    }
}

Fractions MmFluid::State::vol_fracs() const {
    Fractions alpha;
    for (int i = 0; i < Fractions::size(); ++i) {
        if (mass_frac.has(i)) {
            // Если есть массовая концентрация beta_i, значит
            // должна быть определена плотность rho_i !
            alpha[i] = mass_frac[i] * density / densities[i];
        }
    }
    return alpha;
}

MmFluid::State MmFluid::datatype() {
    return {};
}

MmFluid::MmFluid(const phys::MixturePT &mixture)
    : mixture(mixture) {
    m_nf = HLLC::create();
    m_crp_mode = CrpMode::PLIC;
    m_limiter = "MC";
    m_CFL = 0.7;
    m_g = 0.0;
    m_dt     = 1.0e+300;
    m_max_dt = 1.0e+300;
    m_split  = DirSplit::NONE;
}

void MmFluid::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void MmFluid::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
}

void MmFluid::set_method(Fluxes method) {
    if (method == Fluxes::CRP) {
        m_crp_mode = CrpMode::PLIC;

        // Поток по умолчанию
        m_nf = NumFlux::create(Fluxes::HLLC);
    }
    else {
        m_crp_mode = CrpMode::NONE;
        m_nf = NumFlux::create(method);
    }
}

void MmFluid::set_crp_mode(CrpMode mode) {
    if (mode != CrpMode::NONE) {
        m_crp_mode = mode;

        // Поток по умолчанию
        m_nf = NumFlux::create(Fluxes::HLLC);
    }
}

void MmFluid::set_limiter(const std::string& limiter) {
    m_limiter = limiter;
}

void MmFluid::set_splitting(DirSplit splitting) {
    m_split = splitting;
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

void MmFluid::set_max_dt(double dt) {
    m_max_dt = dt;
}

mmf::PState get_current_mm(Cell &cell) {
    return cell(U).get_state();
}

Fractions get_current_a(Cell &cell) {
    return cell(U).vol_fracs();
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

Fractions boundary_value_a(const Fractions &zc, const Vector3d &normal, Boundary flag) {
    return zc;
}

void MmFluid::update(Mesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    //static int counter = 0;

    if (m_split == DirSplit::NONE) {
        integrate(mesh, m_dt, Direction::ANY);
    }
    else if (m_split == DirSplit::SIMPLE) {
        //if (counter % 2 == 0) {
            //std::cout << "integrate X\n";
            integrate(mesh, m_dt, Direction::X);
        //}
        //else {
            //std::cout << "integrate Y\n";
            integrate(mesh, m_dt, Direction::Y);
        //}
        //++counter;
    }
    else if (m_split == DirSplit::STRANG) {
        integrate(mesh, 0.5 * m_dt, Direction::X);
        integrate(mesh, m_dt, Direction::Y);
        integrate(mesh, 0.5 * m_dt, Direction::X);
    }
}

void MmFluid::integrate(Mesh &mesh, double dt, Direction dir) {
    if (m_acc == 1) {
        if (m_crp_mode == CrpMode::MUSCL) {
            // Для MUSCL нужен градиент объемных долей
            fractions_grad(mesh, get_current_a);
        }
        else if (m_crp_mode == CrpMode::PLIC) {
            // Для PLIC нужна реконструкция
            interface_recovery(mesh);
        }

        // Расчет потоков
        fluxes(mesh, dt, dir);
    }
    else {
        compute_grad(mesh, get_current_mm);

        if (m_crp_mode == CrpMode::MUSCL) {
            // Для MUSCL нужен градиент объемных долей
            fractions_grad(mesh, get_current_a);
        }
        else if (m_crp_mode == CrpMode::PLIC) {
            // Для PLIC нужна реконструкция
            interface_recovery(mesh);
        }

        fluxes_stage1(mesh, dt, dir);

        fluxes_stage2(mesh, dt, dir);
    }

    swap(mesh);
}

void MmFluid::compute_dt(Mesh &mesh) {
    double dt = mesh.min([this](Cell &cell) -> double {
        double dt = std::numeric_limits<double>::infinity();

        // скорость звука
        double c = mixture.sound_speed_rP(
                cell(U).density, cell(U).pressure, cell(U).mass_frac,
                {.T0=cell(U).temperature, .rhos=&cell(U).densities});

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

    m_dt = std::min(m_CFL * dt, m_max_dt);
}

namespace {

// Выбор материала для отсечения
int master_component(const mmf::PState &zm, const mmf::PState &zp) {
    // Проверил огромное число вариантов для задачи Triple Point + PLIC.
    // Самый чистый и красивый вариант дает выбор материала по максимальной
    // объемной доле для ячейки против потока (upwind).
    //
    // Но в целом вопрос дискуссионный.

    // Скорость интерфейса
    double vi = 0.5 * (zm.u() + zp.u());
    const mmf::PState &z = vi > 0.0 ? zm : zp;

    int master = 0;
    double alpha_max = 0.0;

    for (int i = 0; i < Fractions::size(); ++i) {
        if (z.mass_frac.has(i)) {
            double alpha = z.alpha(i);

            if (alpha > alpha_max) {
                alpha_max = alpha;
                master = i;
            }
        }
    }
    return master;
}

// Параметр "доля отсечения" от грани.
double alpha_sigma(CrpMode mode, mesh::EuCell &cell, mesh::EuFace &face, int idx, double dt,
                   const mmf::PState &zm, const mmf::PState &zp, double hL, double hR) {

    double a_L = zm.alpha(idx);
    double a_R = zp.alpha(idx);

    switch (mode) {
        case CrpMode::V3: {
            return face_fraction_v3(a_L, a_R);
        }
        case CrpMode::V5: {
            return face_fraction_v5(a_L, a_R);
        }
        case CrpMode::VS: {
            return face_fraction_s(a_L, a_R);
        }
        case CrpMode::MUSCL: {
            // Вытащить градиент данной величины.
            // определить направление движения
            // сделать спуск по градиенту
            return NAN;
        }
        default: {
            // Фактически CrpMode::PLIC
            double vi = 0.5 * (zm.u() + zp.u()); // Скорость интерфейса
            bool upwind = vi > 0.0;
            double CFL = dt * std::abs(vi) / (upwind ? hL : hR);
            double cos = face.normal().dot(upwind ? cell(U).n[idx] : -face.neib(U).n[idx]);
            double a_sig = geom::average_flux(upwind ? a_L : a_R, cos, CFL);
            return math::between(a_sig, a_L, a_R);
        }
    }
}
}

// Посчитать многоматериальный поток с методом CRP с расщеплением относительно материала с номером iA.
// Величина a_sig -- доля деления грани материалом с номером iA.
Flux MmFluid::calc_crp_flux(const PState& zL, const PState& zR, double hL, double hR, int iA, double a_sig, double dt) {
    double aL = zL.alpha(iA);
    double aR = zR.alpha(iA);

    // Либо не содержат iA, либо полностью заполнены iA
    if ((aL == 0.0 || aL == 1.0) && (aR == 0.0 || aR == 1.0)) {
        return m_nf->flux(zL, zR, mixture);
    }

    // Другие случаи пока не рассматриваем a_sig ∈ [a_min, a_max]
    assert(std::min(aL, aR) <= a_sig && a_sig <= std::max(aL, aR));

    double delta_L = (aL < a_sig ? aL / a_sig : (1.0 - aL) / (1.0 - a_sig)) * hL;
    double delta_R = (aR < a_sig ? aR / a_sig : (1.0 - aR) / (1.0 - a_sig)) * hR;

    // aL = 0, aR ∈ (0, 1)
    if (aL == 0.0) {
        auto[zRA, zRB] = zR.split(mixture, iA);
        auto& zLB = zL;

        auto flux_A = m_nf->flux(zLB, zRA, mixture);
        auto flux_B = CrpFlux::inverse(zLB, zRB, zRA, mixture, delta_R, dt);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL = 1, aR ∈ (0, 1)
    if (aL == 1.0) {
        auto[zRA, zRB] = zR.split(mixture, iA);
        auto& zLA = zL;

        auto flux_A = CrpFlux::inverse(zLA, zRA, zRB, mixture, delta_R, dt);
        auto flux_B = m_nf->flux(zLA, zRB, mixture);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL ∈ (0, 1), aR == 0
    if (aR == 0.0) {
        auto[zLA, zLB] = zL.split(mixture, iA);
        auto& zRB = zR;

        auto flux_A = m_nf->flux(zLA, zRB, mixture);
        auto flux_B = CrpFlux::classic(zLA, zLB, zRB, mixture, delta_L, dt);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL ∈ (0, 1), aR == 1
    if (aR == 1.0) {
        auto[zLA, zLB] = zL.split(mixture, iA);
        auto& zRA = zR;

        auto flux_A = CrpFlux::classic(zLB, zLA, zRA, mixture, delta_L, dt);
        auto flux_B = m_nf->flux(zLB, zRA, mixture);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // Две строго смешаные ячейки
    assert((0.0 < aL && aL < 1.0) && (0.0 < aR && aR < 1.0));

    // Расщепление
    auto[zLA, zLB] = zL.split(mixture, iA);
    auto[zRA, zRB] = zR.split(mixture, iA);

    mmf::Flux flux;
    if (aL < aR) {
        auto flux_A = CrpFlux::classic(zLB, zLA, zRA, mixture, delta_L, dt);
        auto flux_B = CrpFlux::inverse(zLB, zRB, zRA, mixture, delta_R, dt);
        flux = a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }
    else {
        auto flux_A = CrpFlux::inverse(zLA, zRA, zRB, mixture, delta_R, dt);
        auto flux_B = CrpFlux::classic(zLA, zLB, zRB, mixture, delta_L, dt);
        flux = a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    return flux;
}

Flux MmFluid::calc_flux(mesh::EuCell& cell, mesh::EuFace& face,
                        const PState& z_L, const PState& z_R,
                        double h_L, double h_R, double dt) {
    // Никакого CRP, обычный поток
    if (m_crp_mode == CrpMode::NONE) {
        return m_nf->flux(z_L, z_R, mixture);
    }

    int idx_L = z_L.mass_frac.index();
    int idx_R = z_R.mass_frac.index();

    // Пара чистых материалов
    if (idx_L >= 0 && idx_R >= 0) {
        return m_nf->flux(z_L, z_R, mixture);
    }

    // Выбрать мастер-материал
    int iA = master_component(z_L, z_R);

    // Объемная доля отсечения
    double a_sig = alpha_sigma(m_crp_mode, cell, face, iA, dt, z_L, z_R, h_L, h_R);

    // Двумерный CRP поток
    return calc_crp_flux(z_L, z_R, h_L, h_R, iA, a_sig, dt);
}

void MmFluid::fluxes(Mesh &mesh, double dt, Direction dir) {
    mesh.for_each([this, dt, dir](Cell &cell) {
        // Объем ячейки
        double V_c = cell.volume();

        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        const Box box = {.vmin={0.0030, 0.0060, 0.0}, .vmax={0.0031, 0.0061, 0.0}};

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces(dir)) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            double V_n;
            if (!face.is_boundary()) {
                z_n = face.neib(U).get_state();
                V_n = face.neib().volume();
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
                V_n = V_c;
            }

            // Значение на грани со стороны ячейки
            PState z_m = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState z_p = z_n.in_local(normal);

            // Численный поток на грани
            double S = face.area();

            double h_L = V_c / S;
            double h_R = V_n / S;

            if (box.inside(cell.center())) {
                std::cout << "stop\n";

            }

            Flux loc_flux = calc_flux(cell, face, z_m, z_p, h_L, h_R, dt);
            loc_flux.to_global(normal);

            if (loc_flux.arr().hasNaN()) {

                Flux loc2 = calc_flux(cell, face, z_m, z_p, h_L, h_R, dt);
            }

            if (box.inside(cell.center())) {
                std::cout << "z_m: " << z_m << "\n";
                std::cout << "z_p: " << z_p << "\n";

                std::cout << "loc_flux: " << loc_flux << "\n";
            }

            // Суммируем поток
            flux.arr() += loc_flux.arr() * S;
        }


        // Консервативный вектор в ячейке
        QState q_c(z_c);

        if (box.inside(cell.center())) {
            std::cout << "\tz_c: " << z_c << "\n";
            std::cout << "\tq_c: " << q_c << "\n";
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (dt / V_c) * flux.arr();

        // Новое значение примитивных переменных
        cell(U).next = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());

        if (box.inside(cell.center())) {
            std::cout << "\tq_c: " << q_c << "\n";
            std::cout << "\tz_h: " << cell(U).next << "\n";
        }
    });
}

void MmFluid::compute_grad(Mesh &mesh, const GetState<PState> &get_state)  {
    mesh.for_each([this, &get_state](Cell &cell) {

        // Для смешаных ячеек в CRP режиме производные равны нулю
        if (m_crp_mode != CrpMode::NONE) {
            bool set_zero = cell(U).mass_frac.index() < 0;
            if (!set_zero) {
                for (auto face: cell.faces()) {
                    if (face.neib(U).mass_frac.index() < 0) {
                        set_zero = true;
                        break;
                    }
                }
            }
            if (set_zero) {
                cell(U).d_dx = PState::Zero();
                cell(U).d_dy = PState::Zero();
                cell(U).d_dz = PState::Zero();
                return;
            }
        }

        auto grad = gradient::LSM<PState>(cell, get_state, boundary_value);
        grad = gradient::limiting<PState>(cell, m_limiter, grad, get_state, boundary_value);

        cell(U).d_dx = grad.x;
        cell(U).d_dy = grad.y;
        cell(U).d_dz = grad.z;
    });
}

void MmFluid::fractions_grad(Mesh& mesh, const GetState<Fractions>& get_state) {
    mesh.for_each([this, &get_state](Cell &cell) {
        auto grad = gradient::LSM<Fractions>(cell, get_current_a, boundary_value_a);
        grad = gradient::limiting<Fractions>(cell, m_limiter, grad, get_current_a, boundary_value_a);

        auto& grad_a = cell(U).grad_a;
        for (int i = 0; i < Fractions::size(); ++i) {
            grad_a[i] = {grad.x[i], grad.y[i], grad.z[i]};
        }
    });
}

void MmFluid::interface_recovery(Mesh &mesh) {
    // Сделаю пока простую схему для квадратов,
    // точка p не восстанавливается, плевать на неё.

    mesh.for_each([](Cell& cell) {
        // Чистый материал, нечего восстанавливать
        if (cell(U).mass_frac.is_pure()) {
            cell(U).n = VectorSet();
            return;
        }

        VectorSet ns;

        Fractions a_c = cell(U).vol_fracs();
        for (auto face: cell.faces()) {
            Vector3d S = face.area() * face.normal();

            // На границе возвращает саму ячейку
            Fractions a_n = face.neib(U).vol_fracs();
            for (int i = 0; i < Fractions::max_size; ++i) {
                if (a_c.has(i)) {
                    ns[i] -= face_fraction(a_c[i], a_n[i]) * S;
                }
            }
        }

        for (int i = 0; i < Fractions::max_size; ++i) {
            if (a_c.has(i)) {
                ns[i].normalize();
            }
        }

        cell(U).n = ns;
    });
}

void MmFluid::fluxes_stage1(Mesh &mesh, double dt, Direction dir)  {
    mesh.for_each([this, dt, dir](Cell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        // Смешаные ячейки продвигаем
        if (z_c.mass_frac.index() < 0) {
            cell(U).half = z_c;
            return;
        }

        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces(dir)) {
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

            auto face_extra = FaceExtra::Direct(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState z_m = face_extra.m(z_c);

            // Восстанавливаем согласованность состояния
            z_m.interpolation_update(mixture);

            // Переводим в локальную систему координат
            z_m.to_local(normal);

            // Численный поток на грани
            Flux loc_flux(z_m);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * face.area();
        }

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).half = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());
    });
}

void MmFluid::fluxes_stage2(Mesh &mesh, double dt, Direction dir)  {
    mesh.for_each([this, dt, dir](Cell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Объем ячейки
        double V_c = cell.volume();

        // Примитивный вектор
        PState z_c = cell(U).get_state();

        // Примитивный вектор (на полушаге)
        PState z_ch = cell(U).half;

        // Переменная для потока (суммирование по промежуточным)
        Flux flux;
        for (auto &face: cell.faces(dir)) {
            // Внешняя нормаль и центр грани
            auto &normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа (на предыдущем и на полушаге)
            PState z_n, z_nh;
            Vector3d neib_c;
            double V_n;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                V_n    = neib.volume();
                z_n  = neib(U).get_state();
                z_nh = neib(U).half;
            }
            else {
                neib_c = face.symm_point(cell_c);
                V_n  = V_c;
                z_n  = boundary_value(z_c, normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell(U).d_dx, cell(U).d_dy, cell(U).d_dz,
                    z_n, neib(U).d_dx, neib(U).d_dy, neib(U).d_dz,
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            PState z_m = face_extra.m(z_ch);
            z_m.interpolation_update(mixture);

            // Интерполяция на грань со стороны соседа
            PState z_p;
            if (!face.is_boundary()) {
                z_p = face_extra.p(z_nh);
                z_p.interpolation_update(mixture);
            }
            else {
                z_p = boundary_value(z_m, normal, face.flag());
            }

            // Переводим в локальную систему координат
            z_m.to_local(normal);
            z_p.to_local(normal);

            double S = face.area();

            double h_L = V_c / S;
            double h_R = V_n / S;

            // Численный поток на грани
            auto loc_flux = calc_flux(cell, face, z_m, z_p, h_L, h_R, dt);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * S;
        }

        // Консервативный вектор в ячейке на прошлом шаге
        QState q_c(z_c);

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).next = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());
    });
}

void MmFluid::swap(Mesh &mesh) {
    mesh.for_each([](Cell &cell) {
        cell(U).set_state(cell(U).next);
    });
}

AmrStorage MmFluid::body(Mesh& mesh, int idx) const {
    int count = 0;
    for (auto cell: mesh) {
        double alpha = cell(U).vol_frac(idx);
        if (std::isnan(alpha) || alpha < 1.0e-12) {
            continue;
        }
        ++count;
    }

    AmrStorage cells(count);

    count = 0;
    for (auto cell: mesh) {
        double alpha = cell(U).vol_frac(idx);

        if (std::isnan(alpha) || alpha < 1.0e-12) {
            continue;
        }

        Vector3d normal = cell(U).n[idx];

        if (alpha < 1.0 - 1.0e-12) {
            if (normal.isZero()) {
                double d = 0.5 * std::sqrt(alpha * cell.volume());
                Quad quad = {
                        cell.center() + Vector3d{-d, -d, 0.0},
                        cell.center() + Vector3d{+d, -d, 0.0},
                        cell.center() + Vector3d{-d, +d, 0.0},
                        cell.center() + Vector3d{+d, +d, 0.0},
                };
                cells[count] = mesh::AmrCell(quad);
            }
            else {
                auto poly = cell.polygon();
                Vector3d point = poly.find_section(normal, alpha);
                auto part = poly.clip(point, normal);
                cells[count] = mesh::AmrCell(part);
            }
        }
        else {
            cells[count] = cell.geom();
        }
        ++count;
    }

    return cells;
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
//            child_state.sync_temperature_energy_rP(mixture, {.T0 = parent(U).t});
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
            child(U).temperature = mixture.temperature_rP(child(U).density, child(U).pressure, child(U).mass_frac, {.T0 = parent(U).temperature});
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
            sum.arr() += QState(child(U).get_state()).arr() * child.volume;
            mean_p += child(U).pressure * child.volume;
            mean_t += child(U).temperature * child.volume;
        }
        sum.arr() /= parent.volume;
        mean_p /= parent.volume;
        mean_t /= parent.volume;
        for (auto &b: sum.mass_frac.m_data)
            if (b < 0)
                b = 0;
        PState state(sum, mixture, mean_p, mean_t, Fractions::NaN());
        parent(U).set_state(state);
//        auto [t, e] = mixture.temperature_energy_rP(parent(U).density, parent(U).p, parent(U).mass_frac, {.T0 = mean_t});
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