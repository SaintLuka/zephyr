#include <zephyr/math/funcs.h>
#include <zephyr/math/solver/mm_fluid.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>

namespace zephyr::math {

using namespace geom;
using namespace mesh;
using namespace mmf;
using utils::threads;

MmFluid::MmFluid(const MixturePT &mixture)
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

MmFluid::Parts MmFluid::add_types(EuMesh& mesh) {
    part.init = mesh.add<PState>("init");
    part.d_dx = mesh.add<PState>("d_dx");
    part.d_dy = mesh.add<PState>("d_dy");
    part.d_dz = mesh.add<PState>("d_dz");
    part.half = mesh.add<PState>("half");
    part.next = mesh.add<PState>("next");

    part.n = mesh.add<VectorSet>("normals");
    part.p = mesh.add<VectorSet>("origins");
    part.grad_a = mesh.add<VectorSet>("vol_grad");
    return part;
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

void MmFluid::update(EuMesh &mesh) {
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

void MmFluid::integrate(EuMesh &mesh, double dt, Direction dir) {
    if (m_acc == 1) {
        if (m_crp_mode == CrpMode::MUSCL) {
            // Для MUSCL нужен градиент объемных долей
            fractions_grad(mesh, part.init);
        }
        else if (m_crp_mode == CrpMode::PLIC) {
            // Для PLIC нужна реконструкция
            interface_recovery(mesh);
        }

        // Расчет потоков
        fluxes(mesh, dt, dir);
    }
    else {
        compute_grad(mesh, part.init);

        if (m_crp_mode == CrpMode::MUSCL) {
            // Для MUSCL нужен градиент объемных долей
            fractions_grad(mesh, part.init);
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

void MmFluid::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell &cell) -> double {
        double dt = std::numeric_limits<double>::max();

        // скорость звука
        const auto&z = cell[part.init];
        double c = mixture.sound_speed_rP(z.density, z.pressure, z.mass_frac,
                                          {.T0=z.temperature, .rhos=&z.densities});

        for (auto &face: cell.faces()) {
            // Нормальная составляющая скорости
            double vn = cell[part.init].velocity.dot(face.normal());

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
double alpha_sigma(CrpMode mode, EuCell &cell, const MmFluid::Parts& part, EuFace &face, int idx, double dt,
                   const PState &zm, const PState &zp, double hL, double hR) {

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
            double cos = face.normal().dot(upwind ? cell[part.n][idx] : -face.neib(part.n)[idx]);
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

Flux MmFluid::calc_flux(EuCell& cell, EuFace& face,
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
    double a_sig = alpha_sigma(m_crp_mode, cell, part, face, iA, dt, z_L, z_R, h_L, h_R);

    // Двумерный CRP поток
    return calc_crp_flux(z_L, z_R, h_L, h_R, iA, a_sig, dt);
}

void MmFluid::fluxes(EuMesh &mesh, double dt, Direction dir) {
    mesh.for_each([this, dt, dir](EuCell &cell) {
        // Объем ячейки
        double V_c = cell.volume();

        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

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
                z_n = face.neib(part.init);
                V_n = face.neib_volume();
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
        cell[part.next] = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());

        if (box.inside(cell.center())) {
            std::cout << "\tq_c: " << q_c << "\n";
            std::cout << "\tz_h: " << cell[part.next] << "\n";
        }
    });
}

void MmFluid::compute_grad(EuMesh &mesh, Storable<PState> U)  {
    mesh.for_each([this, U](EuCell &cell) {
        // Для смешанных ячеек в CRP режиме производные равны нулю
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
                cell[part.d_dx] = PState::Zero();
                cell[part.d_dy] = PState::Zero();
                cell[part.d_dz] = PState::Zero();
                return;
            }
        }

        auto grad = gradient::LSM<PState>(cell, U, boundary_value);
        grad = gradient::limiting<PState>(cell, m_limiter, grad, U, boundary_value);

        cell[part.d_dx] = grad.x;
        cell[part.d_dy] = grad.y;
        cell[part.d_dz] = grad.z;
    });
}

void MmFluid::fractions_grad(EuMesh& mesh, Storable<PState> U) {
    auto get_vol_fracs = [U](EuCell& cell) -> Fractions {
        return cell[U].vol_fracs();
    };

    mesh.for_each([this, &get_vol_fracs](EuCell &cell) {
        auto grad = gradient::LSM<Fractions>(cell, get_vol_fracs, boundary_value_a);
        grad = gradient::limiting<Fractions>(cell, m_limiter, grad, get_vol_fracs, boundary_value_a);

        auto& grad_a = cell[part.grad_a];
        for (int i = 0; i < Fractions::size(); ++i) {
            grad_a[i] = {grad.x[i], grad.y[i], grad.z[i]};
        }
    });
}

void MmFluid::interface_recovery(EuMesh &mesh) {
    // Сделаю пока простую схему для квадратов,
    // точка p не восстанавливается, плевать на неё.

    mesh.for_each([this](EuCell& cell) {
        // Чистый материал, нечего восстанавливать
        if (cell[part.init].mass_frac.is_pure()) {
            cell[part.n] = VectorSet{};
            return;
        }

        VectorSet ns;

        Fractions a_c = cell[part.init].vol_fracs();
        for (auto face: cell.faces()) {
            Vector3d S = face.area() * face.normal();

            // На границе возвращает саму ячейку
            Fractions a_n = face.neib(part.init).vol_fracs();
            for (int i = 0; i < Fractions::max_size; ++i) {
                if (a_c.has(i)) {
                    ns[i] -= face_fraction(a_c[i], a_n[i]) * S;
                }
            }
        }

        for (int i = 0; i < Fractions::max_size; ++i) {
            if (a_c.has(i)) { ns[i].normalize(); }
        }
        cell[part.n] = ns;
    });
}

void MmFluid::fluxes_stage1(EuMesh &mesh, double dt, Direction dir)  {
    mesh.for_each([this, dt, dir](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Смешанные ячейки продвигаем
        if (z_c.mass_frac.index() < 0) {
            cell[part.half] = z_c;
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
                z_n = neib[part.init];
            } else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
            }

            auto face_extra = FaceExtra::Direct(
                    z_c, cell[part.d_dx], cell[part.d_dy], cell[part.d_dz],
                    z_n, neib[part.d_dx], neib[part.d_dy], neib[part.d_dz],
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
        QState q_c{z_c};

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (0.5 * dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell[part.half] = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());
    });
}

void MmFluid::fluxes_stage2(EuMesh &mesh, double dt, Direction dir)  {
    mesh.for_each([this, dt, dir](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Объем ячейки
        double V_c = cell.volume();

        // Примитивный вектор
        PState z_c = cell[part.init];

        // Примитивный вектор (на полушаге)
        PState z_ch = cell[part.half];

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
                z_n  = neib[part.init];
                z_nh = neib[part.half];
            }
            else {
                neib_c = face.symm_point(cell_c);
                V_n  = V_c;
                z_n  = boundary_value(z_c, normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell[part.d_dx], cell[part.d_dy], cell[part.d_dz],
                    z_n, neib[part.d_dx], neib[part.d_dy], neib[part.d_dz],
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
        QState q_c{z_c};

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell[part.next] = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.rhos());
    });
}

void MmFluid::swap(EuMesh &mesh) {
    mesh.swap(part.init, part.next);
}

EuMesh MmFluid::body(EuMesh& mesh, int idx) const {
    throw std::runtime_error("MmFluid::body: Not implemented");
#if 0
    int count = 0;
    for (auto cell: mesh) {
        double alpha = cell(U).vol_frac(idx);
        if (std::isnan(alpha) || alpha < 1.0e-12) {
            continue;
        }
        ++count;
    }

    EuCell cells(count);

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
#endif
}

Distributor MmFluid::distributor() const {
    Distributor distr;

    distr.split = [this](const EuCell &parent, Children &children) {
        PState zp = parent[part.init];

        for (auto child: children) {
            Vector3d dr = child.center() - parent.center();
            PState delta = parent[part.d_dx].arr() * dr.x() +
                           parent[part.d_dy].arr() * dr.y() +
                           parent[part.d_dz].arr() * dr.z();

            PState zc;
            zc.density   = zp.density + delta.density;
            zc.velocity  = zp.velocity + (zp.density / zc.density) * delta.velocity;
            zc.mass_frac = zp.mass_frac.arr() + (zp.density / zc.density) * delta.mass_frac.arr();
            zc.mass_frac.normalize();
            zc.energy    = zp.energy + 0.5 * (zp.velocity - zc.velocity).squaredNorm() + (zp.density / zc.density) * delta.energy;
            zc.pressure = mixture.pressure_re(zc.density, zc.energy, zc.mass_frac, {.P0 = zp.pressure, .T0 = zp.temperature});
            zc.temperature = mixture.temperature_rP(zc.density, zc.pressure, zc.mass_frac, {.T0 = zp.temperature});

            if (zc.is_bad()) {
                std::cerr << "Failed to calc child PState in split\n";
                std::cerr << "Parent PState: " << zp << '\n';
                std::cerr << "Child  center: {" << child.center().transpose() << "}\n";
                std::cerr << "Parent center: {" << parent.center().transpose() << "}\n";
                std::cerr << "Parent grad x: " << parent[part.d_dx] << '\n';
                std::cerr << "Parent grad y: " << parent[part.d_dy] << '\n';
                std::cerr << "Parent grad z: " << parent[part.d_dz] << '\n';
                std::cerr << "Child  PState: " <<  zc << '\n';
                throw std::runtime_error("bad cell");
            }
            child[part.init] = zc;
        }
    };

    distr.merge = [this](const Children &children, EuCell &parent) {
        QState sum;
        double mean_p = 0.0, mean_t = 0.0;
        for (auto child: children) {
            QState qc(child[part.init]);

            sum.arr() += qc.arr() * child.volume();
            mean_p += child[part.init].P() * child.volume();
            mean_t += child[part.init].T() * child.volume();
        }
        sum.arr() /= parent.volume();
        mean_p /= parent.volume();
        mean_t /= parent.volume();

        for (auto &b: sum.mass_frac.m_data) {
            b = std::max(0.0, b);
        }

        PState state(sum, mixture, mean_p, mean_t, Fractions::NaN());
        if (state.is_bad()) {
            std::cerr << "Failed to calc parent PState in merge\n";
            std::cerr << "QState: " << sum << '\n';
            std::cerr << "PState: " << state << "\n";
            throw std::runtime_error("bad cell");
        }

        parent[part.init] = state;
    };

    return distr;
}

void MmFluid::set_flags(EuMesh &mesh) {
    compute_grad(mesh, part.init);

    mesh.for_each([this](EuCell &cell) -> void {
        double rho = cell[part.init].density;
        double p   = cell[part.init].pressure;

        Fractions mass_frac = cell[part.init].mass_frac;
        bool need_split = false;
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            // проверяем большой перепад давлений
            if (std::abs(face.neib(part.init).pressure - p) > 0.3 * abs(p)) {
                need_split = true;
                break;
            }

            if (abs(face.neib(part.init).density - rho) > 0.1 * rho) {
                need_split = true;
                break;
            }

            // проверяем большое различие в долях веществ
            Fractions neib_mass_frac = face.neib(part.init).mass_frac;
            for (int i = 0; i < mixture.size(); i++) {
                if (abs(mass_frac[i] - neib_mass_frac[i]) > 0.01) {
                    need_split = true;
                    break;
                }
            }
        }

        cell.set_flag(need_split ? 1 : -1);
    });
}

}