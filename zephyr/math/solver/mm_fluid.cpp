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
    return cell(U).vol_frac;
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

    if (m_split == DirSplit::NONE) {
        integrate(mesh, m_dt, Direction::ANY);
    }
    else if (m_split == DirSplit::SIMPLE) {
        integrate(mesh, m_dt, Direction::X);
        integrate(mesh, m_dt, Direction::Y);
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
            // Для PLIC нужна, собственно реконструкция
            interface_recovery(mesh);
        }

        // Расчет потоков
        fluxes(mesh, dt, dir);
    }
    else {
        compute_grad(mesh, get_current_mm);

        fluxes_stage1(mesh, dt, dir);

        fluxes_stage2(mesh, dt, dir);
    }

    // Обновляем слои
    swap(mesh);
}

void MmFluid::compute_dt(Mesh &mesh) {
    double dt = mesh.min([this](Cell &cell) -> double {
        double dt = std::numeric_limits<double>::infinity();

        // скорость звука
        double c = mixture.sound_speed_rP(
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

    m_dt = std::min(m_CFL * dt, m_max_dt);
}

using phys::Eos;

struct Point {
    double x, t;
};

struct Char {
    double x, t, S;

    Point& p() {
        return (Point &) x;
    }

    double edge_t() const {
        double ts = (S * t - x) / S;
        if (ts <= 0.0) {
            return std::numeric_limits<double>::infinity();
        } else {
            return ts;
        }
    }
};

// Пересечение характеристик, находятся только пересечения во времени дальше,
// чем время обеих точек исходных характеристик. В обратном случае считается,
// что характеристики пересекаются на бесконечности t = +inf.
Point operator&(const Char& c1, const Char& c2) {
    const double inf = std::numeric_limits<double>::infinity();

    if (c1.S == c2.S) {
        return {.x = 0.5 * (c1.x + c2.x), .t = inf};
    }

    double t = (c2.x - c1.x + c1.S * c1.t - c2.S * c2.t) / (c1.S - c2.S);

    if (t <= std::max(c1.t, c2.t)) {
        return {.x = 0.5 * (c1.x + c2.x), .t = inf};
    }

    double x = 0.5 * (c1.x + c2.x + c1.S * (t - c1.t) + c2.S * (t - c2.t));

    return {.x = x, .t = t};
}

// Классическая задача для CRP. Только два материала в паре ячеек.
// mixture -- список материалов
// iA, iB -- индексы материалов.
// zLA, zLB -- чистые состояния в ячейке слева.
// zRB -- чистое состояние в ячейке справа.
// delta > 0 -- отступ от грани контактной границы в левой ячейке.
// dt -- шаг интегрирования
Flux crp_flux_1(const smf::PState& zLA, const smf::PState& zLB, const smf::PState& zRB,
                const MixturePT& mixture, int iA, int iB, double delta, double dt) {
    auto &eosA = mixture[iA];
    auto &eosB = mixture[iB];

    smf::QState qLB(zLB);
    smf::Flux   fLB(zLB);
    smf::QState Q_R(zRB);
    smf::Flux   F_R(zRB);

    // Характеристики из HLL
    auto[S_0L, S_0R, Q_s1, F_s1] = HLL::wave_config(eosB, qLB, fLB, Q_R, F_R);
    Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};
    Char C_0R = {.x = 0.0, .t = 0.0, .S = S_0R};

    // Первый поток через грань
    Flux F1(F_s1, iB);

    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = zLA.u()};

    // Первое взаимодействие
    Point O1 = C_0 & C_0L;

    // Нет взаимодействия за dt
    if (O1.t >= dt) {
        return F1;
    }

    smf::QState Q_L(zLA);
    smf::Flux   F_L(zLA);

    // Пара характеристик из HLLC
    auto[S_1L, S_1C, S_1R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(eosA, Q_L, F_L, eosB, Q_s1, F_s1);
    Char C_1C = {.x = O1.x, .t = O1.t, S_1C};
    Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

    double tau1 = C_1R.edge_t() / dt;
    if (tau1 >= 1.0) {
        return F1;
    }

    // Второй поток
    Flux F2(F_s2R, iB);

    // Второе взаимодействие
    //Point O2 = C_1R & C_0R;

    double tau2 = C_1C.edge_t() / dt;


    if (tau2 >= 1.0) {
        return tau1 * F1.arr() + (1.0 - tau1) * F2.arr();
    }

    Flux F3(F_s2L, iA);

    /*
    std::cout << "\t\tyep, tau1, tau2: " << tau1 << " " << tau2 << "\n";
    std::cout << "\t   F1: " << F1 << "\n";
    std::cout << "\t   F2: " << F2 << "\n";
    std::cout << "\t   F3: " << F3 << "\n";
    */

    return tau1 * F1.arr() + (tau2 - tau1) * F2.arr() + (1.0 - tau2) * F3.arr();


    /*
    auto[S_2L, S_2R, Q_s3, F_s3] = HLL::wave_config(eosB, Q_s2R, F_s2R, Q_R, F_R);
    Char C_2L = {.x = O2.x, .t = O1.t, S_2L};
    Char C_2R = {.x = O2.x, .t = O1.t, S_2R};

    Point O3 = C_1C & C_2L;

    */
}

Flux crp_flux_2(const smf::PState& zLA, const smf::PState& zRA, const smf::PState& zRB,
                const MixturePT& mixture, int iA, int iB, double delta, double dt) {
    auto& zLA_i = const_cast<smf::PState&>(zLA);
    auto& zRA_i = const_cast<smf::PState&>(zRA);
    auto& zRB_i = const_cast<smf::PState&>(zRB);

    zLA_i.inverse();
    zRA_i.inverse();
    zRB_i.inverse();

    auto flux = crp_flux_1(zRB_i, zRA_i, zLA_i, mixture, iB, iA, delta, dt);
    flux.inverse();

    zLA_i.inverse();
    zRA_i.inverse();
    zRB_i.inverse();

    return flux;

}

static double face_fraction_v3(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    if (a1 + a2 < 1.0) return a_max;
    if (a1 + a2 > 1.0) return a_min;

    return 0.5 * (a1 + a2);
}

static double face_fraction_v5(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    return 0.5 * (a1 + a2);
}

// Формула Серёжкина
static double face_fraction_s(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);

    double a_sig = a_min / (1.0 - (a_max - a_min));

    // Случай a_min = 0, a_max = 1
    if (std::isnan(a_sig)) {
        a_sig = 0.5;
    }

    return  between(a_sig, a_min, a_max);
}

Flux MmFluid::calc_flux(const PState& zL, const PState& zR, double hL, double hR, double a_sig, double dt) {
    if (m_crp_mode == CrpMode::NONE) {
        return m_nf->flux(zL, zR, mixture);
    }

    int iA = 0;
    int iB = 1;

    double aL = zL.vol_frac[iA];
    double aR = zR.vol_frac[iA];

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
        auto flux_B = crp_flux_2(zLB.to_smf(), zRB.to_smf(), zRA.to_smf(), mixture, iB, iA, delta_R, dt);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL = 1, aR ∈ (0, 1)
    if (aL == 1.0) {
        auto[zRA, zRB] = zR.split(mixture, iA);
        auto& zLA = zL;

        auto flux_A = crp_flux_2(zLA.to_smf(), zRA.to_smf(), zRB.to_smf(), mixture, iA, iB, delta_R, dt);
        auto flux_B = m_nf->flux(zLA, zRB, mixture);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL ∈ (0, 1), aR == 0
    if (aR == 0.0) {
        auto[zLA, zLB] = zL.split(mixture, iA);
        auto& zRB = zR;

        auto flux_A = m_nf->flux(zLA, zRB, mixture);
        auto flux_B = crp_flux_1(zLA.to_smf(), zLB.to_smf(), zRB.to_smf(), mixture, iA, iB, delta_L, dt);

        return a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    // aL ∈ (0, 1), aR == 1
    if (aR == 1.0) {
        auto[zLA, zLB] = zL.split(mixture, iA);
        auto& zRA = zR;

        auto flux_A = crp_flux_1(zLB.to_smf(), zLA.to_smf(), zRA.to_smf(), mixture, iB, iA, delta_L, dt);
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
        auto flux_A = crp_flux_1(zLB.to_smf(), zLA.to_smf(), zRA.to_smf(), mixture, iB, iA, delta_L, dt);
        auto flux_B = crp_flux_2(zLB.to_smf(), zRB.to_smf(), zRA.to_smf(), mixture, iB, iA, delta_R, dt);
        flux = a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }
    else {
        auto flux_A = crp_flux_2(zLA.to_smf(), zRA.to_smf(), zRB.to_smf(), mixture, iA, iB, delta_R, dt);
        auto flux_B = crp_flux_1(zLA.to_smf(), zLB.to_smf(), zRB.to_smf(), mixture, iA, iB, delta_L, dt);
        flux = a_sig * flux_A.arr() + (1.0 - a_sig) * flux_B.arr();
    }

    return flux;
}

void MmFluid::fluxes(Mesh &mesh, double dt, Direction dir) {
    mesh.for_each([this, dt, dir](Cell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        double V_c = cell.volume();

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces(dir)) {
            // Внешняя нормаль
            auto &normal = face.normal();

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
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Численный поток на грани
            double S = face.area();
            double hL = V_c / S;
            double hR = V_n / S;

            double a_sig = face_fraction(zm.vol_frac[0], zp.vol_frac[0]);
            bool positive = zm.u() + zp.u() > 0.0;
            double CFL = 0.5 * dt * std::abs(zm.u() + zp.u()) / hL;
            double cos = (positive ? cell(U).n[0] : face.neib(U).n[0]).dot(face.normal());
            if (!positive) {
                cos *= -1.0;
            }
            a_sig = geom::average_flux(positive ? zm.vol_frac[0] : zp.vol_frac[0], cos, CFL);
            a_sig = between(a_sig, zm.vol_frac[0], zp.vol_frac[0]);

            Flux loc_flux = calc_flux(zm, zp, hL, hR, a_sig, dt);
            loc_flux.to_global(normal);

            // Суммируем поток
            flux.arr() += loc_flux.arr() * S;
        }

        // Обновляем значение в ячейке (консервативные переменные)
        q_c.arr() -= (dt / V_c) * flux.arr();

        // Новое значение примитивных переменных
        cell(U).next = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());

        /*
        if (cell(U).next.P() < 0.0) {

            PState zsd(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());

            std::cout << "Nfr: P " << z_c.P() << " " << z_c.T() <<" " << z_c.alpha() << "\n";
            std::cout << "OLOLO SUK BLYAD " << q_c << "\n";
            std::cout << zsd << "\n";
            
            PState zsp(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());
        }
         */
    });
}

void MmFluid::compute_grad(Mesh &mesh, const GetState<PState> &get_state)  {
    mesh.for_each([this, &get_state](Cell &cell) {
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

        cell(U).d_dx.vol_frac = grad.x;
        cell(U).d_dy.vol_frac = grad.y;
        cell(U).d_dz.vol_frac = grad.z;
    });
}

void MmFluid::interface_recovery(Mesh &mesh) {
    // Сделаю пока простую схему для квадратов,
    // точка p не восстанавливается, плевать на неё.

    mesh.for_each([](Cell& cell) {
        VectorSet ns;

        const Fractions& a_c = cell(U).vol_frac;
        for (auto face: cell.faces()) {
            Vector3d S = face.area() * face.normal();

            // На границе возвращает саму ячейку
            const Fractions& a_n = face.neib(U).vol_frac;
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
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        PState z_c = cell(U).get_state();

        // Консервативный вектор в ячейке
        QState q_c(z_c);

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
        q_c.arr() -= (0.5 * dt / cell.volume()) * flux.arr();

        // Значение примитивных переменных на полушаге
        cell(U).half = PState(q_c, mixture, z_c.P(), z_c.T(), z_c.alpha());
    });
}

void MmFluid::fluxes_stage2(Mesh &mesh, double dt, Direction dir)  {
    mesh.for_each([this, dt, dir](Cell &cell) {
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
        for (auto &face: cell.faces(dir)) {
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
        q_c.arr() -= (dt / cell.volume()) * flux.arr();

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
            sum.arr() += QState(child(U).get_state()).arr() * child.volume();
            mean_p += child(U).pressure * child.volume();
            mean_t += child(U).temperature * child.volume();
        }
        sum.arr() /= parent.volume();
        mean_p /= parent.volume();
        mean_t /= parent.volume();
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