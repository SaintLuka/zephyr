#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using namespace zephyr::math;

// ============================================================================
//                               Sod Test
// ============================================================================

SodTest::SodTest() {
    auto eos = IdealGas::create(1.4);
    m_materials += eos;

    rL = 1.0; rR = 0.125;
    pL = 1.0; pR = 0.1;
    uL = 0.0; uR = 0.0;

    x_jump = 0.5;
    finish = 0.2;

    eL = eos->energy_rP(rL, pL);
    eR = eos->energy_rP(rR, pR);

    smf::PState zL(rL, {uL, 0.0, 0.0}, pL, eL);
    smf::PState zR(rR, {uR, 0.0, 0.0}, pR, eR);

    StiffenedGas sg = eos->stiffened_gas(rL, pL);

    exact = RiemannSolver(zL, zR, sg, x_jump);
}

double SodTest::density(const Vector3d& r) const {
    return r.x() < x_jump ? rL : rR;
}

Vector3d SodTest::velocity(const Vector3d& r) const {
    return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
}

double SodTest::pressure(const Vector3d& r) const {
    return r.x() < x_jump ? pL : pR;
}

double SodTest::density_t(const Vector3d &r, double t) const {
    return exact.density(r.x(), t);
}

Vector3d SodTest::velocity_t(const Vector3d &r, double t) const {
    return exact.velocity(r.x(), t) * Vector3d::UnitX();
}

double SodTest::pressure_t(const Vector3d &r, double t) const {
    return exact.pressure(r.x(), t);
}

// ============================================================================
//                           Toro Tests
// ============================================================================

ToroTest::ToroTest(int num, bool multimat, bool adjust_cv) {
    Eos::Ptr eos_L, eos_R;
    if (!multimat) {
        eos_L = IdealGas::create(1.4);
        eos_R = eos_L;
        m_materials += eos_L;
    }
    else {
        eos_L = IdealGas::create(1.3);
        eos_R = IdealGas::create(1.5);
        m_materials += eos_L;
        m_materials += eos_R;
    }

    switch (num) {
        case 1:
            // Тест 1 - это так называемая задача Sod-теста; это очень мягкий тест, и его решение состоит из разрежения слева, контакта и удара справа.
            rL = 1.0;
            uL = 0.75;
            pL = 1.0;
            rR = 0.125;
            uR = 0.0;
            pR = 0.1;
            x_jump = 0.3;
            finish = 0.2;
            break;
        case 2:
            // Тест 2, называемый задачей 123, имеет решение, состоящее из двух сильных разрежений и тривиального стационарного разрыва контакта;
            // давление p∗ очень мало (близко к вакууму), и это может привести к трудностям в итерационной схеме численного нахождения p∗.
            // Тест 2 также полезен при оценке эффективности численных методов для потоков с низкой плотностью.
            rL = 1.0;
            uL = -2.0;
            pL = 0.4;
            rR = 1.0;
            uR = 2.0;
            pR = 0.4;
            x_jump = 0.5;
            finish = 0.15;
            break;
        case 3:
            // Тест 3 - это очень серьезная тестовая задача, решение которой содержит левое разрежение, контакт и удар справа;
            // на самом деле этот тест является левой половиной задачи Вудворда и Колеллы о взрывной волне.
            rL = 1.0;
            uL = 0.0;
            pL = 1000.0;
            rR = 1.0;
            uR = 0.0;
            pR = 0.01;
            x_jump = 0.5;
            finish = 0.012;
            break;
        case 4:
            // Тест 4 представляет собой правую половину задачи Вудворда и Колеллы;
            // его решение содержит удар слева, разрыв контакта и разрежение справа.
            rL = 5.99924;
            uL = 19.5975;
            pL = 460.894;
            rR = 5.99242;
            uR = -6.19633;
            pR = 46.0950;
            x_jump = 0.4;
            finish = 0.035;
            break;
        case 5:
            // Тест 5 состоит из возникающих правого и левого ударов из решения для тестов 3 и 4 соответственно;
            // его решение представляет столкновение этих двух сильных толчков и состоит из удара, направленного влево
            // (очень медленно перемещающегося вправо), разрыва контакта, перемещающегося вправо, и ударной волны, перемещающейся вправо.
            rL = 1.0;
            uL = -19.59745;
            pL = 1000.0;
            rR = 1.0;
            uR = -19.59745;
            pR = 0.01;
            x_jump = 0.8;
            finish = 0.012;
            break;
        case 6:
            rL = 1.4;
            uL = 0.0;
            pL = 1.0;
            rR = 1.0;
            uR = 0.0;
            pR = 1.0;
            x_jump = 0.5;
            finish = 2.0;
            break;
        case 7:
            rL = 1.4;
            uL = 0.1;
            pL = 1.0;
            rR = 1.0;
            uR = 0.1;
            pR = 1.0;
            x_jump = 0.5;
            finish = 2.0;
            break;
        default:
            throw std::runtime_error("Unknown Toro test (num > 7)");
    }

    if (multimat && adjust_cv) {
        eos_L->adjust_cv(rL, pL, 1.0);
        eos_R->adjust_cv(rR, pR, 1.0);
    }

    eL = eos_L->energy_rP(rL, pL);
    eR = eos_R->energy_rP(rR, pR);

    smf::PState zL(rL, {uL, 0.0, 0.0}, pL, eL);
    smf::PState zR(rR, {uR, 0.0, 0.0}, pR, eR);

    StiffenedGas sg_L = eos_L->stiffened_gas(rL, pL);
    StiffenedGas sg_R = eos_R->stiffened_gas(rL, pL);

    exact = RiemannSolver(zL, zR, sg_L, sg_R, x_jump);
}

int ToroTest::index(const Vector3d& r) const {
    return m_materials.single() ? 0 : (r.x() < x_jump ? 0 : 1);
}

double ToroTest::density(const Vector3d& r) const {
    return r.x() < x_jump ? rL : rR;
}

Vector3d ToroTest::velocity(const Vector3d& r) const {
    return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
}

double ToroTest::pressure(const Vector3d& r) const {
    return r.x() < x_jump ? pL : pR;
}

int ToroTest::index_t(const Vector3d& r, double t) const {
    if (m_materials.single()) {
        return 0;
    } else {
        bool left = exact.fraction(r.x(), t) > 0.5;
        return left ? 0 : 1;
    }
}

double ToroTest::density_t(const Vector3d &r, double t) const {
    return exact.density(r.x(), t);
}

Vector3d ToroTest::velocity_t(const Vector3d &r, double t) const {
    return exact.velocity(r.x(), t) * Vector3d::UnitX();
}

double ToroTest::pressure_t(const Vector3d &r, double t) const {
    return exact.pressure(r.x(), t);
}

// ============================================================================
//                               Shock Wave
// ============================================================================

ShockWave::ShockWave(double Ms, double x_jump,  double length, Params params)
    : x_jump(x_jump), length(length) {

    bool inv = Ms < 0.0;
    Ms = std::abs(Ms);

    double gamma = params.gamma;
    auto gas = IdealGas::create(gamma);
    m_materials += gas;

    pR = params.P;
    rR = params.rho;
    uR = 0.0;

    pL = pR * (2 * gamma * Ms * Ms - gamma + 1) / (gamma + 1);
    rL = rR * (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms);
    uL = 2 / Ms * std::sqrt(gamma * pR / rR) * (Ms * Ms - 1) / (gamma + 1);

    speed = Ms * std::sqrt(gamma * pR / rR);
    finish = (length - x_jump) / speed;

    if (inv) {
        uL *= -1.0;
        std::swap(rL, rR);
        std::swap(uL, uR);
        std::swap(pL, pR);
    }
}

double ShockWave::density(const Vector3d& r) const {
    return r.x() < x_jump ? rL : rR;
}

Vector3d ShockWave::velocity(const Vector3d& r) const {
    return {r.x() < x_jump ? uL : uR, 0, 0};
}

double ShockWave::pressure(const Vector3d& r) const {
    return r.x() < x_jump ? pL : pR;
}

double ShockWave::density_t(const Vector3d& r, double t) const {
    return r.x() < x_jump + speed * t ? rL : rR;
}

Vector3d ShockWave::velocity_t(const Vector3d& r, double t) const {
    return {r.x() < x_jump + speed * t ? uL : uR, 0, 0};
}

double ShockWave::pressure_t(const Vector3d& r, double t) const {
    return r.x() < x_jump + speed * t ? pL : pR;
}

// ============================================================================
//                               Shu-Osher Test
// ============================================================================


ShuOsherTest::ShuOsherTest() {
    m_materials += IdealGas::create(1.4);

    rL = 27.0 / 7.0;
    uL = 4.0 * std::sqrt(35.0) / 9.0;
    pL = 31.0 / 3.0;

    rR  = 1.0;
    uR = 0.0;
    pR = 1.0;

    epsilon = 0.2;
}

double ShuOsherTest::density(const Vector3d& r) const {
    return r.x() < x_jump ? rL : rR + epsilon * std::sin(5.0 * r.x());
}

Vector3d ShuOsherTest::velocity(const Vector3d& r) const {
    return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
}

double ShuOsherTest::pressure(const Vector3d& r) const {
    return r.x() < x_jump ? pL : pR;;
}

// ============================================================================
//                               Rarefied Water
// ============================================================================

/// @brief Конструктор
RarefiedWater::RarefiedWater() {
    auto water = StiffenedGas::create("Water");
    auto air   = IdealGas::create("Air");
    m_materials += water;
    m_materials += air;

    double T = 0.0_C;

    rL = 0.9_g_cm3;
    uL = 0.0;
    pL = water->pressure_rT(rL, T);
    eL = water->energy_rP(rL, pL);

    rR = 1.16_kg_m3;
    uR = 0.0;
    pR = air->pressure_rT(rR, T);
    eR = air->energy_rP(rR, pR);

    x_jump = 0.9_cm;
    finish = 5.0_us;
}

int RarefiedWater::index(const Vector3d& r) const {
    return r.x() < x_jump ? 0 : 1;
}


double RarefiedWater::density(const Vector3d& r) const {
    return r.x() < x_jump ? rL : rR;
}

Vector3d RarefiedWater::velocity(const Vector3d& r) const {
    return {r.x() < x_jump ? uL : uR, 0, 0};
}

double RarefiedWater::pressure(const Vector3d& r) const {
    return r.x() < x_jump ? pL : pR;
}

// ============================================================================
//                         Multimat 1D Tests
// ============================================================================


std::tuple<std::string, std::string> matpair(int mat) {
    switch (mat) {
        case 1:  return {"Gas"     , "Gas"     };
        case 2:  return {"HeavyGas", "HeavyGas"};
        case 3:  return {"Liquid"  , "Liquid"  };
        case 4:  return {"Solid"   , "Solid"   };
        case 5:  return {"Gas"     , "HeavyGas"};
        case 6:  return {"Gas"     , "Liquid"  };
        case 7:  return {"Gas"     , "Solid"   };
        case 8:  return {"HeavyGas", "Liquid"  };
        case 9:  return {"HeavyGas", "Solid"   };
        case 10: return {"Liquid"  , "Solid"   };
        default:
            throw std::runtime_error("Unknown materials pair, number " + std::to_string(mat) + ".");
    }
}

Multimat1D::Multimat1D(int num, int mat, int c) {
    // В газодинамике все гоняют одномерные безразмерные тесты.
    // У нас тоже будут практически «безразмерные» тесты. Но для наглядности
    // можно считать, что давление измеряется в атмосферах (10^5 Па),
    // плотности в кг/м^3, скорости в 100√10 ≈ 320 м/c, время в миллисекундах,
    // длина в 1/√10 ≈ 0.32 м, температура кратна 300 К.
    // Выберем несколько материалов для тестов. Характеристики заданы в наших
    // безразмерных переменных. Для газов (P = 1, T = 1, ρ = ρ0).
    // Для жидкости и металла P = 0, T = 1, ρ = ρ0. В общем, ρ0 в таблице это
    // характерная плотность материала, при нормальных безразмерных условиях P = T = 1.
    //
    // Материалы в таблице вдохновлены следующими материалами. Газ — простой воздух,
    // тяжелый газ — фторид серы (SF6 ), жидкость — вода, металл — свинец.
    // Соотношения безразмерных параметров у разных материалов близко к
    // естественным. К примеру, видно, что вода имеет теплоемкость гораздо больше свинца.
    //
    // Скорость звука в тестовых материалах:
    //   Gas      ~ 1.0
    //   HeavyGas ~ 0.5
    //   Liquid   ~ 17.0
    //   Solid    ~ 7.0

    auto[name1, name2] = matpair(mat);

    eos_L = StiffenedGas::create(name1);
    eos_R = StiffenedGas::create(name2);

    width = 0.0;

    if (num == 0) {
        // num: 0; m = 1..10;
        //   c = 0: плотность справа= rho0;
        //   c = 1: плотность справа < rho0;
        //   c = 2: плотность справа > rho0;

        // Простой дозвуковой перенос
        x_cont = 0.2;
        finish = 1.0;
        double speed = 0.5;

        uL1 = uL2 = uR1 = uR2 = speed;
        pL1 = pL2 = pR1 = pR2 = 1.0;
        if (name1 == "Gas") {
            rL1 = rL2 = 1.0;
        }
        else if (name1 == "HeavyGas") {
            rL1 = rL2 = 4.0;
        }
        else if (name1 == "Liquid") {
            rL1 = rL2 = 1000.0;
        }
        else if (name1 == "Solid") {
            rL1 = rL2 = 12500.0;
        }

        if (name2 == "Gas") {
            rR1 = rR2 = c < 1 ? 1.0 : (c < 2 ? 0.5 : 2.0);
        }
        else if (name2 == "HeavyGas") {
            rR1 = rR2 = c < 1 ? 4.0 : (c < 2 ? 1.0 : 6.0);
        }
        else if (name2 == "Liquid") {
            rR1 = rR2 = c < 1 ? 1000.0 : (c < 2 ? 800.0 : 1200.0);
        }
        else if (name2 == "Solid") {
            rR1 = rR2 = c < 1 ? 12500.0 : (c < 2 ? 11000 : 14000);
        }
    }
    else if (num == 1) {
        // num: 0; m = 1..10;
        //   c = 0: x_cont = x_jump;
        //   c = 1: x_cont < x_jump;
        //   c = 2: x_cont > x_jump;

        // Аналог теста Сода
        x_jump = 0.3;
        finish = 0.2;

        x_cont = c == 0 ? 0.3 : (c == 1 ? 0.2 : 0.5);

        // Один материал
        if (name1 == name2) {
            rL1 = rL2 = 1.0;
            uL1 = uL2 = 0.75;
            pL1 = pL2 = 1.0;

            rR1 = rR2 = 0.125;
            uR1 = uR2 = 0.0;
            pR1 = pR2 = 0.1;

            if (name1 == "HeavyGas") {
                rL1 = rL2 = 2.0;
                pL1 = pL2 = 1.5;
            }

            if (name1 == "Liquid") {
                rL1 = rL2 = 1000.0;
                uL1 = uL2 = 10.75;
                pL1 = pL2 = 0.0;

                rR1 = rR2 = 9.0;
                uR1 = uR2 = 0.0;
                pR1 = pR2 = -5;
            }


            if (c == 1) {
                rR1 = rL2;
                uR1 = uL2;
                pR1 = pL2;
            }
            if (c == 2) {
                rL2 = rR1;
                uL2 = uR1;
                pL2 = pR1;
            }
        }
    }

    eL1 = eos_L->energy_rP(rL1, pL1);
    eL2 = eos_L->energy_rP(rL2, pL2);

    eR1 = eos_R->energy_rP(rR1, pR1);
    eR2 = eos_R->energy_rP(rR2, pR2);

    tL1 = eos_L->temperature_rP(rL1, pL1);
    tL2 = eos_L->temperature_rP(rL2, pL2);

    tR1 = eos_R->temperature_rP(rR1, pR1);
    tR2 = eos_R->temperature_rP(rR2, pR2);

    m_materials += eos_L;
    m_materials += eos_R;
}

double Multimat1D::density(const Vector3d &r) const {
    double P = pressure(r);
    double T = temperature(r);
    int idx = index(r);
    return 1.0 / m_materials[idx].volume_PT(P, T);
}

Vector3d Multimat1D::velocity(const Vector3d &r) const {
    double x = r.x();
    double u = x < x_jump ? (x < x_cont ? uL1 : uR1) : (x < x_cont ? uL2 : uR2);
    return {u, 0.0, 0.0};
}

double Multimat1D::pressure(const Vector3d &r) const {
    double x = r.x();
    return x < x_jump ? (x < x_cont ? pL1 : pR1) : (x < x_cont ? pL2 : pR2);
}

double Multimat1D::temperature(const Vector3d &r) const {
    double x = r.x();

    double T1 = tL2;
    double T2 = tR1;

    if (x <= x_cont - width) {
        return x < x_jump ? tL1 : tL2;
    }
    else if (x > x_cont + width) {
        return x < x_jump ? tR1 : tR2;
    }

    return 0.5 * ((T1 + T2) + (T2 - T1) * (std::sin(0.5 * M_PI * (x - x_cont) / width)));
}

double Multimat1D::energy(const Vector3d &r) const {
    double P = pressure(r);
    double T = temperature(r);
    return m_materials[index(r)].energy_PT(P, T);
}

} // namespace zephyr::phys