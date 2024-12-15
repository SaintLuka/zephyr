#include <zephyr/phys/tests/multimat_1D.h>
#include <iostream>

namespace zephyr::phys {

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
}

double Multimat1D::density(const Vector3d &r) const {
    double P = pressure(r);
    double T = temperature(r);
    return 1.0 / get_eos(r)->volume_PT(P, T);
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
    return get_eos(r)->energy_PT(P, T);
}

Eos::Ptr Multimat1D::get_eos(const Vector3d &r) const {
    return r.x() < x_cont ? eos_L : eos_R;
}

double Multimat1D::fraction(const Vector3d& r, int mat) const {
    return r.x() < x_cont ? (mat == 0 ? 1.0 : 0.0) : (mat == 1 ? 1.0 : 0.0);
}

} // namespace zephyr::phys