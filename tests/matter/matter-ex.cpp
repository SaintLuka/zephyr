// @brief Пример работы с моделями материалов

#include <iostream>

#include <zephyr/phys/matter/material.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>

#include <zephyr/phys/matter/materials.h>

using namespace zephyr::phys;

int main() {
    Material air = IdealGas::create("Air");
    air += Conductivity::create("Air");

    // Все функции свойств материала вызываются напрямую

    // Функция уравнения состояния
    double P0 = air.pressure_re(1.0, 1.0);

    // Многие функции УРС возвращают тройки: значения + производные,
    // для этого требуется передать опцию {.deriv = true}
    auto[P, dP_dr, dP_de] = air.pressure_re(1.0, 1.0, {.deriv=true});

    // Теплоемкость при постоянном объеме
    double Cv = air.energy_rT(1.0, 300, {.deriv=true}).dT;

    // Коэффициент теплопроводности
    double kappa = air.kappa();

    // Остальные параметры соответствуют сжимаемой жидкости/газу,
    // нет необходимости добавлять отдельно
    double visc = air.shear_visc(); // = 0.0

    // Также нет необходимости добавлять Plasticity для газа
    double Y  = air.yield();        // = 0.0
    double nu = air.poisson();      // = 0.5

    std::cout << "Pressure:        " << P0 << "\n";
    std::cout << "Pressure:       (" << P << ", " << dP_dr << ", " << dP_de << ")\n";
    std::cout << "Heat capacity:   " << Cv << "\n";
    std::cout << "Conductivity:    " << kappa << "\n";
    std::cout << "Viscosity:       " << visc << "\n";
    std::cout << "Yield stress:    " << Y << "\n";
    std::cout << "Poisson's ratio: " << nu << "\n\n";


    // Создадим вязкую воду
    Material water = StiffenedGas::create("Water");
    water += Viscosity::create("Water");

    std::cout << "Water visc:   " << water.shear_visc() << "\n";

    // Создадим пластичную медь с теплопроводностью
    Material cu = MieGruneisen::create("Cu");
    cu += Conductivity::create("Cu");
    cu += Plasticity::create("Cu");

    std::cout << "Copper yield: " << cu.yield() << "\n\n";

    // Создаем список материалов, можно сразу перечислить всё
    Materials materials = {air, water};

    // или добавлять новые материалы в конец
    materials += cu;

    std::cout << "Air \"density\":  " << materials[0].density() << "\n";
    std::cout << "Water density:  "   << materials[1].density() << "\n";
    std::cout << "Copper density: "   << materials[2].density() << "\n\n";


    // Продвинутые свойства: доступ к компонентам материала

    // Можем вытащить из материала уравнение состояния
    Eos::Ptr air_eos = air.eos();

    std::cout << "Old pressure: " << air_eos->pressure_re(1.0, 1.0) << "\n";

    // Eos -- абстрактный класс, нет прямого обращения к данным
    // air_eos->gamma;

    // Но можно сделать приведение типа (меняем показатель адиабаты)
    // Хотя при замене параметров УрС не гарантируется совместность функций.
    air_eos->cast<IdealGas>().gamma = 2.0;

    std::cout << "New pressure: " << air_eos->pressure_re(1.0, 1.0) << "\n";

    return 0;
}