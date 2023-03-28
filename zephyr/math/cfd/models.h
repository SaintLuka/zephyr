#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/math/vectorization.h>
#include <zephyr/math/cfd/rotate.h>
#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace math {

using namespace geom;

struct QState;

/// @brief Примитивный вектор состояния
struct PState {
    double density;
    Vector3d velocity;
    double pressure;
    double energy;

    PState(const double &density, const Vector3d &velocity,
           const double &pressure, const double &energy);

    PState(const QState& q, const phys::Eos& eos);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d& normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    PState in_local(const Vector3d& normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d& normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    PState in_global(const Vector3d& normal) const;

    VECTORIZE(PState)
};

/// @brief Консеравтивный вектор состояния
struct QState {
    double mass;
    geom::Vector3d momentum;
    double energy;

    QState(const double &mass, const Vector3d &momentum, const double &energy);

    explicit QState(const PState& z);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d& normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    QState in_local(const Vector3d& normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d& normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    QState in_global(const Vector3d& normal) const;

    VECTORIZE(QState)
};

/// @brief Вектор потока
struct Flux {
    double mass;
    Vector3d momentum;
    double energy;

    /// @brief Нулевой поток
    Flux();

    /// @brief Дифференциальный поток по вектору примитивных переменных
    Flux(const PState& z);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d& normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    Flux in_local(const Vector3d& normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d& normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    Flux in_global(const Vector3d& normal) const;

    VECTORIZE(Flux)
};

}
}