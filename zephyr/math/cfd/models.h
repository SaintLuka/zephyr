#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/math/vectorization.h>
#include <zephyr/math/cfd/rotate.h>
#include <zephyr/phys/eos/eos.h>
#include <ostream>

namespace zephyr { namespace math {

using namespace geom;

/// @brief Одноматериальная модель Single Material Fluid
namespace smf {

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

    friend std::ostream &operator<<(std::ostream &os, const PState &state);

    VECTORIZE(PState)
};

/// @brief Консервативный вектор состояния
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

/**
 * @brief Вектор потока
 * @code
 *  mass = rho * u;
 *  momentum.x() = rho * u^2 + p;
 *  momentum.y() = rho * u * v;
 *  momentum.z() = rho * u * w;
 *  energy = u * (E + p);
 * @endcode
 */
struct Flux {
    double mass; ///< rho * u
    Vector3d momentum;
    double energy; ///< u * (rho * (E + 0.5 * velocity^2) + p)

    /// @brief Нулевой поток
    Flux();

    Flux(double mass, const Vector3d &momentum, double energy);

    /// @brief Дифференциальный поток по вектору примитивных переменных
    Flux(const PState& state);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d& normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    Flux in_local(const Vector3d& normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d& normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    Flux in_global(const Vector3d& normal) const;

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux);

    VECTORIZE(Flux)
};

} // namespace smf

} // namespace math
} // namespace zephyr