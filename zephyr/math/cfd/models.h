#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/math/vectorization.h>
#include <zephyr/math/cfd/rotate.h>
#include <zephyr/phys/eos/eos.h>

#include <zephyr/phys/fractions.h>

namespace zephyr::math {

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

    PState(const QState &q, const phys::Eos &eos);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    [[nodiscard]] PState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    [[nodiscard]] PState in_global(const Vector3d &normal) const;

    friend std::ostream &operator<<(std::ostream &os, const PState &state);

    VECTORIZE(PState)
};

/// @brief Консервативный вектор состояния
struct QState {
    double mass;
    geom::Vector3d momentum;
    double energy;

    QState(const double &mass, const Vector3d &momentum, const double &energy);

    explicit QState(const PState &z);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    [[nodiscard]] QState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    [[nodiscard]] QState in_global(const Vector3d &normal) const;

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
    explicit Flux(const PState &state);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    [[nodiscard]] Flux in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    [[nodiscard]] Flux in_global(const Vector3d &normal) const;

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux);

    VECTORIZE(Flux)
};

} // namespace smf

/// @brief Многоматериальная модель Multi Material Fluid
namespace mmf {

using zephyr::phys::Fractions;

struct Component {
    double density;
    double energy;
    double frac;
};

struct QState;

/// @brief Примитивный многоматериальный вектор состояния
struct PState {
    double density;
    Vector3d velocity;
    double pressure;
    double temperature;
    double energy;
    Fractions mass_frac;

    PState(const double &pressure, const double &temperature, const Vector3d &velocity,
           const std::vector<Component> &components);

    PState(const double &density, const Vector3d &velocity,
           const double &pressure, const double &energy, const double &temperature, const Fractions &mass_frac);

    PState(const QState &q, const phys::Eos &eos);

    [[nodiscard]] std::vector<double> get_densities() const;

    [[nodiscard]] std::vector<double> get_energies() const;

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    [[nodiscard]] PState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    [[nodiscard]] smf::PState to_smf() const;

    /// @brief Возвращает вектор состояния в глобальной системе координат
    [[nodiscard]] PState in_global(const Vector3d &normal) const;

    friend std::ostream &operator<<(std::ostream &os, const PState &state);

    VECTORIZE(PState)
};

/// @brief Консервативный многоматериальный вектор состояния
struct QState {
    double mass;
    Vector3d momentum;
    double energy;
    Fractions mass_frac;

    QState(const double &mass, const Vector3d &momentum, const double &energy, const Fractions &mass_frac);

    explicit QState(const PState &state);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    [[nodiscard]] QState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    [[nodiscard]] QState in_global(const Vector3d &normal) const;

    VECTORIZE(QState)
};

/**
 * @brief Многоматериальный вектор потока
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
    Fractions mass_frac;

    /// @brief Нулевой поток
    Flux();

    Flux(double mass, const Vector3d &momentum, double energy, const Fractions &mass_frac);

    /// @brief Дифференциальный поток по вектору примитивных переменных
    explicit Flux(const PState &state);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    [[nodiscard]] Flux in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    [[nodiscard]] Flux in_global(const Vector3d &normal) const;

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux);

    VECTORIZE(Flux)
};

} // namespace mmf

} // namespace zephyr