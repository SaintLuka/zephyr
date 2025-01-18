#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/math/vectorization.h>
#include <zephyr/math/cfd/rotate.h>

#include <zephyr/phys/fractions.h>
#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/phys/matter/mixture_pt.h>

namespace zephyr::math {

using namespace geom;

/// @brief Одноматериальная модель Single Material Fluid
namespace smf {

struct QState;

/// @brief Примитивный вектор состояния
struct PState {
    double   density;  ///< Плотность
    Vector3d velocity; ///< Скорость
    double   pressure; ///< Давление
    double   energy;   ///< Удельная внутренняя энергия (доп)

    /// @brief Инициализация нулями
    PState();

    /// @brief Инициализация с полным заданием параметров
    PState(const double &density, const Vector3d &velocity,
           const double &pressure, const double &energy);

    /// @brief Инициализация из консервативного вектора состояния,
    /// давление определяется с использованием УрС.
    PState(const QState &q, const phys::Eos &eos);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    PState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    PState in_global(const Vector3d &normal) const;

    /// @brief Отражает систему координат
    void inverse();

    // Короткий доступ к полям
    inline const double& rho() const { return density; }
    inline const double& u() const { return velocity.x(); }
    inline const double& v() const { return velocity.y(); }
    inline const double& w() const { return velocity.z(); }
    inline const double& P() const { return pressure; }
    inline const double& e() const { return energy; }

    // Квадрат модуля скорости
    inline double v2() const { return velocity.squaredNorm(); }

    // Полная удельная энергия
    inline double E() const { return energy + 0.5 * v2(); }

    // Проверить корректность
    bool is_bad(const phys::Eos &eos);

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const PState &state);

    VECTORIZE(PState)
};

/// @brief Консервативный вектор состояния
struct QState {
    double   density;   ///< Плотность
    Vector3d momentum;  ///< Плотность импульса:       rho * v
    double   energy;    ///< Плотность полной энергии: rho * (e + 0.5 * v^2)

    /// @brief Инициализация нулями
    QState();

    /// @brief Инициализация с полным заданием параметров
    QState(const double &mass, const Vector3d &momentum, const double &energy);

    /// @brief Преобразование из примитивного вектора состояния
    explicit QState(const PState &z);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    QState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    QState in_global(const Vector3d &normal) const;

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const QState &state);

    VECTORIZE(QState)
};

/// @brief Вектор потока
/// @code
///   mass = rho * u;
///   momentum.x() = rho * u^2 + P;
///   momentum.y() = rho * u * v;
///   momentum.z() = rho * u * w;
///   energy = u * (rho * (e + 0.5 * velocity^2) + P);
/// @endcode
struct Flux {
    double   mass;      ///< Плотность потока массы
    Vector3d momentum;  ///< Плотность потока импульса
    double   energy;    ///< Плотность потока энергии

    /// @brief Нулевой поток
    Flux();

    /// @brief Инициализация с полным заданием параметров
    Flux(double mass, const Vector3d &momentum, double energy);

    /// @brief Дифференциальный поток по вектору примитивных переменных
    explicit Flux(const PState &z);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    Flux in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    Flux in_global(const Vector3d &normal) const;

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const Flux &flux);

    VECTORIZE(Flux)
};

} // namespace smf

/// @brief Многоматериальная модель Multi Material Fluid
namespace mmf {

using zephyr::phys::Fractions;
using zephyr::phys::ScalarSet;
using zephyr::phys::MixturePT;

struct QState;

/// @brief Примитивный многоматериальный вектор состояния
/// Используется в модели равновесной по скорости, давлению
/// и температуре (P-V-T модель).
struct PState {
    double    density;      ///< Плотность смеси
    Vector3d  velocity;     ///< Равновесная корость
    double    pressure;     ///< Равновесное давление
    double    energy;       ///< Внутренняя энергия смеси (доп)
    double    temperature;  ///< Равновесная температура (доп)
    Fractions mass_frac;    ///< Массовые доли компонент
    Fractions vol_frac;     ///< Объемные доли компонент (доп)

    /// @brief Инициализация нулями
    PState();

    /// @brief Инициализация с полным заданием параметров
    PState(double density, const Vector3d &velocity, double pressure, double energy,
           double temperature, const Fractions &mass_frac, const Fractions& vol_frac);

    /// @brief Инициализация из консервативного вектора состояния,
    /// давление, температура и объемные доли определяются из уравнения
    /// состояния смеси (PT - замыкание).
    /// В качестве начальных приближений в следует передавать начальные
    /// приближения для давления (P0), температуры (T0) и объемных долей (alpha).
    PState(const QState &q, const MixturePT &mixture,
           double P0, double T0, const Fractions& alpha);


    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    PState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    PState in_global(const Vector3d &normal) const;

    /// @brief Отражает систему координат
    void inverse();

    // Короткий доступ к полям
    inline const double& rho() const { return density; }
    inline const double& u() const { return velocity.x(); }
    inline const double& v() const { return velocity.y(); }
    inline const double& w() const { return velocity.z(); }
    inline const double& P() const { return pressure; }
    inline const double& T() const { return temperature; }
    inline const double& e() const { return energy; }
    inline const Fractions& alpha() const { return vol_frac; }
    inline const Fractions& beta() const { return mass_frac; }

    /// @brief Плотность компоненты
    double true_density(int idx) const;

    /// @brief Энергия компоненты energy(P, T)
    double true_energy(const MixturePT& mixture, int idx) const;

    /// @brief Преобразование в одноматериальный вектор состояния
    smf::PState to_smf() const;

    /// @brief Выделить чистую компоненту материала с индексом idx
    smf::PState extract(const MixturePT& mixture, int idx) const;

    /// @brief Выделяет из вектора состояния компоненту с материалом с индексом idx
    /// Если вектор не содержит материал с индексом idx, или содержит только его,
    /// то поведение не определено, там всякие NAN'ы будут
    std::pair<mmf::PState, mmf::PState> split(const MixturePT& mixture, int iA) const;

    /// @brief Обновить значения, полученные путем интерполяции.
    /// Восстанавливает согласованность состояния (термодинамических величин).
    /// Нормализует концентрации, выставляет согласованные значения
    /// для энергии и температуры (вычисляются по плотности и давлению)
    void interpolation_update(const phys::MixturePT& mixture);

    bool is_bad() const;

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const PState &state);

    VECTORIZE(PState)
};

/// @brief Консервативный многоматериальный вектор состояния
struct QState {
    double    density;    ///< Плотность
    Vector3d  momentum;   ///< Плотность импульса:       rho * v
    double    energy;     ///< Плотность полной энергии: rho * (e + 0.5 * v^2)
    ScalarSet mass_frac;  ///< Плотности компонент:      rho * beta_i

    /// @brief Инициализация нулями
    QState();

    /// @brief Инициализация с полным заданием параметров
    QState(double density, const Vector3d &momentum, double energy, const ScalarSet &mass_frac);

    /// @brief Преобразование из примитивного вектора состояния
    explicit QState(const PState &z);

    /// @brief Переводит вектор состояния в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в локальной системе координат
    QState in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор состояния в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор состояния в глобальной системе координат
    QState in_global(const Vector3d &normal) const;

    friend std::ostream &operator<<(std::ostream &os, const QState &state);

    VECTORIZE(QState)
};

/// @brief Многоматериальный вектор потока
/// @code
///   mass         = rho * u;
///   momentum.x   = rho * u^2 + P;
///   momentum.y   = rho * u * v;
///   momentum.z   = rho * u * w;
///   energy       = u * (rho * E + P);
///   mass_frac[i] = rho * u * beta[i];
/// @endcode
struct Flux {
    double    density;    ///< rho * u
    Vector3d  momentum;
    double    energy;     ///< u * (rho * (e + 0.5 * v^2) + P)
    ScalarSet mass_frac;  ///< rho * u * beta_i

    /// @brief Нулевой поток
    Flux();

    Flux(double mass, const Vector3d &momentum, double energy, const ScalarSet &mass_frac);

    /// @brief Дифференциальный поток по вектору примитивных переменных
    explicit Flux(const PState &z);

    /// @brief Преобразует одноматериальный поток материала с индексом mat
    /// в многоматериальный
    Flux(const smf::Flux& flux, int mat);

    /// @brief Переводит вектор потока в локальную систему координат
    void to_local(const Vector3d &normal);

    /// @brief Возвращает вектор потока в локальной системе координат
    Flux in_local(const Vector3d &normal) const;

    /// @brief Переводит вектор потока в глобальную систему координат
    void to_global(const Vector3d &normal);

    /// @brief Возвращает вектор потока в глобальной системе координат
    Flux in_global(const Vector3d &normal) const;

    /// @brief Отражает систему координат
    void inverse();

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux);

    bool is_bad() const;

    VECTORIZE(Flux)
};

} // namespace mmf

} // namespace zephyr::math