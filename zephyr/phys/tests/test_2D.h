#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/boundary.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/tests/ivp.h>

namespace zephyr::phys {

/// @class Абстрактный класс для двумерных тестов.
/// Только в прямоугольнике?
class Test2D : public IVP {
public:
    // Границы области

    virtual double xmin() const = 0;

    virtual double xmax() const = 0;

    virtual double ymin() const = 0;

    virtual double ymax() const = 0;

    /// @brief Граничные условия для прямоугольника
    virtual Rectangle::Boundaries boundaries() const {
        return {.left=Boundary::ZOE, .right=Boundary::ZOE,
                .bottom=Boundary::ZOE, .top=Boundary::ZOE};
    }

    /// @brief Получить название теста
    std::string name() const override { return "Test2D"; }
};

/// @brief Двумерный распад разрыва.
/// Набор тестов взят из статьи?
class Riemann2D : public Test2D {
public:
    double x_jump;      ///< Положение разрыва
    double y_jump;      ///< Положение разрыва

    double finish;      ///< Конечный момент времени
    
    double r1, r2, r3, r4;  ///< Плотность
    double u1, u2, u3, u4;  ///< Скорость x
    double v1, v2, v3, v4;  ///< Скорость y
    double p1, p2, p3, p4;  ///< Давление


    /// @brief Конструктор
    explicit Riemann2D(int test_case);

    /// @brief Название теста
    std::string name() const final { return "2D Riemann Test"; }

    // Границы области

    double xmin() const final { return 0.0; }

    double xmax() const final { return 1.0; };

    double ymin() const final { return 0.0; };

    double ymax() const final { return 1.0; };

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    double density(const Vector3d &vec) const final;

    Vector3d velocity(const Vector3d &vec) const final;

    double pressure(const Vector3d &vec) const final;
};

/// @brief Неустойчивость Рихтмайера -- Мешкова,
/// пока что одноматериальная, постановка придумана мной
class RichtmyerMeshkov : public Test2D {
public:
    double finish;       ///< Конечный момент времени
    double x_contact;    ///< Положение контакта
    double x_jump;       ///< Положение ударной волны
    double y_width;      ///< Толщина трубки
    double amplitude;    ///< Амплитуда возмущения
    int    n_peaks;      ///< Число возмущений

    double pL, rL, uL;  ///< За фронтом УВ
    double pR, rR, uR;  ///< Перед фронтом УВ
    double p0, r0, u0;  ///< Область за контактом


    /// @brief Конструктор
    /// @param Ms Число Маха
    explicit RichtmyerMeshkov(double Ms = 3.0, bool multimat = false);

    /// @brief Название теста
    std::string name() const final {
        return "Richtmyer-Meshkov instability";
    }

    // Границы области

    double xmin() const final { return 0.0; }

    double xmax() const final { return 10.0; }

    double ymin() const final { return -0.5 * y_width; }

    double ymax() const final { return +0.5 * y_width; }

    Rectangle::Boundaries boundaries() const final {
        return {.left=Boundary::ZOE, .right=Boundary::WALL,
                .bottom=Boundary::WALL, .top=Boundary::WALL};
    }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные условия

    /// @brief Номер области
    /// 0 -- за фронтом ударной волны
    /// 1 -- между УВ и контактом
    /// 2 -- Невозмущенная область справа от контакта
    int region(const Vector3d& r) const;

    int index(const Vector3d& r) const final;

    double density(const Vector3d &r) const final;

    Vector3d velocity(const Vector3d &r) const final;

    double pressure(const Vector3d &r) const final;
};

} // namespace zephyr::phys