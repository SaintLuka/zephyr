#pragma once

#include <zephyr/geom/surface/solid_body.h>

namespace zephyr::geom {

/// @brief Некоторая геометрия
class SolidBody3D : public SolidBody {
public:
    /// @brief Базовый конструктор
    SolidBody3D();

    /// @brief Простейший конструктор с указанием положения центра
    explicit SolidBody3D(const Vector3d& center);

    /// @brief Виртуальный деструктор
    ~SolidBody3D() override = default;

    /// @brief Поворот в локальной системе координат
    /// @param axis Ось вращения
    /// @param angle Угол поворота против часовой стрелки
    void rotation_local(const Vector3d& axis, double angle);

    /// @brief Поворот в глобальной системе координат
    /// @param axis Ось вращения
    /// @param angle Угол поворота против часовой стрелки
    void rotation_global(const Vector3d& axis, double angle);

    /// @brief Поворот относительно точки
    /// @param axis Ось вращения
    /// @param angle Угол поворота против часовой стрелки
    /// @param c Точка, вокруг которой выполняется вращение
    void rotation_relative(const Vector3d& axis, double angle, const Vector3d& c);

    /// @brief Триангуляция поверхности
    /// @param n_elements Максимальное число треугольных элементов
    std::vector<std::array<Vector3d, 3>> triangulation(int n_elements) const;

    /// @brief Объемная доля тела внутри ячейки сетки (наследуется)
    using SolidBody::volume_fraction;

    /// @brief Объем тела внутри ячейки (наследуется)
    using SolidBody::volume_inside;

protected:
    /// @brief Триангуляция поверхности в локальной системе координат
    /// @param n_elements Максимальное число треугольных элементов
    virtual std::vector<std::array<Vector3d, 3>> local_triangulation(int n_elements) const = 0;
};

/// @brief Шар
class BodyBall : public SolidBody3D {
public:
    /// @brief Шар в начале координат
    /// @param radius Радиус шара
    explicit BodyBall(double radius);

    /// @brief Шар с указанием центра
    /// @param radius Радиус шара
    /// @param center Центр шара
    BodyBall(double radius, const Vector3d& center);

protected:
    /// @brief Инициализация m_inside
    void init_inside();

    /// @brief Триангуляция поверхности в локальной системе координат
    /// @param n_elements Максимальное число треугольных элементов
    std::vector<std::array<Vector3d, 3>> local_triangulation(int n_elements) const override;

    double m_radius;  ///< Радиус шара
};

/// @brief Куб
class BodyCube : public SolidBody3D {
public:
    /// @brief Куб в центре координат
    /// @param length Длина стороны куба
    explicit BodyCube(double length);

    /// @brief Куб с центром в определенной точке
    /// @param length Длина стороны куба
    /// @param center Центр куба
    BodyCube(double length, const Vector3d& center);

protected:
    /// @brief Инициализация m_inside
    void init_inside();

    /// @brief Триангуляция поверхности в локальной системе координат
    /// @param n_elements Максимальное число треугольных элементов (не используется)
    std::vector<std::array<Vector3d, 3>> local_triangulation(int n_elements) const override;

    double m_length; ///< Длина стороны
};

} // namespace zephyr::geom