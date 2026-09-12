#pragma once

#include <zephyr/geom/surface/solid_body.h>

namespace zephyr::geom {

/// @brief Некоторая геометрия
class SolidBody2D : public SolidBody {
public:
    /// @brief Базовый конструктор
    SolidBody2D();

    /// @brief Простейший конструктор с указанием положения центра
    explicit SolidBody2D(const Vector3d& center);

    /// @brief Виртуальный деструктор
    ~SolidBody2D() override = default;

    /// @brief Поворот в локальной системе координат
    /// @param angle Угол поворота против часовой стрелки
    void rotation_local(double angle);

    /// @brief Поворот в глобальной системе координат
    /// @param angle Угол поворота против часовой стрелки
    void rotation_global(double angle);

    /// @brief Поворот относительно точки
    /// @param angle Угол поворота против часовой стрелки
    /// @param c Точка, вокруг которой выполняется вращение
    void rotation_relative(double angle, const Vector3d& c);

    /// @brief Многоугольник, формирующий границы тела
    /// @param n_points Максимальное число точек
    std::vector<Vector3d> outline(int n_points) const;

    /// @brief Объемная доля тела внутри ячейки сетки (наследуется)
    using SolidBody::volume_fraction;

    /// @brief Объем тела внутри ячейки (наследуется)
    using SolidBody::volume_inside;

protected:
    /// @brief Многоугольник, формирующий границы тела в локальной системе координат
    /// @param n_points Максимальное число точек
    virtual std::vector<Vector3d> local_outline(int n_points) const = 0;
};

/// @brief Тело в виде бесконечной вертикальной полосы
class BodyStrip : public SolidBody2D {
public:
    /// @brief Полоса в центре координат
    /// @param width Ширина полосы
    explicit BodyStrip(double width);

    /// @brief Полоса с центром в определенной точке
    /// @param width Ширина полосы
    /// @param center Центр квадрата
    BodyStrip(double width, const Vector3d& center);

protected:
    /// @brief Инициализация m_inside
    void init_inside();

    /// @brief Схематический многоугольник, описывающий полосу
    std::vector<Vector3d> local_outline(int n_points) const final;

    double m_width; ///< Ширина полосы
};

/// @brief Круглое тело
class BodyDisk : public SolidBody2D {
public:
    /// @brief Круг в начале координат
    /// @param radius Радиус круга
    explicit BodyDisk(double radius);

    /// @brief Круг с указанием центра
    /// @param radius Радиус круга
    /// @param center Центр круга
    BodyDisk(double radius, const Vector3d& center);

    /// @brief Оптимизированная версия для круга
    double volume_fraction(const mesh::Cell& cell, double eps) const final;

    /// @brief Оптимизированная версия для круга
    double volume_inside(const mesh::Cell& cell, double eps) const final;

protected:
    /// @brief Инициализация m_inside
    void init_inside();

    /// @brief Многоугольник вписанный в окружность
    std::vector<Vector3d> local_outline(int n_points) const override;

    double m_radius;  ///< Радиус круга
};

/// @brief Квадратное тело
class BodySquare : public SolidBody2D {
public:
    /// @brief Квадрат в центре координат
    /// @param length Длина стороны квадрата
    explicit BodySquare(double length);

    /// @brief Квадрат с центром в определенной точке
    /// @param length Длина стороны квадрата
    /// @param center Центр квадрата
    BodySquare(double length, const Vector3d& center);

protected:
    /// @brief Инициализация m_inside
    void init_inside();

    /// @brief Массив из четырех вершин квадрата
    std::vector<Vector3d> local_outline(int n_points) const override;

    double m_length; ///< Длина стороны
};

} // namespace zephyr::geom