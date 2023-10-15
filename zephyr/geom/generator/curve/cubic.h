#pragma once

#include <vector>

#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

class BaseVertex;

/// @brief Кривая в виде кубического сплайна от параметра.
/// То есть функции X(t), Y(t) являются кубическими сплайнами от параметра
/// t на отрезке [-1, 1]. Сплайн задается путем определения четырех точек
/// в узлах t = -1, -1/3, 1/3, 1.
class Cubic : public Curve {
private:
    using BaseVertex_Ref = const std::shared_ptr<BaseVertex> &;

public:
    /// @brief Конструктор сплайна по четырем точкам
    Cubic(const Vector3d &v1, const Vector3d &v2,
          const Vector3d &v3, const Vector3d &v4);

    /// @brief Создать указатель на сплайн по четырем точкам
    static Curve::Ptr create(const Vector3d &v1, const Vector3d &v2,
                             const Vector3d &v3, const Vector3d &v4);

    /// @brief Создать указатель на сплайн по четырем базисным точкам
    static Curve::Ptr create(BaseVertex_Ref v1, BaseVertex_Ref v2,
                             BaseVertex_Ref v3, BaseVertex_Ref v4);

    /// @brief Координата X
    double get_X(double t) const final;

    /// @brief Координата Y
    double get_Y(double t) const final;

    /// @brief Точка на кривой
    Vector3d get_V(double t) const final;

    /// @brief Касательная
    Vector3d get_T(double t) const final;

    /// @brief Нормаль
    Vector3d get_N(double t) const final;

    /// @brief Проекция на кривую
    Vector3d projection(const Vector3d &v) const final;

private:
    /// @brief Первая производная
    double dX(double t) const;

    /// @brief Первая производная
    double dY(double t) const;

    /// @brief Вторая производная
    double d2X(double t) const;

    /// @brief Вторая производная
    double d2Y(double t) const;

    /// @brief Производная расстояния от точки v до кривой
    double obj_func(const Vector3d &v, double t) const;

    /// @brief Вторая производная расстояния от точки v до кривой
    double obj_func_deriv(const Vector3d &v, double t) const;

    Vector3d v1;
    Vector3d v2;
    Vector3d v3;
    Vector3d v4;

    /// @brief Набор точек на сплайне
    std::vector<Vector3d> m_points;
};

} // namespace zephyr::geom::generator
