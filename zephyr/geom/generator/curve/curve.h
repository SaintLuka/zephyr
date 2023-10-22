#pragma once

#include <memory>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/boundary.h>

namespace zephyr::geom::generator {

/// @brief Абстрактный класс, представляющий параметрическую кривую
/// X(t), Y(t) на плоскости. Параметр пробегает отрезок [-1, 1].
class Curve {
public:
    using Ptr = std::shared_ptr<Curve>;
    using Ref = const std::shared_ptr<Curve> &;

    /// @brief Конструктор
    Curve();

    /// @brief Координата X на параметрической кривой
    virtual double get_X(double t) const = 0;

    /// @brief Координата Y на параметрической кривой
    virtual double get_Y(double t) const = 0;

    /// @brief Точка на параметрической кривой
    virtual Vector3d get_V(double t) const;

    /// @brief Касательная к кривой по направлению возрастания параметра
    virtual Vector3d get_T(double t) const;

    /// @brief Нормаль к кривой, правый перпендикуляр к значению get_T(t).
    /// Если при обходе по кривой область остается слева, тогда нормаль
    /// является внешней нормалью области.
    virtual Vector3d get_N(double t) const;

    /// @brief Проекция точки v на кривую. Под проекцией понимается ближайшая
    /// к v точка кривой с параметром на отрезке [-1, 1].
    virtual Vector3d projection(const Vector3d &v) const = 0;

    /// @brief Установить флаг граничных условий
    void set_boundary(Boundary flag);

    /// @brief Получить флаг граничных условий
    Boundary boundary() const;

protected:
    Boundary m_boundary;
};

} // namespace zephyr::geom::generator
