#pragma once

#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

class BaseVertex;

/// @brief Прямая линия с параметрическим заданием X(t), Y(t).
/// Полная прямая получается при пробегании t от -1 до 1 (используется
/// отображение в виде тангенса).
class Line : public Curve {
private:
    using BaseVertex_Ref = const std::shared_ptr<BaseVertex> &;

public:
    /// @brief Конструктор по двум точкам
    Line(const Vector3d &v1, const Vector3d &v2);

    /// @brief Созадние указателя по двум точкам
    static Curve::Ptr create(const Vector3d &v1, const Vector3d &v2);

    /// @brief Создание указателя по двум базисным точкам
    static Curve::Ptr create(BaseVertex_Ref v1, BaseVertex_Ref v2);

    /// @brief Координата X
    double get_X(double t) const final;

    /// @brief Координата Y
    double get_Y(double t) const final;

    /// @brief Точка на прямой
    Vector3d get_V(double t) const final;

    /// @brief Касательная
    Vector3d get_T(double t) const final;

    /// @brief Нормаль
    Vector3d get_N(double t) const final;

    /// @brief Проекция точки на прямую
    Vector3d projection(const Vector3d &v) const final;

private:
    Vector3d v1;
    Vector3d v2;
};

} // namespace zephyr::geom::generator
