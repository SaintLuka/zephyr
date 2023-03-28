#pragma once

#include <zephyr/mesh/generator/curve/curve.h>

namespace zephyr { namespace mesh { namespace generator {

class BaseVertex;

/// @brief Окружность, заданная параметрически
class Circle : public Curve {
private:
    using BaseVertex_Ref = const std::shared_ptr<BaseVertex> &;

public:
    /// @brief Окружность по радиусу и центру
    explicit Circle(double R, const Vector3d &v0 = {0.0, 0.0, 0.0});

    /// @brief Создать указатель на окружность по радиусу и центру
    static Curve::Ptr create(double R, const Vector3d &v0 = {0.0, 0.0, 0.0});

    /// @brief Создать указатель на окружность по трем базисным точкам окружности
    /// @throw Исключение, если точки принадлежат одной прямой
    static Curve::Ptr create(BaseVertex_Ref v1, BaseVertex_Ref v2, BaseVertex_Ref v3);

    /// @brief Координата X
    double get_X(double t) const final;

    /// @brief Координата Y
    double get_Y(double t) const final;

    /// @brief Точка окружности
    Vector3d get_V(double t) const final;

    /// @brief Касательная к окружности
    Vector3d get_T(double t) const final;

    /// @brief Нормаль к окружности
    Vector3d get_N(double t) const final;

    /// @brief Проекция точки на окружность
    Vector3d projection(const Vector3d &v) const final;

private:
    double R;
    Vector3d v0;
};

} // namespace generator
} // namespace mesh
} // namespace zephyr
