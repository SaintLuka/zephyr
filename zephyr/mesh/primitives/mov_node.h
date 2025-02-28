#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/mesh/primitives/element.h>
#include <zephyr/mesh/primitives/adjacent.h>

namespace zephyr::mesh {

/// @class Узел подвижной сетки
class MovNode : public Element {
    using Vector3d = zephyr::geom::Vector3d;
    using Boundary = zephyr::geom::Boundary;
public:
    /// @brief Максимальное число смежных ячеек
    static const int max_neibs = 13;

    Boundary boundary;  ///< Граничное условие
    Vector3d coords;    ///< Положение узла
    Vector3d shift;     ///< Смещение узла при следующем вызове move()

    /// @brief Индексы смежных ячеек
    std::array<Adjacent, max_neibs> neibs;

    /// @brief Основной конструктор
    MovNode(int r = -1, int idx = 0) :
        Element(r, idx),
        boundary(Boundary::UNDEFINED) {
    };

    /// @brief Узел лежит на границе?
    inline bool is_boundary() const {
        return boundary != Boundary::ORDINARY &&
               boundary != Boundary::PERIODIC &&
               boundary != Boundary::UNDEFINED;
    }

    /// @brief Является ли узел актуальным?
    inline bool is_actual() const {
        return boundary != Boundary::UNDEFINED;
    }

    /// @return 'true', если узел не актуален не актуальна
    inline bool is_undefined() const {
        return boundary == Boundary::UNDEFINED;
    }

    /// @brief Установить неопределенную грань
    inline void set_undefined() {
        boundary = Boundary::UNDEFINED;
    }

    /// @brief Переместить узел
    void move() {
        coords += shift;
        shift = {0.0, 0.0, 0.0};
    }

    inline double x() const { return coords.x(); }

    inline double y() const { return coords.y(); }

    inline double z() const { return coords.z(); }
};

} // namespace zephyr::mesh