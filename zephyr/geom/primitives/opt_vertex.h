#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Трехмерный вектор-столбец из Eigen с дополнительными функциями.
/// Приставка Opt от слова опциональный. Поскольку данная вершина может
/// бать определена или нет.
class OptVertex : public Vector3d {
public:

    OptVertex() = default;

    template <class... Args>
    OptVertex(Args&&... args) : Vector3d(std::forward<Args>(args)...) { }

    //template<typename OtherDerived>
    //OptVertex(const Eigen::MatrixBase<OtherDerived> &other) {
    //    *reinterpret_cast<zephyr::geom::Vector3d *>(this) = other;
    //}

    template<typename OtherDerived>
    OptVertex &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        *reinterpret_cast<zephyr::geom::Vector3d *>(this) = other;
        return *this;
    }

    /// @brief Проверка первого элемента на NAN
    inline bool is_undefined() const {
        return x() != x();
    }

    /// @brief Проверка первого элемента на NAN
    inline bool is_actual() const {
        return x() == x();
    }

    /// @brief Сделать вершину неопределенной
    inline void set_undefined() {
        x() = y() = z() = NAN;
    }
};

} // namespace zephyr::geom