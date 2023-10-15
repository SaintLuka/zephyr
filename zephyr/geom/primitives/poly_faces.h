#pragma once

#include <zephyr/geom/primitives/poly_face.h>

namespace zephyr::geom {

/// @brief Список граней ячейки
class PolyFaces {
public:
    /// @brief Максимальное число граней
    static const int max_size = 24;

    /// @brief Актуальное число граней
    int n_faces;

    /// @brief Конструктор по умолчанию
    PolyFaces() = default;

private:
    /// @brief Массив граней ячейки
    std::array<PolyFace, max_size> m_list;
};

} // namespace zephyr::geom