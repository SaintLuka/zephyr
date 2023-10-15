#pragma once

#include <zephyr/geom/primitives/base_face.h>

namespace zephyr::geom {

/// @class Грань подвижной ячейки, содержит до 8 узлов
class PolyFace : public BaseFace<8> {
public:
    /// @brief Максимальное число узлов грани
    static const int max_nodes = 8;

    int n_nodes;  ///< Актуальное число узлов грани

    /// @brief Конструктор по умолчанию
    PolyFace() = default;
};

} // namespace zephyr::geom