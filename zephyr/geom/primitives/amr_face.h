#pragma once

#include <string>

#include <zephyr/geom/primitives/base.h>
#include <zephyr/geom/primitives/base_face.h>
#include <zephyr/geom/primitives/amr_vertices.h>

namespace zephyr::geom {

struct AmrFace;

/// @brief Найти центр простой грани.
/// @param face Ссылка на грань.
/// @param verts Соответствующий список вершин.
template <int dim>
static Vector3d face_center(const AmrFace& face, const AmrVertices& verts);

/// @class Грань AMR ячейки, содержит не более 4 вершин
/// (в трехмерном случае).
class AmrFace : public BaseFace<4> {
public:

    /// @brief Конструктор по умолчанию
    AmrFace() = default;

    /// @brief Центр грани
    template <int dim>
    inline Vector3d center(const AmrVertices& verts) const {
        if (dim < 3) {
            return 0.5 * (verts[vertices[0]] + verts[vertices[1]]);
        } else {
            return 0.25 * (verts[vertices[0]] + verts[vertices[1]] +
                           verts[vertices[2]] + verts[vertices[3]]);
        }
    }

    /// @brief Центр грани
    inline Vector3d center(const AmrVertices& verts, int dim) const {
	    return dim < 3 ? center<2>(verts) : center<3>(verts);
    }
};

/// @brief Центр грани
template <int dim>
inline Vector3d face_center(const AmrFace& face, const AmrVertices& verts) {
    if (dim < 3) {
        return 0.5 * (verts[face.vertices[0]] + verts[face.vertices[1]]);
    } else {
        return 0.25 * (verts[face.vertices[0]] + verts[face.vertices[1]] +
                       verts[face.vertices[2]] + verts[face.vertices[3]]);
    }
}

} // namespace zephyr::geom