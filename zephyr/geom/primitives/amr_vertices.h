#pragma once

#include <zephyr/geom/maps.h>
#include <zephyr/geom/polygon.h>

namespace zephyr::geom {

class AmrVertices : public SqCube {
public:

    /*
    /// @brief Конструктор по умолчанию
    //AmrVertices() = default;

    /// @brief Конструктор по угловым точкам
    AmrVertices(const Vector3d &v000, const Vector3d &v002,
                const Vector3d &v020, const Vector3d &v022,
                const Vector3d &v200, const Vector3d &v202,
                const Vector3d &v220, const Vector3d &v222) 
                : SqCube(v000, v002, v020, v022, 
                         v200, v202, v220, v222) { };

    /// @brief Конструктор по полному набору узлов
    AmrVertices(const Vector3d &v000, const Vector3d &v001, const Vector3d &v002,
                const Vector3d &v010, const Vector3d &v011, const Vector3d &v012,
                const Vector3d &v020, const Vector3d &v021, const Vector3d &v022,
                const Vector3d &v100, const Vector3d &v101, const Vector3d &v102,
                const Vector3d &v110, const Vector3d &v111, const Vector3d &v112,
                const Vector3d &v120, const Vector3d &v121, const Vector3d &v122,
                const Vector3d &v200, const Vector3d &v201, const Vector3d &v202,
                const Vector3d &v210, const Vector3d &v211, const Vector3d &v212,
                const Vector3d &v220, const Vector3d &v221, const Vector3d &v222)
                : SqCube(v000, v001, v002, v010, v011, v012, v020, v021, v022,
                         v100, v101, v102, v110, v111, v112, v120, v121, v122,
                         v200, v201, v202, v210, v211, v212, v220, v221, v222) { };
                         */

    /// @brief Конструктор по угловым точкам
    AmrVertices(const Cube &cube) : SqCube(cube) { };

    /// @brief Конструктор по полному набору точек
    AmrVertices(const SqCube &cube) : SqCube(cube) { };

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    AmrVertices(const Quad& quad) : SqCube(quad) { };

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    AmrVertices(const SqQuad& quad) : SqCube(quad) { };

    /// @brief Конструктор по полигону
    AmrVertices(const Polygon &poly) {
        int max_size = verts.size();
        int mid_size = std::min(poly.size(), max_size);

        for (int i = 0; i < mid_size; ++i) {
            verts[i] = poly[i];
        }
        for (int i = mid_size; i < max_size; ++i) {
            verts[i] = {NAN, NAN, NAN};
        }
    };

    /// @brief Возвращает количество актуальных вершин.
    /// Функция используется для неадаптивных ячеек, чтобы узнать
    /// количество вершин примитива.
    inline int count() const {
        int max_size = verts.size();
        for (int i = 3; i < max_size; ++i) {
            if (verts[i].hasNaN()) {
                return i;
            }
        }
        return max_size;
    }
};

} // namespace zephyr::geom