#pragma once

#include <zephyr/geom/maps.h>
#include <zephyr/geom/polygon.h>

namespace zephyr::geom {

class BVertices : public SqCube {
public:
    ///@brief Максимальное число вешин
    static const int max_count = 27;

    /*
    /// @brief Конструктор по умолчанию
    //BVertices() = default;

    /// @brief Конструктор по угловым точкам
    BVertices(const Vector3d &v000, const Vector3d &v002,
                const Vector3d &v020, const Vector3d &v022,
                const Vector3d &v200, const Vector3d &v202,
                const Vector3d &v220, const Vector3d &v222) 
                : SqCube(v000, v002, v020, v022, 
                         v200, v202, v220, v222) { };

    /// @brief Конструктор по полному набору узлов
    BVertices(const Vector3d &v000, const Vector3d &v001, const Vector3d &v002,
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
    BVertices(const Cube &cube) : SqCube(cube) { };

    /// @brief Конструктор по полному набору точек
    BVertices(const SqCube &cube) : SqCube(cube) { };

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    BVertices(const Quad& quad) : SqCube(quad) { };

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    BVertices(const SqQuad& quad) : SqCube(quad) { };

    /// @brief Конструктор по полигону
    BVertices(const Polygon &poly);

    /// @brief Найти индекс вершины в списке
    /// @param v Целевая вершина
    /// @param eps Точность (близость вершин)
    /// @return -1, если вершина не найдена
    int find(const Vector3d& v, double eps) const;

    /// @brief Возвращает количество актуальных вершин.
    /// Функция используется для неадаптивных ячеек, чтобы узнать
    /// количество вершин примитива.
    inline int count() const {
        for (int i = 3; i < max_count; ++i) {
            if (verts[i].hasNaN()) {
                return i;
            }
        }
        return max_count;
    }
};

static_assert(sizeof(BVertices) == sizeof(SqCube));
static_assert(BVertices::max_count == sizeof(BVertices) / (sizeof(Vector3d)));

} // namespace zephyr::geom