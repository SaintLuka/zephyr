#pragma once

#include <zephyr/geom/vertices.h>
#include <zephyr/geom/faces.h>
#include <zephyr/geom/amr_data.h>

#include <zephyr/geom/geom.h>

namespace zephyr { namespace geom {

/// @struct Геометрические данные ячейки
class Cell {
public:
    short    dim;       ///< Размерность ячейки
    AmrData  amr;       ///< Данные адаптации
    Vector3d coords;    ///< Барицентр ячейки
    double   size;      ///< Линейный размер ячейки
    Vertices vertices;  ///< Вершины ячейки
    Faces    faces;     ///< Список граней

    /// @brief Двумерная простая
    explicit Cell(const ShortList2D &verts);

    /// @brief Двумерный полигон
    explicit Cell(const VerticesList &verts);

    /// @brief Двумерная криволинейная
    explicit Cell(const LargeList2D &verts);

    /// @brief Трехмерная простая
    explicit Cell(const ShortList3D &verts);

    /*
    /// @brief Строит родительскую ячейку по набору дочерних ячеек,
    /// дочерние ячейки должны быть упорядочены также, как и вершины
    /// в ячейке (указания в _ascii.h)
    /// @param children Массив итераторов дочерних ячеек
    explicit Cell(const std::array<Storage::iterator, 4> &children);

    explicit Cell(const std::array<Storage::iterator, 8> &children);
     */

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    double volume() const {
        return size * (dim < 3 ? size : size * size);
    }


private:

    /// Перенести вершины из списка в тип type::_vertices_
    void setup_vertices(const ShortList2D& vlist);

    /// Перенести вершины из списка в тип type::_vertices_
    void setup_vertices(const LargeList2D& vlist);

    /// Перенести вершины из списка в тип type::_vertices_
    void setup_vertices(const ShortList3D& vlist);

    void build2D(const ShortList2D &verts);

    void build2D(const VerticesList &verts);

    void build2D(const LargeList2D &verts);

    void build3D(const ShortList3D &verts);
};

} // namespace geom
} // namespace zephyr