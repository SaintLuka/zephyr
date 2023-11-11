#pragma once

#include <zephyr/geom/primitives/element.h>
#include <zephyr/geom/primitives/amr_vertices.h>
#include <zephyr/geom/primitives/amr_faces.h>

namespace zephyr::geom {

/// @class Обязательные данные ячейки сетки
class AmrCell : public Element {
public:
    // Геометрия ячейки

    int         dim;       ///< Размерность ячейки
    bool        adaptive;  ///< Адаптивная ячейка?
    bool        linear;    ///< Линейная ячейка?
    Vector3d    center;    ///< Барицентр ячейки
    double      size;      ///< Линейный размер ячейки
    AmrVertices vertices;  ///< Вершины ячейки
    AmrFaces    faces;     ///< Список граней

    // Данные AMR

    int b_idx;  ///< Индекс среди базовых ячеек
    int z_idx;  ///< Индекс ячейки на z-кривой
    int level;  ///< Уровень адаптации (0 для базовой)
    int flag;   ///< Желаемый флаг адаптации


    /// @brief Двумерная простая
    explicit AmrCell(const Quad &quad);

    /// @brief Двумерная криволинейная
    explicit AmrCell(const SqQuad &quad);

    /// @brief Трехмерная простая
    explicit AmrCell(const Cube &cube);

    /// @brief Трехмерная криволинейная ячейка
    explicit AmrCell(const SqCube &cube);

    /// @brief Двумерная полигональная ячейка. Не адаптивная ячейка, может
    /// представлять четырехугольник, но вершины упорядочены иначе.
    explicit AmrCell(const Polygon& poly);

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    double volume() const;

    /// @brief Скорпировать вершины в полигон (двумерные ячейки)
    /// Для нелинейных AMR-ячеек возвращает до 8 граней.
    PolygonS<8> polygon() const;

    /// @brief Вывести информацию о ячейке
    void print_info() const;

    /// @brief Вывести информацию о ячейке в виде python скрипта
    /// для визуализации
    void visualize(std::string filename) const;

    /// @brief Проверить базовую геометрию ячейки
    /// @return -1 для плохой ячейки
    int check_geometry() const;

    /// @brief Проверить ориентацию граней
    /// @return -1 для плохой ячейки
    int check_base_face_orientation() const;

    /// @brief Проверить порядок вершин
    /// @return -1 для плохой ячейки
    int check_base_vertices_order() const;

    /// @brief Проверить сложные грани
    /// @return -1 для плохой ячейки
    int check_complex_faces() const;

};

} // namespace zephyr::geom