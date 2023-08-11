#pragma once

#include <zephyr/geom/faces.h>
#include <zephyr/geom/geom.h>

namespace zephyr { namespace geom {

/// @struct Обязательные данные ячейки сетки
class Cell {
public:
    // Данные о ячейке на сетке

    int rank;   ///< Ранг процесса
    int index;  ///< Индекс ячейки в Storage (>= 0)

    // Данные AMR

    int b_idx;  ///< Индекс среди базовых ячеек
    int z_idx;  ///< Индекс ячейки на z-кривой
    int next;   ///< Индекс новой ячейки (в алгоритмах)
    int level;  ///< Уровень адаптации (0 для базовой)
    int flag;   ///< Желаемый флаг адаптации

    // Геометрия ячейки

    int      dim;       ///< Размерность ячейки
    Vector3d coords;    ///< Барицентр ячейки
    double   size;      ///< Линейный размер ячейки
    Vertices vertices;  ///< Вершины ячейки
    Faces    faces;     ///< Список граней

    /// @brief Конструктор по умолчанию, инициализирует
    /// нулями не геометрические данные ячейки
    Cell();

    /// @brief Двумерная простая
    explicit Cell(const ShortList2D &verts);

    /// @brief Двумерный полигон
    explicit Cell(const VerticesList &verts);

    /// @brief Двумерная криволинейная
    explicit Cell(const LargeList2D &verts);

    /// @brief Трехмерная простая
    explicit Cell(const ShortList3D &verts);

    /// @brief Скопировать геометрию ячейки в другую
    void copy_to(Cell& cell) const;

    /*
    /// @brief Строит родительскую ячейку по набору дочерних ячеек,
    /// дочерние ячейки должны быть упорядочены также, как и вершины
    /// в ячейке (указания в _ascii.h)
    /// @param children Массив итераторов дочерних ячеек
    explicit Cell(const std::array<Storage::iterator, 4> &children);

    explicit Cell(const std::array<Storage::iterator, 8> &children);
     */

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    double volume() const;

    /// @brief Устанавливает rank = -1 (ячейка вне сетки)
    void set_undefined();

    /// @brief Актуальная ячейка? (rank >= 0)
    bool is_actual() const;

    /// @brief Ячейка к удалению (rank < 0)
    bool is_undefined() const;

    /// @brief Вывести информацию о ячейке
    void print_info() const;

    /// @brief Вывести информацию о ячейке в виде python скрипта
    /// для визуализации
    void visualize() const;

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

private:

    /// Перенести вершины из списка в тип Vertices
    void setup_vertices(const ShortList2D& vlist);

    /// Перенести вершины из списка в тип Vertices
    void setup_vertices(const LargeList2D& vlist);

    /// Перенести вершины из списка в тип Vertices
    void setup_vertices(const ShortList3D& vlist);

    void build2D(const ShortList2D &verts);

    void build2D(const VerticesList &verts);

    void build2D(const LargeList2D &verts);

    void build3D(const ShortList3D &verts);
};

} // namespace geom
} // namespace zephyr