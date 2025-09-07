#pragma once

#include <memory>
#include <vector>
#include <set>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/cell_type.h>

namespace zephyr::mesh {
class AmrCell;
class MovNode;
class MovCell;
}

namespace zephyr::geom {


/// @brief Узел сетки общего вида
class GNode {
public:
    /// @brief Умный указатель на узел
    using Ptr = std::shared_ptr<GNode>;
    using Ref = const std::shared_ptr<GNode> &;

    int index = -1; ///< Индекс в массиве
    Vector3d v;     ///< Положение вершины

    /// @brief Конструктор
    explicit GNode(const Vector3d &v);

    /// @brief Конструктор
    explicit GNode(double x, double y, double z = 0.0);

    /// @brief Создание умного указателя
    static GNode::Ptr create(const Vector3d &v);

    /// @brief Создание умного указателя
    static GNode::Ptr create(double x, double y, double z = 0.0);

    /// @brief Установить граничное условие
    void add_boundary(Boundary flag);

    /// @brief Установить граничное условие
    void set_boundaries(std::set<Boundary> flags);

    /// @brief Граничное условие на грани, образованной несколькими вершинами
    /// @details Ищет пересечение флагов граничных условий
    static Boundary face_boundary(const std::vector<GNode::Ptr> &vs);

    /// @brief Список индексов ячеек, которые являются общими для набора
    /// вершин. Если vs образуют грань, тогда возвращаемый вектор имеет от
    /// размер от 0 до 2.
    /// @details Фактически ищет пересечения в множествах смежных ячеек
    static std::vector<int> shared_cells(const std::vector<GNode::Ptr> &vs);

    /// @brief Добавить соседа
    void add_neib_cell(int idx);

private:
    /// @brief Множество граничных условий, для узла может быть более одного
    /// граничного условия
    std::set<Boundary> m_bounds;

    /// @brief Список ячеек, которые содержат данную вершину
    std::set<int> m_neibs;
};

/// @brief Грань - индексы вершин в массиве в ячейке
using GFace = std::vector<int>;

class GCell {
public:

    int index = -1;

    static GCell triangle(const std::array<GNode::Ptr, 3>& nodes);

    static GCell quad(const std::array<GNode::Ptr, 4>& nodes);

    static GCell polygon(std::vector<GNode::Ptr> nodes);

    static GCell tetra(const std::array<GNode::Ptr, 4>& nodes);

    static GCell pyramid(const std::array<GNode::Ptr, 5>& nodes);

    static GCell wedge(const std::array<GNode::Ptr, 6>& nodes);

    static GCell hexagedron(const std::array<GNode::Ptr, 8>& nodes);


    CellType type() const;

    GNode& node(int idx);

    const GNode& node(int idx) const;

    int adjacent(int side) const;

    int adjacent(const std::vector<GNode::Ptr>& face) const;

    void add_self_to_nodes();

    void find_neibs();

    int n_nodes() const;

    int n_faces() const;

    std::vector<GNode::Ptr> face_nodes(int idx) const;

    Boundary boundary(const std::vector<GNode::Ptr>& face) const;


    int& neib(int idx);

    int neib(int idx) const;



private:
    GCell();

    GCell(CellType type);

    CellType m_type;
    std::vector<GNode::Ptr> m_nodes;
    std::vector<GFace>      m_faces;
    std::vector<int>        m_neibs;
};


/// @brief Сетка общего вида, которую выдают сеточные генераторы.
/// Непосредственно в расчетах не используется, после создания преобразуется
/// в специализированный сеточный класс.
class Grid {
public:
    Grid();

    void set_axial(bool axial = true);

    int n_nodes() const;

    int n_cells() const;

    void reserve_nodes(int size);

    void reserve_cells(int size);

    void operator +=(GNode::Ref node);

    void operator +=(const GCell& cell);


    GNode::Ptr node(int idx);

    GNode::Ref node(int idx) const;

    GCell& cell(int idx);

    const GCell& cell(int idx) const;

    void setup_adjacency();

    // Интересный код
    //mesh::AmrCell amr_cell(int idx) const;

    void assume_structured(int nx, int ny, int nz = 1);

    bool is_structured() const;

    int nx() const;

    int ny() const;

    int nz() const;

private:
    std::vector<GNode::Ptr> m_nodes;
    std::vector<GCell> m_cells;

    ///
    bool m_axial = false;

    /// Поля для структурированной сетки
    bool structured = false;
    int m_nx, m_ny, m_nz;
};

} // namespace zephyr::geom

