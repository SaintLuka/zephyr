#pragma once

#include <zephyr/mesh/cell.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

// forward declaration
class Cell;
class Mesh;

class Face_Iter;

/// @brief Итератор по инцидентным ячейкам
class IncCell_Iter final {
    RawIncident* incident_;   ///< Список инцидентных ячеек
    index_t      inc_index_;  ///< Индекс в списке инцидентных

    RawCells* local_cells_;   ///< Локальные ячейки
    RawCells* ghost_cells_;   ///< Ячейки с других процессов

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    IncCell_Iter(RawIncident* incident, index_t inc_index,
                   RawCells* locals, RawCells* ghosts = nullptr)
        : incident_(incident),
          inc_index_(inc_index),
          local_cells_(locals),
          ghost_cells_(ghosts) { }

    /// @brief Ссылка на грань при разыменовании
    Cell operator*() const {
        if (incident_->ghost[inc_index_] < 0) {
            return Cell(local_cells_, incident_->index[inc_index_]);
        }
        else {
            return Cell(ghost_cells_, incident_->ghost[inc_index_]);
        }
    }

    /// @brief Перейти к следующей инцидентной ячейке
    IncCell_Iter &operator++() {
        ++inc_index_; return *this;
    }

    /// @brief Сравнение итераторов
    bool operator!=(const IncCell_Iter &cell_iter) const{
        return inc_index_ != cell_iter.inc_index_;
    }
};


/// @brief Интерфейс для итераций по инцидентным ячейкам
class IncidentCells final {
    RawIncident* incident_;
    index_t inc_begin_;
    index_t inc_end_;
    RawCells* locals_;
    RawCells* ghosts_;

public:
    IncidentCells(RawIncident *incident, index_t node_idx,
                    RawCells *locals, RawCells* ghosts);

    IncCell_Iter begin() const {
        return {incident_, inc_begin_, locals_, ghosts_};
    }

    IncCell_Iter end() const {
        return {incident_, inc_end_, locals_, ghosts_};
    }

    int size() const { return inc_end_ - inc_begin_; }
};

/// @brief Узел сетки
/// @ingroup raw-mesh
class Node final {
    friend class Node_Iter;

    using Vector3d = geom::Vector3d;

private:
    RawNodes* nodes_{nullptr};  //< Указатель на хранилище узлов
    index_t   index_{-1};       //< Индекс узла

    /// @brief Массивы нужны для прохода по инцидентным ячейкам
    RawCells* locals_{nullptr};
    RawCells* ghosts_{nullptr};

public:
    /// @brief Конструктор по умолчанию (иногда требует TBB)
    Node() = default;

    Node(RawNodes* nodes, index_t index, RawCells* locals = nullptr, RawCells* ghosts = nullptr)
        : nodes_(nodes), index_(index), locals_(locals), ghosts_(ghosts) { }

    /// @{ @name Характеристики узла

    /// @brief Ранг, которому принадлежит ячейка
    int rank() const;

    /// @brief Реальное положение в массиве
    index_t id() const;

    /// @brief Текущий индекс ячейки в массиве
    index_t index() const;

    /// @brief Индекс новой ячейки (в алгоритмах)
    index_t next() const;

    /// @brief Установить ранг ячейки
    void set_rank(int rank) const;

    /// @brief Хранилище, которому принадлежит узел
    const RawNodes& nodes() const { return *nodes_; }

    /// @}

    /// @{ @name Геометрия ячейки
        
    /// @brief Координата узла (неявное приведение)
    operator Vector3d& () { return nodes_->coord[index_]; }

    /// @brief Координата узла (неявное приведение)
    operator const Vector3d& () const { return nodes_->coord[index_]; }
    
    /// @brief Координата узла
    const Vector3d& coord() const { return nodes_->coord[index_]; }
    
    double x() const { return nodes_->coord[index_].x(); }
    double y() const { return nodes_->coord[index_].y(); }
    double z() const { return nodes_->coord[index_].z(); }

    /// @}

    /// @{ @name Данные узла

    /// @brief Ссылка на данные узла
    template <typename T>
    T& operator[](Storable<T> var);

    /// @brief Ссылка на данные узла
    template <typename T>
    std::span<T> operator[](Storable<T[]> var);

    /// @brief Ссылка на данные узла
    template <typename T, size_t N>
    std::span<T, N> operator[](Storable<T[N]> var);

    /// @brief Константная ссылка на данные узла
    template <typename T>
    const T& operator[](Storable<T> var) const;

    /// @brief Константная ссылка на данные узла
    template <typename T>
    std::span<const T> operator[](Storable<T[]> var) const;

    /// @brief Константная ссылка на данные узла
    template <typename T, size_t N>
    std::span<const T, N> operator[](Storable<T[N]> var) const;

    /// @brief Скопировать данные в другой узел
    void copy_data_to(Node& dst_node) const;

    /// @}

    /// @{ @name Инцидентные ячейки

    /// @brief Число инцидентных ячеек
    int n_incident() const;

    /// @brief Итератор по инцидентным ячейкам
    IncidentCells incident() const;

    /// @}
};

/// @brief Итератор по узлам из Mesh или RawNodes
class Node_Iter final {
private:
    Node node_{};  ///< Реальная ячейка

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = Node;
    using pointer    = Node*;
    using reference  = Node&;

    /// @brief Конструктор по умолчанию (иногда требует TBB)
    Node_Iter() = default;

    /// @brief Конструктор как у узла
    Node_Iter(RawNodes *nodes, index_t index, RawCells* locals = nullptr, RawCells *ghosts = nullptr)
            : node_{nodes, index, locals, ghosts} { }

    /// @brief Ссылка на узел при разыменовании
    Node &operator*() { return node_; }

    /// @brief Ссылка на узел при разыменовании
    const Node &operator*() const { return node_; }

    /// @brief Инкремент
    Node_Iter &operator++() {
        ++node_.index_;
        return *this;
    }

    /// @brief Декремент
    Node_Iter &operator--() {
        --node_.index_;
        return *this;
    }

    /// @brief Итератор через step
    Node_Iter &operator+=(index_t step) {
        node_.index_ += step;
        return *this;
    }

    /// @brief Итератор через step
    Node_Iter operator+(index_t step) const {
        return {node_.nodes_,
                node_.index_ + step,
                node_.locals_,
                node_.ghosts_};
    }

    /// @brief Оператор доступа как для указателя (random access iterator)
    Node operator[](index_t offset) const {
        return {node_.nodes_,
                node_.index_ + offset,
                node_.locals_,
                node_.ghosts_};
    }

    /// @brief Расстояние между двумя ячейками
    index_t operator-(const Node_Iter &cell) const {
        return node_.index_ - cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator<(const Node_Iter &cell) const {
        return node_.index_ < cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator!=(const Node_Iter &cell) const {
        return node_.index_ != cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator==(const Node_Iter &cell) const {
        return node_.index_ == cell.node_.index_;
    }
};

/// @brief Для итераций по RawNodes, locals = ghosts = nullptr,
/// поэтому проход по инцидентным ячейкам невозможен
inline Node_Iter begin(RawNodes& nodes) {
    return {&nodes, 0, nullptr, nullptr};
}

/// @brief Для итераций по RawNodes, locals = ghosts = nullptr,
/// поэтому проход по инцидентным ячейкам невозможен
inline Node_Iter end(RawNodes& nodes) {
    return {&nodes, nodes.n_nodes(), nullptr, nullptr};
}

class NodesRange {
private:
    RawNodes* nodes_;
    RawCells* local_cells_{nullptr};
    RawCells* ghost_cells_{nullptr};

public:
    NodesRange(RawNodes* nodes, RawCells* local_cells, RawCells* ghost_cells)
        : nodes_{nodes}, local_cells_{local_cells}, ghost_cells_{ghost_cells} {
    }

    Node_Iter begin() const {
        return Node_Iter(nodes_, 0, local_cells_, ghost_cells_);
    }

    Node_Iter end() const {
        return Node_Iter(nodes_, nodes_->n_nodes(), local_cells_, ghost_cells_);
    }
};

// ================================================================================================
//                                     inline функции итераторов
// ================================================================================================

inline IncidentCells::IncidentCells(RawIncident *incident,
    index_t node_idx, RawCells *locals, RawCells* ghosts)
    : incident_(incident), locals_(locals), ghosts_(ghosts) {

    z_assert(incident, "IncidentCells nullptr incident");
    z_assert(locals, "IncidentCells nullptr locals");
    z_assert(ghosts, "IncidentCells nullptr ghosts");
    z_assert(node_idx < incident->offsets.size() - 1, "IncidentCells node_idx out of range");

    inc_begin_ = incident->offsets[node_idx];
    inc_end_ = inc_begin_;

    index_t inc_max = incident->offsets[node_idx + 1];
    while (incident->role[inc_end_] >= 0 && inc_end_ < inc_max) {
        ++inc_end_;
    }
}

// ================================================================================================
//                                     inline функции Node
// ================================================================================================

inline int Node::rank() const { return nodes_->rank[index_]; }

inline index_t Node::id() const { return index_; }

inline index_t Node::index() const { return nodes_->index[index_]; }

inline index_t Node::next() const { return nodes_->next[index_]; }

inline void Node::set_rank(int rank) const { nodes_->rank[index_] = rank; }

template <typename T>
T& Node::operator[](Storable<T> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
std::span<T> Node::operator[](Storable<T[]> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<T, N> Node::operator[](Storable<T[N]> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
const T& Node::operator[](Storable<T> var) const { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
std::span<const T> Node::operator[](Storable<T[]> var) const { return nodes_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<const T, N> Node::operator[](Storable<T[N]> var) const { return nodes_->data.get_val<T>(var, index_); }

inline void Node::copy_data_to(Node &dst_node) const {
    nodes_->copy_data(index_, dst_node.nodes_, dst_node.index_);
}

inline int Node::n_incident() const { return nodes_->incident.max_count(index_); }

inline IncidentCells Node::incident() const {
    return IncidentCells(&nodes_->incident, index_, locals_, ghosts_);
}

} // namespace zephyr::mesh