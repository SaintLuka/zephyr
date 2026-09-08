#pragma once

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

// forward declaration
class EuCell;
class EuMesh;

class EuFace_Iter;

/// @brief Итератор по инцидентным ячейкам
class EuIncCell_Iter final {
    AmrIncident* incident_;   ///< Список инцидентных ячеек
    index_t      inc_index_;  ///< Индекс в списке инцидентных

    AmrCells* local_cells_;   ///< Локальные ячейки
    AmrCells* ghost_cells_;   ///< Ячейки с других процессов

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    EuIncCell_Iter(AmrIncident* incident, index_t inc_index,
                   AmrCells* locals, AmrCells* ghosts = nullptr)
        : incident_(incident),
          inc_index_(inc_index),
          local_cells_(locals),
          ghost_cells_(ghosts) { }

    /// @brief Ссылка на грань при разыменовании
    EuCell operator*() const {
        if (incident_->ghost[inc_index_] < 0) {
            return EuCell(local_cells_, incident_->index[inc_index_]);
        }
        else {
            return EuCell(ghost_cells_, incident_->ghost[inc_index_]);
        }
    }

    /// @brief Перейти к следующей инцидентной ячейке
    EuIncCell_Iter &operator++() {
        ++inc_index_; return *this;
    }

    /// @brief Сравнение итераторов
    bool operator!=(const EuIncCell_Iter &cell_iter) const{
        return inc_index_ != cell_iter.inc_index_;
    }
};


/// @brief Интерфейс для итераций по инцидентным ячейкам
class EuIncidentCells final {
    AmrIncident* incident_;
    index_t inc_begin_;
    index_t inc_end_;
    AmrCells* locals_;
    AmrCells* ghosts_;

public:
    EuIncidentCells(AmrIncident *incident, index_t node_idx,
                    AmrCells *locals, AmrCells* ghosts);

    EuIncCell_Iter begin() const {
        return {incident_, inc_begin_, locals_, ghosts_};
    }

    EuIncCell_Iter end() const {
        return {incident_, inc_end_, locals_, ghosts_};
    }

    int size() const { return inc_end_ - inc_begin_; }
};

class EuNode_Iter;

/// @brief Узел сетки
/// @ingroup euler-mesh
class EuNode final {
    friend class EuNode_Iter;

    using Vector3d = geom::Vector3d;

private:
    AmrNodes* nodes_{nullptr};  //< Указатель на хранилище узлов
    index_t   index_{-1};       //< Индекс узла

    /// @brief Массивы нужны для прохода по инцидентным ячейкам
    AmrCells* locals_{nullptr};
    AmrCells* ghosts_{nullptr};

public:
    /// @brief Конструктор по умолчанию (иногда требует TBB)
    EuNode() = default;    

    EuNode(AmrNodes* nodes, index_t index, AmrCells* locals = nullptr, AmrCells* ghosts = nullptr)
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
    const AmrNodes& nodes() const { return *nodes_; }

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
    void copy_data_to(EuNode& dst_node) const;

    /// @}

    /// @{ @name Инцидентные ячейки

    /// @brief Число инцидентных ячеек
    int n_incident() const;

    /// @brief Итератор по инцидентным ячейкам
    EuIncidentCells incident() const;

    /// @}
};

/// @brief Итератор по узлам из EuMesh или AmrNodes
class EuNode_Iter final {
private:
    EuNode node_{};  ///< Реальная ячейка

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = EuNode;
    using pointer    = EuNode*;
    using reference  = EuNode&;

    /// @brief Конструктор по умолчанию (иногда требует TBB)
    EuNode_Iter() = default;

    /// @brief Конструктор как у узла
    EuNode_Iter(AmrNodes *nodes, index_t index, AmrCells* locals = nullptr, AmrCells *ghosts = nullptr)
            : node_{nodes, index, locals, ghosts} { }

    /// @brief Ссылка на узел при разыменовании
    EuNode &operator*() { return node_; }

    /// @brief Ссылка на узел при разыменовании
    const EuNode &operator*() const { return node_; }

    /// @brief Инкремент
    EuNode_Iter &operator++() {
        ++node_.index_;
        return *this;
    }

    /// @brief Декремент
    EuNode_Iter &operator--() {
        --node_.index_;
        return *this;
    }

    /// @brief Итератор через step
    EuNode_Iter &operator+=(index_t step) {
        node_.index_ += step;
        return *this;
    }

    /// @brief Итератор через step
    EuNode_Iter operator+(index_t step) const {
        return {node_.nodes_,
                node_.index_ + step,
                node_.locals_,
                node_.ghosts_};
    }

    /// @brief Оператор доступа как для указателя (random access iterator)
    EuNode operator[](index_t offset) const {
        return {node_.nodes_,
                node_.index_ + offset,
                node_.locals_,
                node_.ghosts_};
    }

    /// @brief Расстояние между двумя ячейками
    index_t operator-(const EuNode_Iter &cell) const {
        return node_.index_ - cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator<(const EuNode_Iter &cell) const {
        return node_.index_ < cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator!=(const EuNode_Iter &cell) const {
        return node_.index_ != cell.node_.index_;
    }

    /// @brief Оператор сравнения
    bool operator==(const EuNode_Iter &cell) const {
        return node_.index_ == cell.node_.index_;
    }
};

/// @brief Для итераций по AmrNodes, locals = ghosts = nullptr,
/// поэтому проход по инцидентным ячейкам невозможен
inline EuNode_Iter begin(AmrNodes& nodes) {
    return {&nodes, 0, nullptr, nullptr};
}

/// @brief Для итераций по AmrNodes, locals = ghosts = nullptr,
/// поэтому проход по инцидентным ячейкам невозможен
inline EuNode_Iter end(AmrNodes& nodes) {
    return {&nodes, nodes.n_nodes(), nullptr, nullptr};
}

// ================================================================================================
//                                     inline функции итераторов
// ================================================================================================

inline EuIncidentCells::EuIncidentCells(AmrIncident *incident,
    index_t node_idx, AmrCells *locals, AmrCells* ghosts)
    : incident_(incident), locals_(locals), ghosts_(ghosts) {

    z_assert(incident, "EuIncidentCells nullptr incident");
    z_assert(locals, "EuIncidentCells nullptr locals");
    z_assert(ghosts, "EuIncidentCells nullptr ghosts");
    z_assert(node_idx < incident->offsets.size() - 1, "EuIncidentCells node_idx out of range");

    inc_begin_ = incident->offsets[node_idx];
    inc_end_ = inc_begin_;

    index_t inc_max = incident->offsets[node_idx + 1];
    while (incident->role[inc_end_] >= 0 && inc_end_ < inc_max) {
        ++inc_end_;
    }
}

// ================================================================================================
//                                     inline функции EuNode
// ================================================================================================

inline int EuNode::rank() const { return nodes_->rank[index_]; }

inline index_t EuNode::id() const { return index_; }

inline index_t EuNode::index() const { return nodes_->index[index_]; }

inline index_t EuNode::next() const { return nodes_->next[index_]; }

inline void EuNode::set_rank(int rank) const { nodes_->rank[index_] = rank; }

template <typename T>
T& EuNode::operator[](Storable<T> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
std::span<T> EuNode::operator[](Storable<T[]> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<T, N> EuNode::operator[](Storable<T[N]> var) { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
const T& EuNode::operator[](Storable<T> var) const { return nodes_->data.get_val<T>(var, index_); }

template <typename T>
std::span<const T> EuNode::operator[](Storable<T[]> var) const { return nodes_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<const T, N> EuNode::operator[](Storable<T[N]> var) const { return nodes_->data.get_val<T>(var, index_); }

inline void EuNode::copy_data_to(EuNode &dst_node) const {
    nodes_->copy_data(index_, dst_node.nodes_, dst_node.index_);
}

inline int EuNode::n_incident() const { return nodes_->incident.max_count(index_); }

inline EuIncidentCells EuNode::incident() const {
    return EuIncidentCells(&nodes_->incident, index_, locals_, ghosts_);
}

} // namespace zephyr::mesh