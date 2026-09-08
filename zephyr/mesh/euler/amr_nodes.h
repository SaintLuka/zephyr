#pragma once

#include <vector>

#include <zephyr/utils/range.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/index.h>

namespace zephyr::mesh {

class AmrCells;
class AmrVerts;

using role_t = std::int8_t;

/// @brief Индексы инцидентных ячеек
/// @details Короткое объяснение. Если инцидентная ячейка
///             с этого процесса  |  с другого процесса
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < ghosts.size
class AmrIncident final {
public:
    static constexpr int max_amr_count_2D = 5;
    static constexpr int max_amr_count_3D = 10;

    static constexpr int max_amr_count(int dim) {
        return dim < 3 ? max_amr_count_2D : max_amr_count_3D;
    }


    /// @brief Эта структура имеет формат CSR, даны смещения
    std::vector<index_t> offsets = {0};

    /// @brief Роль узла внутри ячейки (индекс внутри ячейки)
    std::vector<role_t> role;

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    std::vector<int> rank;

    /// @brief Индекс ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    std::vector<index_t> index;

    /// @brief Индекс смежной ячейки в массиве ghosts (или -1, если соседняя
    /// ячейка с данного процесса).
    std::vector<index_t> ghost;


    /// @brief Пустые массивы по умолчанию
    AmrIncident() = default;

    /// @brief Количество значений/записей
    index_t n_values() const { return rank.size(); }

    /// @brief Очистить массивы
    void clear();

    /// @brief Расширить массивы по числу граней
    void resize(index_t n_nodes, index_t n_values);

    /// @brief Расширить массивы по числу граней
    void reserve(index_t n_nodes, index_t n_values);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Увеличить размер под массив для AMR узлов
    void resize_amr(index_t n_nodes, int dim);

    /// @brief Число инцидентных ячеек, для адаптивной ячейки может быть меньше
    /// max_count, для неструктурированной ячейки (полигон или многогранник
    /// совпадает с max_count).
    int count(index_t inode) const;

    /// @brief Максимальное число граней для ячейки
    int max_count(index_t inode) const {
        return offsets[inode + 1] - offsets[inode];
    }

    /// @brief Полный диапазон соседних ячеек (могут встречаться неактуальные)
    range_t<index_t> range(index_t inode) const {
        return std::views::iota(offsets[inode], offsets[inode + 1]);
    }

    /// @brief Является ли грань актуальной?
    bool is_actual(index_t inc) const { return role[inc] >= 0; }

    /// @return 'true', если грань не актуальна
    bool is_undefined(index_t inc) const { return role[inc] < 0; }

    /// @brief Установить неопределенную грань
    void set_undefined(index_t inc) { role[inc] = -1; }

    /// @brief Локальная соседняя ячейка?
    bool is_local(index_t inc) const { return ghost[inc] < 0; }

    /// @brief Удаленная соседняя ячейка?
    bool is_ghost(index_t inc) const { return ghost[inc] >= 0; }

    /// @brief Получить хранилище ячеек, в котором находится сосед, а также
    /// индекс соседа в данном хранилище
    template <class SomeArray>
    std::tuple<const SomeArray &, index_t> get_neib(index_t inc,
            const SomeArray &locals, const SomeArray &ghosts) const {
        if (ghost[inc] < 0) {
            return {locals, index[inc]};
        } else {
            return {ghosts, ghost[inc]};
        }
    }

    /// @brief Расход памяти
    memory_t memory_usage() const;
};

/// @brief Список уникальных узлов, дополняет класс AmrCells
class AmrNodes final {
    // aliases inside class
    using Vector3d = geom::Vector3d;

public:
    std::vector<int>     rank;    ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> next;    ///< Новый индекс в хранилище (в алгоритмах с перестановками)
    std::vector<index_t> index;   ///< Глобальный индекс узла в локальном хранилище
                                  /// (< 0 для неопределенных узлов, узлов на удаление)

    /// @brief Координаты узлов
    std::vector<Vector3d> coord;

    /// @brief Списки инцидентных ячеек
    AmrIncident incident;

    /// @brief Данные узлов
    Storage data;


    /// @brief Конструктор по умолчанию
    AmrNodes() = default;

    /// @brief Пустое хранилище узлов?
    bool empty() const { return coord.empty(); }

    /// @brief Число уникальных узлов
    index_t size() const { return coord.size(); }

    /// @brief Число уникальных узлов
    index_t n_nodes() const { return coord.size(); }

    /// @brief Размер списка инцидентных
    index_t n_incident() const { return incident.n_values(); }

    /// @brief Очистить массивы
    void clear();

    /// @brief Расширить массивы по числу узлов
    void resize(index_t n_nodes, index_t n_incident);

    /// @brief Расширить массивы по числу узлов
    void reserve(index_t n_nodes, index_t n_incident);

    /// @brief Увеличить размер под массив для AMR узлов
    void resize_amr(index_t n_nodes, int dim);

    void shrink_to_fit();

    // {vert.index, nodes}
    using Incomplete = std::tuple<AmrVerts, AmrNodes>;

    template <bool complete>
    static Incomplete generate(const AmrCells& cells);

    void setup_for(AmrCells& cells);

    /// @{ @name Топологические свойства узлов

    /// @brief Актуальный узел?
    bool is_actual(index_t in) const { return index[in] >= 0; }

    /// @brief Узел к удалению
    bool is_undefined(index_t in) const { return index[in] < 0; }

    /// @brief Устанавливает index = -1 (узел вне сетки)
    void set_undefined(index_t in) { index[in] = -1; }

    /// @}

    /// @brief Скопировать все данные с индекса from на индекс to.
    void copy_data(index_t from, index_t to);

    /// @brief Скопировать все данные целиком с индекса from,
    /// в хранилище dst на индекс to
    void copy_data(index_t from, AmrNodes* dst, index_t to) const;

    void copy_geom(index_t ic, AmrNodes& nodes,
        index_t jc, index_t inc_offset) const;

    /// @brief Расход памяти
    memory_t memory_usage() const;

    /// @brief Проверить согласованность размеров
    int check_sizes() const;

    /// @brief Проверка уникальных узлов для однопроцессорной версии
    int check_nodes(const AmrCells& locals) const;

    /// @brief Проверка уникальных узлов в MPI версии
    int check_nodes(const AmrCells& locals, const AmrCells& ghosts, const AmrNodes& ghost_nodes) const;
};

extern template AmrNodes::Incomplete AmrNodes::generate<true> (const AmrCells& cells);
extern template AmrNodes::Incomplete AmrNodes::generate<false>(const AmrCells& cells);

} // namespace zephyr::mesh