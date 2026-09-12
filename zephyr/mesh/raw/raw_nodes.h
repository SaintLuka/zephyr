#pragma once

#include <vector>

#include <zephyr/utils/range.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/index.h>

namespace zephyr::mesh {

class RawCells;
class RawVerts;

using role_t = std::int8_t;

/// @brief Индексы инцидентных ячеек
/// @details Короткое объяснение. Если инцидентная ячейка
///             с этого процесса  |  с другого процесса
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < ghosts.size
class RawIncident final {
public:
    /// @brief Эта структура имеет формат CSR, даны смещения
    std::vector<index_t> offsets = {0};

    /// @brief Роль узла внутри ячейки (индекс внутри ячейки).
    /// Значение role = -1, для неактуальных записей.
    std::vector<role_t> role;

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    std::vector<int> rank;

    /// @brief Индекс ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    std::vector<index_t> index;

    /// @brief Индекс смежной ячейки в массиве ghosts (или -1, если
    /// соседняя ячейка с данного процесса).
    std::vector<index_t> ghost;


    /// @brief Пустые массивы по умолчанию
    RawIncident() = default;

    /// @brief Количество значений/записей
    index_t n_values() const { return static_cast<index_t>(rank.size()); }

    /// @brief Очистить массивы
    void clear();

    /// @brief Расширить массивы
    void resize(index_t n_nodes, index_t n_values);

    /// @brief Увеличить размер под массив для AMR узлов
    void resize_amr(index_t n_nodes, int dim);

    /// @brief Расширить буферы массивов
    void reserve(index_t n_nodes, index_t n_values);

    /// @brief Увеличить размер под массив для AMR узлов
    void reserve_amr(index_t n_nodes, int dim);

    /// @brief Сжать буферы массивов до актуальных размеров
    void shrink_to_fit();

    /// @brief Максимальное число инцидентных ячеек для AMR-сетки
    static constexpr int max_incident_amr(int dim) {
        return dim < 3 ? 5 : 10;
    }

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

    /// @brief Является ли запись актуальной?
    bool is_actual(index_t inc) const { return role[inc] >= 0; }

    /// @return 'true', если запись не актуальна
    bool is_undefined(index_t inc) const { return role[inc] < 0; }

    /// @brief Установить неопределенную инцидентную
    void set_undefined(index_t inc) { role[inc] = -1; }

    /// @brief Локальная инцидентная ячейка?
    bool is_local(index_t inc) const { return ghost[inc] < 0; }

    /// @brief Удаленная инцидентная ячейка?
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

/// @brief Список уникальных узлов, дополняет класс RawCells
class RawNodes final {
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
    RawIncident incident;

    /// @brief Данные узлов
    Storage data;


    /// @brief Конструктор по умолчанию
    RawNodes() = default;

    /// @brief Пустое хранилище узлов?
    bool empty() const { return coord.empty(); }

    /// @brief Число уникальных узлов
    index_t n_nodes() const { return static_cast<index_t>(coord.size()); }

    /// @brief Размер списка инцидентных ячеек
    index_t n_incident() const { return incident.n_values(); }

    /// @brief Очистить массивы
    void clear();

    /// @brief Расширить массивы по числу узлов
    void resize(index_t n_nodes, index_t n_incident);

    /// @brief Увеличить размер под массив для AMR узлов
    void resize_amr(index_t n_nodes, int dim);

    /// @brief Расширить массивы по числу узлов
    void reserve(index_t n_nodes, index_t n_incident);

    /// @brief Увеличить размер под массив для AMR узлов
    void reserve_amr(index_t n_nodes, int dim);

    /// @brief Сжать буферы массивов под актуальные размеры
    void shrink_to_fit();

    /// @{ @name Топологические свойства узлов

    /// @brief Актуальный узел?
    bool is_actual(index_t inode) const { return index[inode] >= 0; }

    /// @brief Узел к удалению
    bool is_undefined(index_t inode) const { return index[inode] < 0; }

    /// @brief Устанавливает index = -1 (узел вне сетки)
    void set_undefined(index_t inode) { index[inode] = -1; }

    /// @}

    /// @brief Скопировать все данные с индекса from на индекс to.
    void copy_data(index_t from, index_t to);

    /// @brief Скопировать все данные целиком с индекса from,
    /// в хранилище dst на индекс to
    void copy_data(index_t from, RawNodes* dst, index_t to) const;

    void copy_geom(index_t ic, RawNodes& nodes,
                   index_t jc, index_t inc_offset) const;

    /// @brief Кортеж {vert.index, nodes}
    using Incomplete = std::tuple<RawVerts, RawNodes>;

    /// @brief Сгенерировать уникальные узлы для множества ячеек
    template <bool complete>
    static Incomplete generate(const RawCells& cells);

    /// @brief Установить уникальные узлы для множества ячеек
    void setup_for(RawCells& cells);

    /// @brief Расход памяти
    memory_t memory_usage() const;

    /// @brief Проверить согласованность размеров
    int check_sizes() const;

    /// @brief Проверка уникальных узлов для однопроцессорной версии
    int check_nodes(const RawCells& locals) const;

    /// @brief Проверка уникальных узлов в MPI версии
    int check_nodes(const RawCells& locals, const RawCells& ghosts, const RawNodes& ghost_nodes) const;
};

extern template RawNodes::Incomplete RawNodes::generate<true> (const RawCells& cells);
extern template RawNodes::Incomplete RawNodes::generate<false>(const RawCells& cells);

} // namespace zephyr::mesh