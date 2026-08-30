#pragma once

#include <zephyr/configuration.h>

#ifdef ZEPHYR_MPI

#include <vector>

#include <zephyr/mesh/euler/router.h>
#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/mesh/euler/amr_nodes.h>

namespace zephyr::mesh {

/// @brief Поддерживает согласованные обменные слои
class Tourism final {
public:
    /// @{ @name Инициализация обменных слоёв

    /// @brief Тривиальный конструктор
    Tourism() = default;

    /// @brief Синхронизует тип сетки и данных с основным хранилищем
    void init_types(const AmrCells& cells);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Добавить тип данных в border/ghosts
    template <typename T>
    Storable<T> add_cell_data(const std::string& name);

    /// @brief Добавить векторный тип данных в border/ghosts
    template<typename T>
    Storable<T> add_cell_data(const std::string& name, int count);

    /// @brief Поменять местами два типа в хранилище
    template<typename T>
    void swap_cell_data(Storable<T> var1, Storable<T> var2);

    /// @brief Построить обменные слои (border и ghosts).
    /// @param cells
    /// @param nodes Локальное хранилище уникальных узлов, в списке
    /// incident должны быть корректно указаны индексы смежности
    /// (incident.rank, incident.index). В local_nodes изменяются индексы
    /// incident.ghost.
    ///
    /// Супер важное замечание, в cells.verts.index должны быть индексы
    /// только локальные! Никаких глобальных.

    void update(AmrCells& cells, AmrNodes& nodes);

    /// @}

    /// @{ name Обменные операции

    /// @brief Отправка/получение сеточных данных
    /// @param cells Локальное хранилище ячеек
    /// @param vars Набор переменных типа Storable<T>
    template <typename... Args>
    void sync(const AmrCells& cells, Args&&... vars);

    /// @brief Отправка/получение сеточных данных
    /// @param cells Локальное хранилище ячеек
    template <MpiTag tag>
    void sync(const AmrCells& cells);

    /// @brief Отправка/получение сеточных данных
    /// @param nodes Локальное хранилище узлов
    template <MpiTag tag>
    void sync(const AmrNodes& nodes);

    /// @brief Скопировать данные из cells в border_cells_ (индексы, геометрия, ...)
    template <MpiTag tag>
    void prepare(const AmrCells& cells) {
        throw std::runtime_error("prepare<" + to_string(tag) + "> is not implemented");
    }

    /// @brief Скопировать данные из nodes в border_nodes_ (индексы, геометрия, ...)
    template <MpiTag tag>
    void prepare(const AmrNodes& nodes) {
        throw std::runtime_error("prepare<" + to_string(tag) + "> is not implemented");
    }

    /// @brief Асинхронная отправка сеточных данных (индексы, геометрия, ...)
    template <MpiTag tag>
    Requests isend() {
        throw std::runtime_error("isend<" + to_string(tag) + ">() is not implemented");
    }

    /// @brief Асинхронное получение сеточных данных (индексы, геометрия, ...)
    template <MpiTag tag>
    Requests irecv() {
        throw std::runtime_error("irecv<" + to_string(tag) + ">() is not implemented");
    }

    /// @}

    /// @{ @name get-функции

    /// @brief Ссылка на ghost-слой ячеек
    AmrCells& ghost_cells() { return ghost_cells_; }

    /// @brief Ссылка на ghost-слой ячеек
    const AmrCells& ghost_cells() const { return ghost_cells_; }

    /// @brief Ссылка на ghost-слой узлов
    AmrNodes& ghost_nodes() { return ghost_nodes_; }

    /// @brief Ссылка на ghost-слой узлов
    const AmrNodes& ghost_nodes() const { return ghost_nodes_; }

    /// @brief Ссылка на border-слой ячеек
    AmrCells& border_cells() { return border_cells_; }

    /// @brief Ссылка на border-слой ячеек
    const AmrCells& border_cells() const { return border_cells_; }

    /// @brief Ссылка на border-слой ячеек
    AmrNodes& border_nodes() { return border_nodes_; }

    /// @brief Ссылка на border-слой ячеек
    const AmrNodes& border_nodes() const { return border_nodes_; }

    /// @brief Маршрутизатор при обмене ячейками
    const Router& cell_router() const { return cell_router_; }

    /// @brief Индексы ячеек для отправки (из locals)
    const std::vector<index_t>& border_indices() const {
        return border_cells_indices_;
    }

    /// @}

    /// @{ @name Специальные функции

    /// @brief Долго объяснять...
    // Выставляет next у border и ghosts (получает после отправки),
    // расширяет все массивы. Делает корректные router для пересылок,
    // но сами слои не заполняет, только resize.
    // Также выставляет border_cells_indices_ для полностью адаптированной
    // сетки, это делается по массиву next внутри расширенного locals.
    // Портит index у border-ячеек, там кодируются дочерние ячейки.
    // На border слое должны быть предварительно выставлены флаги.
    // И в целом border-слой должен быть составлен правильно.
    template<int dim>
    void setup_positions(const std::vector<index_t>& locals_next);

    /// @brief Переслать геометрию ячеек locals -> ghosts
    void send_geometry(const AmrCells& cells);

    /// @brief Восстановить индексы adj.index для локальных ячеек
    void restore_indices(AmrCells& cells) const;

    /// @brief Изменить border хранилище под текущий Router
    void resize_border_cells();
    void resize_border_nodes();

    /// @brief Расширить border хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_border_cells();
    void extend_border_nodes();

    /// @brief Изменить ghosts хранилище под текущий Router
    void resize_ghost_cells();
    void resize_ghost_nodes();

    /// @brief Расширить ghosts хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_ghost_cells();
    void extend_ghost_nodes();

    /// @}

private:
    // Копирует данные из local_cells_ в border_cells_
    template <typename T>
    void prepare(const AmrCells& cells, Storable<T> var);

    // ------------------------------------------------------------------------
    //                             Функция update()
    // ------------------------------------------------------------------------

    // Создает border-слой узлов, считает количество пересылаемых элементов,
    // копирует геометрию узлов из nodes в созданный border-слой узлов.
    void build_border_nodes(const AmrNodes& nodes);

    // Создает border-слой ячеек, считает количество пересылаемых элементов,
    // копирует геометрию ячеек из cells в созданный border-слой ячеек.
    void build_border_cells(const AmrCells& cells, const AmrNodes& nodes);

    // Считает количество примитивов для пересылки, заполняет величины
    // send_count в роутерах для узлов.
    // Требует выставленные nodes.incident.rank.
    void fill_nodes_send_count(const AmrNodes& nodes);

    // Заполняет массивы индексов пересылаемых узлов.
    // Соответствует функции fill_nodes_send_count(nodes).
    void fill_border_nodes_indices(const AmrNodes& nodes);

    // Считает количество примитивов для пересылки, заполняет величины
    // send_count в роутерах для ячеек. Рассматривается окрестность Неймана.
    // Требует выставленные faces.adjacent.rank.
    void fill_cells_send_count(const AmrCells& cells);

    // Заполняет массивы индексов пересылаемых ячеек.
    // Соответствует функции fill_cells_send_count(cells).
    void fill_border_cells_indices(const AmrCells& cells);

    // Считает количество примитивов для пересылки, заполняет величины
    // send_count в роутерах для ячеек. Рассматривается окрестность Мура.
    // Требует построенные ghost_nodes_, в них incident.rank, также
    // у cells.verts должны быть выставлены индексы (ghost, index).
    void fill_cells_send_count(const AmrCells& cells, const AmrNodes& nodes);

    // Заполняет массивы индексов пересылаемых ячеек.
    // Соответствует функции fill_cells_send_count(cells, nodes).
    void fill_border_cells_indices(const AmrCells& cells, const AmrNodes& nodes);

    // Скопировать геометрию из nodes в border-узлы
    void prepare_nodes_geometry(const AmrNodes& nodes);

    // Скопировать геометрию из cells в border-ячейки
    void prepare_cells_geometry(const AmrCells& cells);

    // Запаковать и отправить геометрию узлов
    void sync_nodes_geometry();

    // Запаковать и отправить геометрию ячеек
    void sync_cells_geometry();

    // Построить с нуля border-слой, ghost-слой и выполнить обмен
    void build_ghost_nodes(const AmrNodes &nodes);

    // Построить с нуля border-слой, ghost-слой и выполнить обмен
    void build_ghost_cells(const AmrCells& cells, const AmrNodes& nodes);

    // Найти индекс ячейки в массиве ghost_cells_ со значениями (rank, index).
    // Предполагается, что ghost_cells_ упорядочены по rank и index.
    index_t find_ghost_cell(int rank, index_t index) const;

    // Найти индекс узла в массиве ghost_nodes_ со значениями (rank, index).
    // Предполагается, что ghost_nodes_ упорядочены по rank и index.
    index_t find_ghost_node(int rank, index_t index) const;

    // Найти окрестность Неймана. Выставляет на гранях индекс adjacent.ghost
    // для ячеек, у которых есть ghost-сосед.
    void find_connections_neumann(AmrCells& cells, int rank) const;

    // Выставляет у дублирующихся вершин правильные индексы verts.ghost
    void find_connections_verts(AmrVerts &verts, int rank) const;

    // Найти окрестность Мура. Выставляет индексы incident.ghost у узлов и
    // adjacent.ghost у ячеек. В отличие от окрестности Неймана выставляет
    // данные индексы также у узлов и ячеек, которые сами находятся
    // в ghost-слоях. Это позволяет в дальнейшем переходить от ghost-ячеек
    // ко всем соседям, которые есть на данном процессе.
    void find_connections_moore(AmrCells& cells, AmrNodes& nodes, int rank);

    // ------------------------------------------------------------------------
    //                          AMR-приколы
    // ------------------------------------------------------------------------

    // Установить индекс NEXT у border-ячеек
    template <int dim>
    std::vector<index_t> setup_border_next();

    // Выставить финальные значения border_cells_indices_
    template <int dim>
    void update_border_indices(const std::vector<index_t>& locals_next);

private:
    // Уникальные индексы border-ячеек по возрастанию
    // std::vector<index_t> unique_border_indices_;

    /// Индексы ячеек, которые составляют хранилище border_cells_.
    std::vector<index_t> border_cells_indices_;

    // Хранилище для ячеек на отправку. Ячейки, которые отправляются на один
    // процесс, располагаются сплошным блоком. Ячейка может быть включена
    // в массив дважды, если отправляется нескольким процессам.
    AmrCells border_cells_;
    AmrCells ghost_cells_;

    // Маршрутизаторы для отправки примитивов из border_cells_.
    Router cell_router_;
    Router face_router_;
    Router vert_router_;

    // Индексы уникальных узлов, которые составляют border_nodes_.
    std::vector<index_t> border_nodes_indices_;

    // Хранилище уникальных узлов на отправку. Узлы, которые отправляются
    // на один процесс, располагаются сплошным блоком. Узел может быть включен
    // в массив дважды, если отправляется нескольким процессам.
    AmrNodes border_nodes_;

    // Хранилище для уникальных узлов, получаемых с других процессов.
    // Узлы упорядочены по рангу, а затем по возрастанию index (всегда ?).
    AmrNodes ghost_nodes_;

    // Маршрутизаторы для отправки примитивов из border_nodes
    Router node_router_;
    Router inct_router_;
};

// ============================================================================
//                        Реализации шаблонных функций
// ============================================================================

template<typename T>
Storable<T> Tourism::add_cell_data(const std::string& name) {
    auto res1 = border_cells_.data.add<T>(name);
    auto res2 = ghost_cells_.data.add<T>(name);
    if (res1 != res2) {
        throw std::runtime_error("Tourism error: different types in border and ghost cells #1");
    }
    return res1;
}

template<typename T>
Storable<T> Tourism::add_cell_data(const std::string& name, int count) {
    auto res1 = border_cells_.data.add<T>(name, count);
    auto res2 = ghost_cells_.data.add<T>(name, count);
    if (res1 != res2) {
        throw std::runtime_error("Tourism error: different types in border and ghost cells #2");
    }
    return res1;
}

template<typename T>
void Tourism::swap_cell_data(Storable<T> var1, Storable<T> var2) {
    border_cells_.data.swap<T>(var1, var2);
    ghost_cells_.data.swap<T>(var1, var2);
}

template <typename T>
void Tourism::prepare(const AmrCells& cells, Storable<T> var) {
    const utils::Buffer& src = cells.data[var];
          utils::Buffer& dst = border_cells_.data[var];

    for (size_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        src.copy_data(border_cells_indices_[ic], dst, ic);
    }
}

template <typename... Args>
void Tourism::sync(const AmrCells& cells, Args&&... vars) {
    static_assert(sizeof...(Args) > 0, "Tourism::sync, zero arguments");
    soa::assert_storable<Args...>();
    
    // Отправить и дождаться одну переменную
    auto sync_one = [&](auto&& var) {
        prepare(cells, var);
        const utils::Buffer& src = border_cells_.data[var];
        auto send_req = cell_router_.isend(src, static_cast<MpiTag>(var.tag()));

        utils::Buffer& dst = ghost_cells_.data[var];
        auto recv_req = cell_router_.irecv(dst, static_cast<MpiTag>(var.tag()));

        send_req.wait();
        recv_req.wait();
    };

    ( sync_one(vars), ... );
}




template <> inline
void Tourism::prepare<MpiTag::RANK>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        border_cells_.rank[ic] = locals.rank[border_cells_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::NEXT>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        border_cells_.next[ic] = locals.next[border_cells_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::INDEX>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        border_cells_.index[ic] = locals.index[border_cells_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::FLAG>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        border_cells_.flag[ic] = locals.flag[border_cells_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::NODE_RANK>(const AmrNodes& locals) {
    for (size_t in = 0; in < border_nodes_indices_.size(); ++in) {
        border_nodes_.rank[in] = locals.rank[border_nodes_indices_[in]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::NODE_INDEX>(const AmrNodes& locals) {
    for (size_t in = 0; in < border_nodes_indices_.size(); ++in) {
        border_nodes_.index[in] = locals.index[border_nodes_indices_[in]];
    }
}

template <> inline
Requests Tourism::isend<MpiTag::RANK>() {
    return cell_router_.isend(border_cells_.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::NEXT>() {
    return cell_router_.isend(border_cells_.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::isend<MpiTag::INDEX>() {
    return cell_router_.isend(border_cells_.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::isend<MpiTag::FLAG>() {
    return cell_router_.isend(border_cells_.flag, MpiTag::FLAG);
}
template <> inline
Requests Tourism::isend<MpiTag::NODE_RANK>() {
    return node_router_.isend(border_nodes_.rank, MpiTag::NODE_RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::NODE_INDEX>() {
    return node_router_.isend(border_nodes_.index, MpiTag::NODE_INDEX);
}

template <> inline
Requests Tourism::irecv<MpiTag::RANK>() {
    return cell_router_.irecv(ghost_cells_.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::NEXT>() {
    return cell_router_.irecv(ghost_cells_.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::irecv<MpiTag::INDEX>() {
    return cell_router_.irecv(ghost_cells_.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::irecv<MpiTag::FLAG>() {
    return cell_router_.irecv(ghost_cells_.flag, MpiTag::FLAG);
}
template <> inline
Requests Tourism::irecv<MpiTag::NODE_RANK>() {
    return node_router_.irecv(ghost_nodes_.rank, MpiTag::NODE_RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::NODE_INDEX>() {
    return node_router_.irecv(ghost_nodes_.index, MpiTag::NODE_INDEX);
}

template <MpiTag tag>
void Tourism::sync(const AmrCells& cells) {
    prepare<tag>(cells);
    auto send_req = isend<tag>();
    auto recv_req = irecv<tag>();

    send_req.wait();
    recv_req.wait();
}

template <MpiTag tag>
void Tourism::sync(const AmrNodes& nodes) {
    prepare<tag>(nodes);
    auto send_req = isend<tag>();
    auto recv_req = irecv<tag>();

    send_req.wait();
    recv_req.wait();
}

extern template std::vector<index_t> Tourism::setup_border_next<2>();
extern template std::vector<index_t> Tourism::setup_border_next<3>();

extern template void Tourism::setup_positions<2>(const std::vector<index_t>&);
extern template void Tourism::setup_positions<3>(const std::vector<index_t>&);

extern template void Tourism::update_border_indices<2>(const std::vector<index_t>&);
extern template void Tourism::update_border_indices<3>(const std::vector<index_t>&);

} // namespace zephyr::mesh

#endif