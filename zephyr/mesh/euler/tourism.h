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

    /// @brief Синхронизует типы с основным хранилищем
    void init_types(const AmrCells& locals);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Добавить тип данных в border/ghosts
    template <typename T>
    Storable<T> add(const std::string& name);

    /// @brief Добавить векторный тип данных в border/ghosts
    template<typename T>
    Storable<T> add(const std::string& name, int count);

    /// @brief Поменять местами два типа в хранилище
    template<typename T>
    void swap(Storable<T> var1, Storable<T> var2);

    /// @brief Построить обменные слои (border и ghosts).
    /// @param locals Локальное хранилище ячеек, на гранях должны быть корректно
    /// указаны индексы смежности (adjacent.rank, adjacent.index).
    /// В хранилище locals изменяются индексы смежности adjacent.ghost.
    void update(AmrCells& locals);

    /// @}

    /// @{ name Обменные операции

    /// @brief Синхронная передача данных
    /// @param vars Набор переменных типа Storable<T>
    template <typename... Args>
    void sync(const AmrCells& locals, Args&&... vars);

    /// @brief Синхронная передача сеточных данных
    template <MpiTag tag>
    void sync(const AmrCells& locals);

    /// @brief Перенести сеточные данные из locals в border_
    template <MpiTag tag>
    void prepare(const AmrCells& locals) {
        throw std::runtime_error("prepare<" + to_string(tag) + "> is not implemented");
    }

    /// @brief Асинхронная отправка сеточных данных
    template <MpiTag tag>
    Requests isend() {
        throw std::runtime_error("isend<" + to_string(tag) + ">() is not implemented");
    }

    /// @brief Асинхронное получение сеточных данных
    template <MpiTag tag>
    Requests irecv() {
        throw std::runtime_error("irecv<" + to_string(tag) + ">() is not implemented");
    }

    /// @}

    /// @{ @name get-функции

    /// @brief Ссылка на ghost-слой ячеек
    AmrCells& ghosts() { return ghosts_; }

    /// @brief Ссылка на ghost-слой ячеек
    const AmrCells& ghosts() const { return ghosts_; }

    /// @brief Ссылка на ghost-слой узлов
    AmrNodes& ghosts_nodes() { return ghosts_nodes_; }

    /// @brief Ссылка на ghost-слой узлов
    const AmrNodes& ghosts_nodes() const { return ghosts_nodes_; }

    /// @brief Ссылка на border-слой
    AmrCells& border() { return border_; }

    /// @brief Ссылка на border-слой
    const AmrCells& border() const { return border_; }

    /// @brief Маршрутизатор при обмене ячейками
    const Router& cell_router() const { return cell_router_; }

    /// @brief Индексы ячеек для отправки (из locals)
    const std::vector<index_t>& border_indices() const {
        return border_indices_;
    }

    /// @}

    /// @{ @name Специальные функции

    /// @brief Долго объяснять...
    // Выставляет next у border и ghosts (получает после отправки),
    // расширяет все массивы. Делает корректные router для пересылок,
    // но сами слои не заполняет, только resize.
    // Также выставляет border_indices_ для полностью адаптированной
    // сетки, это делается по массиву next внутри расширенного locals.
    // Портит index у border-ячеек, там кодируются дочерние ячейки.
    // На border слое должны быть предварительно выставлены флаги.
    // И в целом border-слой должен быть составлен правильно.
    template<int dim>
    void setup_positions(const std::vector<index_t>& locals_next);

    /// @brief Переслать геометрию ячеек locals -> ghosts
    void send_geometry(const AmrCells& locals);

    /// @brief Восстановить индексы adj.index для локальных ячеек
    void restore_indices(AmrCells& locals) const;

    /// @brief Изменить border хранилище под текущий Router
    void resize_border();

    /// @brief Изменить ghosts хранилище под текущий Router
    void resize_ghosts();

    /// @brief Расширить border хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_border();

    /// @brief Расширить ghosts хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_ghosts();

    /// @}

private:
    // Скопировать геометрию в border-ячейки
    void prepare_geometry(const AmrCells& locals);

    // Запаковать и отправить геометрию ячеек
    void sync_geometry();

    // Установить индекс NEXT у border-ячеек
    template<int dim>
    std::vector<index_t> setup_border_next();

    // Построить обменный border-слой, также подсчитывает количество пересылаемых
    // элементов, копирует геометрию ячеек из locals в созданный border слой.
    void build_border(const AmrCells& locals);

    // Считает количество примитивов для пересылки, заполняет величины
    // send_count во всех Router, но не recv_count
    void fill_send_count(const AmrCells& locals);

    // Заполняет массивы индексов пересылаемых ячеек
    void fill_indices(const AmrCells& locals);

    void pack_border_indices();

    void unpack_border_indices();

    void find_connections(AmrCells& locals, int rank) const;

    void unpack_ghost_indices();

    // ========================================================================
    //            Выставить финальные значения border_indices_
    // ========================================================================
    template<int dim>
    void update_border_indices(const std::vector<index_t>& locals_next);



    // Копирует данные из locals в border-слой
    template <typename T>
    void prepare(const AmrCells& locals, Storable<T> var);

    // Уникальные индексы border-ячеек по возрастанию
    // std::vector<index_t> unique_border_indices_;

    // Индексы ячеек, которые составляют хранилище border_
    std::vector<index_t> border_indices_;

    // Хранилище для ячеек на отправку. Ячейки, которые отправляются на один
    // процесс, располагаются сплошным блоком. Ячейка может быть включена
    // в массив дважды, если отправляется нескольким процессам.
    AmrCells border_;
    AmrCells ghosts_;

    // Маршрутизаторы для отправки примитивов из border
    Router cell_router_;
    Router face_router_;
    Router vert_router_;

    // Индексы уникальных узлов, которые составляют border_nodes.
    std::vector<index_t> border_nodes_indices_;

    // Хранилище уникальных узлов на отправку. Узлы, которые отправляются
    // на один процесс, располагаются сплошным блоком. Узел может быть включен
    // в массив дважды, если отправляется нескольким процессам.
    AmrNodes border_nodes_;
    AmrNodes ghosts_nodes_;

    // Маршрутизаторы для отправки примитивов из border_nodes
    Router node_router_;
    Router inct_router_;
};

// ============================================================================
//                        Реализации шаблонных функций
// ============================================================================

template<typename T>
Storable<T> Tourism::add(const std::string& name) {
    auto res1 = border_.data.add<T>(name);
    auto res2 = ghosts_.data.add<T>(name);
    if (res1 != res2) {
        throw std::runtime_error("EuMesh error: bad add<T> #1");
    }
    return res1;
}

template<typename T>
Storable<T> Tourism::add(const std::string& name, int count) {
    auto res1 = border_.data.add<T>(name, count);
    auto res2 = ghosts_.data.add<T>(name, count);
    if (res1 != res2) {
        throw std::runtime_error("EuMesh error: bad add<T> #1");
    }
    return res1;
}

template<typename T>
void Tourism::swap(Storable<T> var1, Storable<T> var2) {
    border_.data.swap<T>(var1, var2);
    ghosts_.data.swap<T>(var1, var2);
}

template <typename T>
void Tourism::prepare(const AmrCells& locals, Storable<T> var) {
    const utils::Buffer& src = locals .data[var];
          utils::Buffer& dst = border_.data[var];

    for (size_t ic = 0; ic < border_indices_.size(); ++ic) {
        src.copy_data(border_indices_[ic], dst, ic);
    }
}

template <typename... Args>
void Tourism::sync(const AmrCells& locals, Args&&... vars) {
    static_assert(sizeof...(Args) > 0, "Tourism::sync, zero arguments");
    soa::assert_storable<Args...>();
    
    // Отправить и дождаться одну переменную
    auto sync_one = [&](auto&& var) {
        prepare(locals, var);
        const utils::Buffer& src = border_.data[var];
        auto send_req = cell_router_.isend(src, static_cast<MpiTag>(var.tag()));

        utils::Buffer& dst = ghosts_.data[var];
        auto recv_req = cell_router_.irecv(dst, static_cast<MpiTag>(var.tag()));

        send_req.wait();
        recv_req.wait();
    };

    ( sync_one(vars), ... );
}




template <> inline
void Tourism::prepare<MpiTag::RANK>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_indices_.size(); ++ic) {
        border_.rank[ic] = locals.rank[border_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::NEXT>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_indices_.size(); ++ic) {
        border_.next[ic] = locals.next[border_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::INDEX>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_indices_.size(); ++ic) {
        border_.index[ic] = locals.index[border_indices_[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::FLAG>(const AmrCells& locals) {
    for (size_t ic = 0; ic < border_indices_.size(); ++ic) {
        border_.flag[ic] = locals.flag[border_indices_[ic]];
    }
}

template <> inline
Requests Tourism::isend<MpiTag::RANK>() {
    return cell_router_.isend(border_.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::NEXT>() {
    return cell_router_.isend(border_.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::isend<MpiTag::INDEX>() {
    return cell_router_.isend(border_.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::isend<MpiTag::FLAG>() {
    return cell_router_.isend(border_.flag, MpiTag::FLAG);
}

template <> inline
Requests Tourism::irecv<MpiTag::RANK>() {
    return cell_router_.irecv(ghosts_.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::NEXT>() {
    return cell_router_.irecv(ghosts_.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::irecv<MpiTag::INDEX>() {
    return cell_router_.irecv(ghosts_.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::irecv<MpiTag::FLAG>() {
    return cell_router_.irecv(ghosts_.flag, MpiTag::FLAG);
}

template <MpiTag tag>
void Tourism::sync(const AmrCells& locals) {
    prepare<tag>(locals);
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