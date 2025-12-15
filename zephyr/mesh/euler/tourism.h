#pragma once

#include <zephyr/configuration.h>

#ifdef ZEPHYR_MPI

#include <vector>

#include <zephyr/mesh/euler/router.h>
#include <zephyr/mesh/euler/amr_cells.h>

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

    /// @brief Добавить тип данных в border/aliens
    template<typename T>
    auto add(const std::string& name);

    /// @brief Поменять местами два типа в хранилище
    template<typename T>
    void swap(Storable<T> var1, Storable<T> var2);

    /// @brief Построить обменные слои (border и aliens).
    /// @param locals Локальное хранилище ячеек, на гранях должны быть корректно
    /// указаны индексы смежности (adjacent.rank, adjacent.index).
    /// В хранилище locals изменяются индексы смежности adjacent.alien.
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

    /// @brief Перенести сеточные данные из locals в m_border
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

    /// @brief Ссылка на alien-слой
    AmrCells& aliens() { return m_aliens; }

    /// @brief Ссылка на alien-слой
    const AmrCells& aliens() const { return m_aliens; }

    /// @brief Ссылка на border-слой
    AmrCells& border() { return m_border; }

    /// @brief Ссылка на border-слой
    const AmrCells& borders() const { return m_border; }

    /// @brief Маршрутизатор при обмене ячейками
    const Router& cell_router() const { return m_cell_router; }

    /// @brief Индексы ячеек для отправки (из locals)
    const std::vector<index_t>& border_indices() const {
        return m_border_indices;
    }

    /// @}

    /// @{ @name Специальные функции

    /// @brief Долго объяснять...
    // Выставляет next у border и aliens (получает после отправки),
    // расширяет все массивы. Делает корректные router для пересылок,
    // но сами слои не заполняет, только resize.
    // Также выставляет m_border_indices для полностью адаптированной
    // сетки, это делается по массиву next внутри расширенного locals.
    // Портит index у border-ячеек, там кодируются дочерние ячейки.
    // На border слое должны быть предварительно выставлены флаги.
    // И в целом border-слой должен быть составлен правильно.
    template<int dim>
    void setup_positions(const std::vector<index_t>& locals_next);

    /// @brief Переслать геометрию ячеек locals -> aliens
    void send_geometry(const AmrCells& locals);

    /// @brief Восстановить индексы adj.index для локальных ячеек
    void restore_indices(AmrCells& locals) const;

    /// @brief Изменить border хранилище под текущий Router
    void resize_border();

    /// @brief Изменить aliens хранилище под текущий Router
    void resize_aliens();

    /// @brief Расширить border хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_border();

    /// @brief Расширить aliens хранилище под текущий Router
    /// (может только увеличить размеры)
    void extend_aliens();

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

    void unpack_aliens_indices();

    // ========================================================================
    //            Выставить финальные значения m_border_indices
    // ========================================================================
    template<int dim>
    void update_border_indices(const std::vector<index_t>& locals_next);



    // Копирует данные из locals в border-слой
    template <typename T>
    void prepare(const AmrCells& locals, Storable<T> var);

    // Уникальные индексы border-ячеек по возрастанию
    // std::vector<index_t> m_unique_border_indices;

    // Индексы ячеек, которые составляют хранилище m_border
    std::vector<index_t> m_border_indices;

    // Хранилище для ячеек на отправку. Ячейки, которые отправляются на один
    // процесс, располагаются сплошным блоком. Ячейка может быть включена
    // в массив дважды, если отправляется нескольким процессам.
    AmrCells m_border;
    AmrCells m_aliens;

    // Маршрутизаторы для отправки примитивов из m_border
    Router m_cell_router;
    Router m_face_router;
    Router m_node_router;
};

// ============================================================================
//                        Реализации шаблонных функций
// ============================================================================

/// @brief Добавить тип данных в border/aliens
template<typename T>
auto Tourism::add(const std::string& name) {
    auto res1 = m_border.data.add<T>(name);
    auto res2 = m_aliens.data.add<T>(name);
    if (res1 != res2) {
        throw std::runtime_error("EuMesh error: bad add<T> #1");
    }
    return res1;
}

/// @brief Поменять местами два типа в хранилище
template<typename T>
void Tourism::swap(Storable<T> var1, Storable<T> var2) {
    m_border.data.swap<T>(var1, var2);
    m_aliens.data.swap<T>(var1, var2);
}

template <typename T>
void Tourism::prepare(const AmrCells& locals, Storable<T> var) {
    const utils::Buffer& src = locals  .data[var];
          utils::Buffer& dst = m_border.data[var];

    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        src.copy_data(m_border_indices[ic], dst, ic);
    }
}

template <typename... Args>
void Tourism::sync(const AmrCells& locals, Args&&... vars) {
    static_assert(sizeof...(Args) > 0, "Tourism::sync, zero arguments");
    soa::assert_storable<Args...>();
    
    // Отправить и дождаться одну переменную
    auto sync_one = [&](auto&& var) {
        prepare(locals, var);
        const utils::Buffer& src = m_border.data[var];
        auto send_req = m_cell_router.isend(src, static_cast<MpiTag>(var.tag()));

        utils::Buffer& dst = m_aliens.data[var];
        auto recv_req = m_cell_router.irecv(dst, static_cast<MpiTag>(var.tag()));

        send_req.wait();
        recv_req.wait();
    };

    ( sync_one(vars), ... );
}




template <> inline
void Tourism::prepare<MpiTag::RANK>(const AmrCells& locals) {
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        m_border.rank[ic] = locals.rank[m_border_indices[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::NEXT>(const AmrCells& locals) {
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        m_border.next[ic] = locals.next[m_border_indices[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::INDEX>(const AmrCells& locals) {
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        m_border.index[ic] = locals.index[m_border_indices[ic]];
    }
}

template <> inline
void Tourism::prepare<MpiTag::FLAG>(const AmrCells& locals) {
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        m_border.flag[ic] = locals.flag[m_border_indices[ic]];
    }
}

template <> inline
Requests Tourism::isend<MpiTag::RANK>() {
    return m_cell_router.isend(m_border.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::NEXT>() {
    return m_cell_router.isend(m_border.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::isend<MpiTag::INDEX>() {
    return m_cell_router.isend(m_border.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::isend<MpiTag::FLAG>() {
    return m_cell_router.isend(m_border.flag, MpiTag::FLAG);
}

template <> inline
Requests Tourism::irecv<MpiTag::RANK>() {
    return m_cell_router.irecv(m_aliens.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::NEXT>() {
    return m_cell_router.irecv(m_aliens.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::irecv<MpiTag::INDEX>() {
    return m_cell_router.irecv(m_aliens.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::irecv<MpiTag::FLAG>() {
    return m_cell_router.irecv(m_aliens.flag, MpiTag::FLAG);
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