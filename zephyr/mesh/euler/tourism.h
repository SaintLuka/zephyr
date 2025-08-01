#pragma once

#include <zephyr/configuration.h>

#ifdef ZEPHYR_MPI

#include <vector>
#include <thread>

#include <zephyr/mesh/euler/router.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh {

using zephyr::utils::Storable;

/// @brief Поддерживает согласованные обменные слои
class Tourism final {
public:
    /// @brief Тривиальный конструктор
    Tourism() = default;

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Синхронизует типы с основным хранилищем
    void init_types(const AmrCells& locals);

    /// @brief Добавить тип данных на border-слой ячеек
    template<typename T>
    auto add(const std::string& name) {
        return m_border.data.add<T>(name);
    }

    /// @brief Добавить тип данных на border-слой ячеек
    template<typename T>
    auto swap(Storable<T> var1, Storable<T> var2) {
        return m_border.data.swap<T>(var1, var2);
    }


    /// @brief Построить обменные слои, на гранях ячеек должны быть корректно
    /// заданы величины (adjacent.rank, adjacent.index)
    /// @details Сначала выполняет build_border, затем пересылает геометрию
    /// соседним процессам и формирует aliens-слой.
    void build_aliens(AmrCells& locals, AmrCells& aliens);

    /// @brief Неполное восстановление обменных слоев, корректно выставляются
    /// размеры слоев, матрицы пересылок, но геометрические данные копируются
    /// в border и aliens только частично.
    void build_aliens_basic(AmrCells& locals, AmrCells& aliens);


    /// @brief Синхронная передача данных
    /// @param vars Набор переменных типа Storable<T>
    template <typename... Args>
    void sync(const AmrCells& locals, AmrCells& aliens, Args&&... vars);

    /// @brief Синхронная передача сеточных данных
    template <MpiTag tag>
    void sync(const AmrCells& locals, AmrCells& aliens);


    // Копирует сеточные данные из locals в m_border
    template <MpiTag tag>
    void prepare(const AmrCells& locals) {
        throw std::runtime_error("prepare<" + to_string(tag) + "> is not implemented");
    }

    // Вызывает isend для пересылки сеточных данных
    template <MpiTag tag>
    Requests isend() {
        throw std::runtime_error("isend<" + to_string(tag) + ">() is not implemented");
    }

    // Вызывает irecv для пересылки сеточных данных
    template <MpiTag tag>
    Requests irecv(AmrCells& aliens) {
        throw std::runtime_error("irecv<" + to_string(tag) + ">() is not implemented");
    }

private:
    // Считает количество примитивов для пересылки, заполняет величины
    // send_count во всех Router, но не recv_count
    void fill_send_count(const AmrCells& locals);

    // Заполняет массивы индексов пересылаемых ячеек
    void fill_indices(const AmrCells& locals);

    // Построить обменный border-слой, также подсчитывает количество пересылаемых
    // элементов, копирует геометрию ячеек из locals в созданный border слой.
    void build_border(const AmrCells& locals);

    // Построить обменный border-слой, также подсчитывает количество пересылаемых
    // элементов, копирует геометрию ячеек из locals в созданный border слой.
    void build_border_basic(const AmrCells& locals);


    // Копирует данные из locals в m_border
    template <typename T>
    void prepare(const AmrCells& locals, Storable<T> var);

private:
    // Уникальные индексы border-ячеек по возрастанию
    std::vector<index_t> m_unique_border_indices;

    // Индексы ячеек, которые составляют хранилище m_border
    std::vector<index_t> m_border_indices;

    // Хранилилще для ячеек на отправку. Ячейки, которые отправляются на один
    // процесс, располагаются сплошным блоком. Ячейка может быть включена
    // в массив дважды, если отправляется нескольким процессам.
    AmrCells m_border;

    // Маршрутизаторы для отправки примитивов из m_border
    Router m_cell_route;
    Router m_face_route;
    Router m_node_route;
};



template <typename T>
void Tourism::prepare(const AmrCells& locals, Storable<T> var) {
    auto src = locals  .data[var];
    auto dst = m_border.data[var];

    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        std::memcpy(dst + ic, src + m_border_indices[ic], sizeof(T));
    }
}

template <typename... Args>
void Tourism::sync(const AmrCells& locals, AmrCells& aliens, Args&&... vars) {
    static_assert(sizeof...(Args) > 0, "Tourism::sync, zero arguments");
    utils::assert_all_storable<Args...>();
    
    // Отправить и дождаться одну переменную
    auto sync_one = [&](auto&& var) {
        prepare(locals, var);
        const auto *src = m_border.data[var];
        auto send_req = m_cell_route.isend(src, MpiTag(var.tag()));

        auto *dst = aliens.data[var];
        auto recv_req = m_cell_route.irecv(dst, MpiTag(var.tag()));

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
    return m_cell_route.isend(m_border.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::NEXT>() {
    return m_cell_route.isend(m_border.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::isend<MpiTag::INDEX>() {
    return m_cell_route.isend(m_border.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::isend<MpiTag::FLAG>() {
    return m_cell_route.isend(m_border.flag, MpiTag::FLAG);
}

template <> inline
Requests Tourism::irecv<MpiTag::RANK>(AmrCells& aliens) {
    return m_cell_route.irecv(aliens.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::NEXT>(AmrCells& aliens) {
    return m_cell_route.irecv(aliens.next, MpiTag::NEXT);
}
template <> inline
Requests Tourism::irecv<MpiTag::INDEX>(AmrCells& aliens) {
    return m_cell_route.irecv(aliens.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::irecv<MpiTag::FLAG>(AmrCells& aliens) {
    return m_cell_route.irecv(aliens.flag, MpiTag::FLAG);
}

template <MpiTag tag>
void Tourism::sync(const AmrCells& locals, AmrCells& aliens) {
    prepare<tag>(locals);
    auto send_req = isend<tag>();
    auto recv_req = irecv<tag>(aliens);

    send_req.wait();
    recv_req.wait();
}

} // namespace zephyr::mesh

#endif