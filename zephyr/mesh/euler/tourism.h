#pragma once

#include <vector>
#include <thread>
#include <future>

#include <zephyr/mesh/euler/router.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh {

using namespace zephyr::utils;
using zephyr::utils::mpi;

#ifdef ZEPHYR_MPI

/// @brief Поддерживает согласованные обменные слои
class Tourism final {
public:
    Tourism();

    /// @brief Синхронизует типы с основным хранилищем
    void init_types(const AmrCells& locals);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    /// @brief Должен быть готов m_border, осуществляет копирование данных
    /// из хранилища locals в массив m_border.
    template <typename T>
    void prepare(const AmrCells& locals, Storable<T> val);

    /// @brief Аналогично предыдущему, но запускается асинхронно
    template <typename T>
    std::future<void> iprepare(const AmrCells& locals, Storable<T> val);

    /// @brief Перенести данные из locals в border-слой
    template <MpiTag tag>
    void prepare(const AmrCells& locals) {
        throw std::runtime_error("prepare<tag>() is not implemented");
    }

    template <MpiTag tag>
    std::future<void> iprepare(const AmrCells& locals);

    /// @brief Асинхронная отправка, требует завершения prepare,
    /// то есть locals должны быть уже перенесены в m_border
    template <typename T>
    Requests isend(Storable<T> val);

    /// @brief Специализации isend для пересылки сеточных данных
    template <MpiTag tag>
    Requests isend() { throw std::runtime_error("isend<tag>() is not implemented"); }


    /// @brief Асинхронное получение
    template <typename T>
    Requests irecv(AmrCells& aliens, Storable<T> val);

    /// @brief Специализации irecv для пересылки сеточных данных
    template <MpiTag tag>
    Requests irecv(AmrCells& aliens) {
        throw std::runtime_error("irecv<tag>() is not implemented");
    }


    /// @brief Отправить и получить (синхронная операция)
    template <typename T>
    void sync(const AmrCells& locals, AmrCells& aliens, Storable<T> val);


    // Построить обменные слои, сначала выполняет build_border, затем пересылает
    // геометрию различным процессам и формирует alien-слои.
    void build_aliens(AmrCells& locals, AmrCells& aliens);

public:
    // Построить обменный border-слой, также подсчитывает количество пересылаемых
    // элементов, копирует геометрию ячеек из locals в созданный border слой.
    void build_border(AmrCells& locals);


    std::vector<std::vector<int>> m_border_indices;

    AmrCells m_border;

    Router cell_route;
    Router face_route;
    Router node_route;
};

template <typename T>
void Tourism::prepare(const AmrCells& locals, Storable<T> val) {
    const T* src = locals.data(val);
    T* dst = m_border.data(val);

    int counter = 0;
    for(auto& indices: m_border_indices){
        for (index_t ic: indices) {
            std::memcpy(dst + counter, src + ic, sizeof(T));
            ++counter;
        }
    }
}

template <typename T>
std::future<void> Tourism::iprepare(const AmrCells& locals, Storable<T> val) {
    return std::async(std::launch::async,
                      [this, &locals, val]() { this->prepare(locals, val); });
}

template <> inline
void Tourism::prepare<MpiTag::RANK>(const AmrCells& locals) {
    int counter = 0;
    for(auto& indices: m_border_indices){
        for (index_t ic: indices) {
            m_border.rank[counter] = locals.rank[ic];
            ++counter;
        }
    }
}

template <> inline
void Tourism::prepare<MpiTag::INDEX>(const AmrCells& locals) {
    int counter = 0;
    for(auto& indices: m_border_indices){
        for (index_t ic: indices) {
            m_border.index[counter] = locals.index[ic];
            ++counter;
        }
    }
}

template <> inline
void Tourism::prepare<MpiTag::FLAG>(const AmrCells& locals) {
    int counter = 0;
    for(auto& indices: m_border_indices){
        for (index_t ic: indices) {
            m_border.flag[counter] = locals.flag[ic];
            ++counter;
        }
    }
}

template <MpiTag tag> inline
std::future<void> Tourism::iprepare(const AmrCells& locals) {
    return std::async(std::launch::async,
                      [this, &locals]() { this->prepare<tag>(locals); });
}

template <typename T>
Requests Tourism::isend(Storable<T> val) {
    const T* src = m_border.data(val);
    return cell_route.isend(src, MpiTag(val.tag()));
}

template <> inline
Requests Tourism::isend<MpiTag::RANK>() {
    return cell_route.isend(m_border.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::isend<MpiTag::INDEX>() {
    return cell_route.isend(m_border.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::isend<MpiTag::FLAG>() {
    return cell_route.isend(m_border.flag, MpiTag::FLAG);
}

template <typename T>
Requests Tourism::irecv(AmrCells& aliens, Storable<T> val) {
    T* dst = aliens.data(val);
    return cell_route.irecv(dst, MpiTag(val.tag()));
}

template <> inline
Requests Tourism::irecv<MpiTag::RANK>(AmrCells& aliens) {
    return cell_route.irecv(aliens.rank, MpiTag::RANK);
}
template <> inline
Requests Tourism::irecv<MpiTag::INDEX>(AmrCells& aliens) {
    return cell_route.irecv(aliens.index, MpiTag::INDEX);
}
template <> inline
Requests Tourism::irecv<MpiTag::FLAG>(AmrCells& aliens) {
    return cell_route.irecv(aliens.flag, MpiTag::FLAG);
}

template <typename T>
void Tourism::sync(const AmrCells& locals, AmrCells& aliens, Storable<T> val) {
    prepare(locals, val);
    auto data_send = isend(val);
    auto data_recv = irecv(aliens, val);

    data_send.wait();
    data_recv.wait();
}

#endif

} // namespace zephyr::mesh