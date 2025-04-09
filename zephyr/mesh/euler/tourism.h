#pragma once
#include <vector>

#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/euler/amr_storage.h>

namespace zephyr::mesh {

using namespace zephyr::utils;
using zephyr::utils::mpi;

/// @brief Поля данных для отправки
enum class Post : int {
    FULL = 0,
    DATA = 2,
    FLAG = 1,
    // ...
};

#ifdef ZEPHYR_MPI

/// @brief Поддерживает согласованные обменные слои
class Tourism final {
public:
    Tourism();

    /// @brief Синхронизует размеры типов с основным хранилищем,
    /// также инициализирует MPI-типы для пересылок?
    void init_types(const AmrStorage& locals);

    /// @brief Сжать массивы до актуальных размеров
    void shrink_to_fit();

    mpi::datatype get_mpi_type() const;

    // Должен быть готов m_border, отправляет m_border -> aliens
    void send(const AmrStorage& locals, Post post = Post::FULL);
    
    // Получает aliens
    void recv(AmrStorage& aliens, Post post = Post::FULL);

    // Построить обменный слой (список border)
    void build_aliens(AmrStorage& locals, AmrStorage& aliens);

private:
    // Построить обменный слой (список border)
    void build_border(AmrStorage& locals);

    AmrStorage m_border;

    mpi::datatype m_item_mpi_type = nullptr;

    std::vector<int> m_count_to_send;
    std::vector<int> m_count_to_recv;

    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

    std::vector<std::vector<int>> m_border_indices;

    std::vector<MPI_Request> m_requests_send;
    std::vector<MPI_Request> m_requests_recv;
};

#endif

} // namespace zephyr::mesh