#pragma once
#include <vector>
#include <zephyr/mesh/euler/amr_storage.h>
#include <zephyr/utils/mpi.h>

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

class Tourism {
public:
    Tourism() :
#ifdef ZEPHYR_MPI
        m_requests_recv(mpi::size()),
        m_requests_send(mpi::size()),
#endif
	    m_count_to_send(mpi::size()),
	    m_count_to_recv(mpi::size()),
	    m_send_offsets(mpi::size(), 0),
	    m_recv_offsets(mpi::size(), 0),
        m_border_indices(mpi::size(), std::vector<int>())
    {}
    ~Tourism();

    void init_mpi_type(int size);
    MPI_Datatype get_mpi_type() const;

    // Должно вызываться в начале build_aliens
    // Зануляет переменные, которые нужно занулить
    void reset();

    // Должен быть готов: m_tourism
    // Отправляет m_tourism.m_border -> aliens
    void send(const AmrStorage& locals, Post post = Post::FULL);
    
    // Получает aliens
    void recv(AmrStorage& aliens, Post post = Post::FULL);

    // 
    void build_border(AmrStorage& locals, AmrStorage& aliens);

    MPI_Datatype m_item_mpi_type = nullptr;
private:
    std::vector<int> m_count_to_send;
    std::vector<int> m_count_to_recv;

    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

    std::vector<std::vector<int>> m_border_indices;

    AmrStorage m_border;
#ifdef ZEPHYR_MPI 
    std::vector<MPI_Request> m_requests_send;
    std::vector<MPI_Request> m_requests_recv;
    std::vector<MPI_Status> m_status;

#endif
};

} // namespace zephyr::mesh