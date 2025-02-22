#pragma once
#include <vector>
#include <zephyr/mesh/euler/amr_storage.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

using zephyr::utils::mpi;

/// @class Просто забавное название)
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
    ~Tourism() = default;

    // Должно вызываться в начале build_aliens
    // Зануляет переменные, которые нужно занулить
    void reset() {
        m_recv_offsets[0] = 0;
        m_send_offsets[0] = 0;

        for(int r = 0; r<m_border_indices.size(); ++r)
            m_border_indices[r].clear();
    }

    // Должен быть готов: m_tourism
    // Отправляет m_tourism.m_border -> aliens
    void exchange_start(const AmrStorage& locals){
#ifdef ZEPHYR_MPI
        const int size = mpi::size();
        const int rank = mpi::rank();

        m_requests_send.resize(size);
        m_requests_recv.resize(size);
        std::fill(m_requests_send.begin(), m_requests_send.end(), MPI_REQUEST_NULL);    
        std::fill(m_requests_recv.begin(), m_requests_recv.end(), MPI_REQUEST_NULL);

        int temp_border_it = 0;
        for(auto& border_indices : m_border_indices){
            for(auto cell_index : border_indices){
                // [?]            
                memcpy(m_border[temp_border_it].ptr(), locals[cell_index].ptr(), locals.itemsize());
                ++temp_border_it;
            }
        }
        for (int i = 0; i < size; ++i) {
            if (i != rank && m_count_to_send[i] != 0) {            
                MPI_Isend(
                    m_border.item(0).ptr() + m_send_offsets[i], m_count_to_send[i], MPI_BYTE,
                    i, 0, mpi::comm(), &m_requests_send[i]
                );
            }
        }
#endif
    }

    void exchange_end(AmrStorage& aliens) {
#ifdef ZEPHYR_MPI
        const int size = mpi::size();
        const int rank = mpi::rank();

        for (int i = 0; i < size; ++i) {
            if (i != rank && m_count_to_recv[i] != 0) {
                MPI_Irecv(
                    aliens.item(0).ptr() + m_recv_offsets[i], m_count_to_recv[i], MPI_BYTE, 
                    i, 0, mpi::comm(), &m_requests_recv[i]
                );
            }
        }

        if(m_requests_send[rank] != MPI_REQUEST_NULL){
            MPI_Wait(&m_requests_send[rank], MPI_STATUSES_IGNORE);
            MPI_Wait(&m_requests_recv[rank], MPI_STATUSES_IGNORE);
        }

        MPI_Barrier(mpi::comm());
#endif
    }

    std::vector<int> m_count_to_send;
    std::vector<int> m_count_to_recv;

    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

#ifdef ZEPHYR_MPI 
    std::vector<MPI_Request> m_requests_send;
    std::vector<MPI_Request> m_requests_recv;
#endif

    std::vector<std::vector<int>> m_border_indices;

    AmrStorage m_border;

#ifdef ZEPHYR_MPI
    std::vector<MPI_Status> m_status;
#endif
};

} // namespace zephyr::mesh