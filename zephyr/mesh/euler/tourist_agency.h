#pragma once
#include <vector>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

using zephyr::utils::mpi;

/// @class Просто забавное название)
class TouristAgency {
public:
    TouristAgency() :
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

    // Должно вызываться в начале build_aliens
    void reset() {
        m_recv_offsets[0] = 0;
        m_send_offsets[0] = 0;

        for(int r = 0; r<m_border_indices.size(); ++r)
            m_border_indices[r].clear();
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