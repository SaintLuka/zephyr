#include "tourism.h"

namespace zephyr::mesh {

void Tourism::reset() {
    m_recv_offsets[0] = 0;
    m_send_offsets[0] = 0;

    for(int r = 0; r<m_border_indices.size(); ++r)
        m_border_indices[r].clear();
}

void Tourism::send(const AmrStorage& locals, Post post) {
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

void Tourism::recv(AmrStorage& aliens, Post post) {
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

} // namespace zephyr::mesh