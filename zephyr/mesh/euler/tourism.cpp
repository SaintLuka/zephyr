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

    m_aliens.resize(m_tourism.m_recv_offsets[size - 1] + m_tourism.m_count_to_recv[size - 1]);

	for(int i = 0; i < size; ++i){
		m_tourism.m_count_to_send[i] *= m_aliens.itemsize();
		m_tourism.m_count_to_recv[i] *= m_aliens.itemsize();
		m_tourism.m_send_offsets[i] *= m_aliens.itemsize();
		m_tourism.m_recv_offsets[i] *= m_aliens.itemsize();
	}

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

void Tourism::build_border(){
#ifdef ZEPHYR_MPI
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Заполняем m_border_indices
    for (auto cell: *this) {
        // build alien можно вызвать для не совсем нормальной сетки
        if (cell.geom().is_undefined()) {
            continue;
        }
        for (auto face: cell.faces()) {
            //if(mpi::master())
            //	printf("face.adjacent().rank: %d\n", face.adjacent().rank);
            auto& border_indices = m_tourism.m_border_indices[face.adjacent().rank];
            if(face.adjacent().rank != rank && (border_indices.empty() || border_indices.back() != cell.geom().index))
                border_indices.push_back(cell.geom().index);
        }
    }

    // Заполняем m_count_to_send
    for(int r = 0; r < size; ++r)
        m_tourism.m_count_to_send[r] = m_tourism.m_border_indices[r].size();
    // Заполняем m_send_offsets
    for(int r = 1; r < size; ++r){
        m_tourism.m_send_offsets[r] = m_tourism.m_send_offsets[r - 1] + m_tourism.m_count_to_send[r - 1];
    }

    // Отправляем m_count_to_send -> получаем m_count_to_recv
    mpi::all_to_all(m_tourism.m_count_to_send, m_tourism.m_count_to_recv);

    // Заполняем m_recv_offsets
    for(int r = 1; r < size; ++r){
        m_tourism.m_recv_offsets[r] = m_tourism.m_recv_offsets[r - 1] + m_tourism.m_count_to_recv[r - 1];
    }
    // Заполняем m_border
    int border_size = m_tourism.m_send_offsets[size - 1] + m_tourism.m_count_to_send[size - 1];
    m_tourism.m_border = m_locals;
    m_tourism.m_border.resize(border_size);

    int temp_border_it = 0;
    for(auto& border_indices : m_tourism.m_border_indices){
        for(auto cell_index : border_indices){
            // [?]
            memcpy(m_tourism.m_border[temp_border_it].ptr(), m_locals[cell_index].ptr(), m_locals.itemsize());

            ++temp_border_it;
        }
    }
#endif
}

} // namespace zephyr::mesh