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

    for (int i = 0; i < size; ++i) {
        if (i != rank && m_count_to_send[i] > 0) {
            MPI_Wait(&m_requests_send[i], MPI_STATUSES_IGNORE);
        }
    }
    for (int i = 0; i < size; ++i) {
        if (i != rank && m_count_to_recv[i] > 0) {
            MPI_Wait(&m_requests_recv[i], MPI_STATUSES_IGNORE);
        }
    }
#endif
}

void Tourism::build_border(AmrStorage& locals, AmrStorage& aliens){
#ifdef ZEPHYR_MPI
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Заполняем m_border_indices
    for (auto& cell: locals) {
        // build alien можно вызвать для не совсем нормальной сетки
        if (cell.is_undefined())
            continue;
        for (auto& face: cell.faces) {
            if(face.is_undefined()) 
                continue;
            auto& border_indices = m_border_indices[face.adjacent.rank];
            if(face.adjacent.rank != rank && (border_indices.empty() || border_indices.back() != cell.index))
                border_indices.push_back(cell.index);
        }
    }

    // Заполняем m_count_to_send
    for(int r = 0; r < size; ++r)
        m_count_to_send[r] = m_border_indices[r].size();
    // Заполняем m_send_offsets
    for(int r = 1; r < size; ++r){
        m_send_offsets[r] = m_send_offsets[r - 1] + m_count_to_send[r - 1];
    }

    // Отправляем m_count_to_send -> получаем m_count_to_recv
    mpi::all_to_all(m_count_to_send, m_count_to_recv);

    // Заполняем m_recv_offsets
    for(int r = 1; r < size; ++r){
        m_recv_offsets[r] = m_recv_offsets[r - 1] + m_count_to_recv[r - 1];
    }
    // Заполняем m_border
    int border_size = m_send_offsets[size - 1] + m_count_to_send[size - 1];
    // [!] штука ниже мне не нравится, но как сделать иначе?
    m_border = locals;
    m_border.resize(border_size);

    int temp_border_it = 0;
    for(auto& border_indices : m_border_indices){
        for(auto cell_index : border_indices){
            // [?]
            memcpy(m_border[temp_border_it].ptr(), locals[cell_index].ptr(), locals.itemsize());

            ++temp_border_it;
        }
    }

    aliens.resize(m_recv_offsets[size - 1] + m_count_to_recv[size - 1]);

	for(int i = 0; i < size; ++i){
		m_count_to_send[i] *= aliens.itemsize();
		m_count_to_recv[i] *= aliens.itemsize();
		m_send_offsets[i] *= aliens.itemsize();
		m_recv_offsets[i] *= aliens.itemsize();
	}
#endif
}

} // namespace zephyr::mesh