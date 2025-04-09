#include <zephyr/mesh/euler/tourism.h>

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI

Tourism::Tourism() :
        m_requests_recv(mpi::size()),
        m_requests_send(mpi::size()),
        m_count_to_send(mpi::size()),
        m_count_to_recv(mpi::size()),
        m_send_offsets(mpi::size(), 0),
        m_recv_offsets(mpi::size(), 0),
        m_border_indices(mpi::size(), std::vector<int>())
{}

void Tourism::shrink_to_fit() {
    m_requests_recv.shrink_to_fit();
    m_requests_send.shrink_to_fit();
    m_count_to_send.shrink_to_fit();
    m_count_to_recv.shrink_to_fit();
    m_send_offsets.shrink_to_fit();
    m_recv_offsets.shrink_to_fit();
    m_border_indices.shrink_to_fit();
}

void Tourism::init_types(const AmrStorage& locals) {
    m_border = locals.same();
    m_item_mpi_type = mpi::datatype::contiguous(locals.itemsize());
}

mpi::datatype Tourism::get_mpi_type() const {
    return m_item_mpi_type;
}

void Tourism::send(const AmrStorage& locals, Post post) {
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
                m_border.item(m_send_offsets[i]).ptr(), m_count_to_send[i], m_item_mpi_type,
                i, 0, mpi::comm(), &m_requests_send[i]
            );
        }
    }
}

void Tourism::recv(AmrStorage& aliens, Post post) {
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Расширим буфер, если необходимо
    aliens.resize(m_recv_offsets[size - 1] + m_count_to_recv[size - 1]);

    for (int i = 0; i < size; ++i) {
        if (i != rank && m_count_to_recv[i] != 0) {
            MPI_Irecv(
                aliens.item(m_recv_offsets[i]).ptr(), m_count_to_recv[i], m_item_mpi_type, 
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
}

void Tourism::build_border(AmrStorage& locals){
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Очистить списки и смещения
    m_recv_offsets[0] = 0;
    m_send_offsets[0] = 0;

    for (auto &border_indices : m_border_indices)
        border_indices.clear();

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
    m_border.resize(border_size);

    int temp_border_it = 0;
    for(auto& border_indices : m_border_indices){
        for(auto cell_index : border_indices){
            // [?]
            memcpy(m_border[temp_border_it].ptr(), locals[cell_index].ptr(), locals.itemsize());

            ++temp_border_it;
        }
    }
}

void Tourism::build_aliens(AmrStorage& locals, AmrStorage& aliens) {
    const int rank = mpi::rank();

    build_border(locals);

    // Отправляем locals
    send(locals);

    for (auto &cell : locals) {
        for (auto &face : cell.faces) {
            if (face.is_undefined())
                continue;

            if (face.adjacent.rank == rank)
                face.adjacent.alien = -1;
        }
    }

    // Получам в aliens
    recv(aliens);

    int al_it = 0;
    for (auto &cell: aliens) {
        for (auto &face : cell.faces) {
            if (face.is_undefined())
                continue;

            if (face.adjacent.index != -1 && face.adjacent.rank == rank) {
                //printf("I: %d\n", face.adjacent.index);
                auto &curr_cell = locals[face.adjacent.index];
                for (auto &l_face : curr_cell.faces) {
                    if (l_face.adjacent.rank == cell.rank && l_face.adjacent.index == cell.index) {
                        l_face.adjacent.alien = al_it;
                        break;
                    }
                }
            }
        }
        ++al_it;
    }

    // [?] если m_border_indices[r].size() == m_count_to_send[r], то зачем второе вообще нужно?
}

#endif

} // namespace zephyr::mesh