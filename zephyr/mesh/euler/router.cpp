#include <zephyr/mesh/euler/router.h>

namespace zephyr::mesh {

using zephyr::utils::mpi;

inline std::vector<index_t> accumulate(const std::vector<index_t> &arr) {
    std::vector<index_t> res(arr.size());
    res[0] = 0;
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = res[i - 1] + arr[i - 1];
    }
    return res;
}

inline index_t sum(const std::vector<index_t> &arr) {
    return std::accumulate(arr.begin(), arr.end(), index_t{0});
}

inline std::ostream &operator<<(std::ostream &os, const std::vector<index_t> &arr) {
    os << "[";
    for (size_t i = 0; i < arr.size() - 1; ++i) {
        os << arr[i] << ", ";
    }
    if (!arr.empty()) {
        os << arr.back() << "]";
    }
    return os;
}

#ifdef ZEPHYR_MPI
void Requests::wait() const {
    for (auto r: m_requests) {
        if (r != MPI_REQUEST_NULL) {
            MPI_Wait(&r, MPI_STATUS_IGNORE);
        }
    }
}

Router::Router() {
    m_size = mpi::size();
    m_send_count  = std::vector<index_t>(m_size, 0);
    m_send_offset = std::vector<index_t>(m_size, 0);
    m_recv_count  = std::vector<index_t>(m_size, 0);
    m_recv_offset = std::vector<index_t>(m_size, 0);
}

void Router::set_send_count(const std::vector<index_t> &send_count) {
    assert(m_size == send_count.size());

    m_send_count  = send_count;
    m_send_offset = accumulate(send_count);
}

void Router::set_recv_count(const std::vector<index_t> &recv_count) {
    assert(m_size == recv_count.size());

    m_recv_count  = recv_count;
    m_recv_offset = accumulate(recv_count);

    m_send_recv.clear();
}

void Router::fill_partial() {
    // На получение
    mpi::all_to_all(m_send_count, m_recv_count);

    // Посчитать смещения
    m_recv_offset = accumulate(m_recv_count);
}

void Router::fill_complete() {
    // Полные обмены числами (все со всеми)
    m_send_recv.resize(m_size * m_size);

    // Полные обмены числами (все со всеми)
    MPI_Allgather(m_send_count.data(), m_size, mpi::type<index_t>(),
                  m_send_recv.data(), m_size, mpi::type<index_t>(),
                  mpi::comm());

    // Соберем массив recv_count
    m_recv_count.resize(m_size);
    for (int r = 0; r < m_size; ++r) {
        m_recv_count[r] = get(r, mpi::rank());
    }

    // Посчитать смещения
    m_recv_offset = accumulate(m_recv_count);
}

index_t Router::get(int i, int j) const {
    return m_send_recv[m_size * i + j];
}

index_t Router::operator()(int i, int j) const {
    return get(i, j);
}

index_t Router::send_buffer_size() const {
    return m_send_offset.back() + m_send_count.back();
}

index_t Router::recv_buffer_size() const {
    return m_recv_offset.back() + m_recv_count.back();
}

void Router::print() const {
    if (complete()) {
        print_complete();
    } else {
        print_partial();
    }
}

void Router::print_partial() const {
    std::cout << "Rank " << mpi::rank() << ". send count: " << m_send_count << "\n";
    std::cout << "         recv count: " << m_recv_count << "\n";
}

void Router::print_complete() const {
    int n = 7;
    std::cout << "from \\ to |";
    for (int i = 0; i < m_size; ++i) {
        std::cout << std::setw(n) << i << " |";
    }
    std::cout << "\n";

    for (int i = 0; i < m_size; ++i) {
        std::cout << "   " << i << "      |";
        for (int j = 0; j < m_size; ++j) {
            std::cout << std::setw(n) << get(i, j) << " |";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh