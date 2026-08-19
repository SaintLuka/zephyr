#include <numeric>
#include <iomanip>

#include <zephyr/mesh/euler/router.h>

namespace zephyr::mesh {

using utils::mpi;

inline std::vector<index_t> accumulate(const std::vector<index_t> &arr) {
    std::vector<index_t> res(arr.size());
    res[0] = 0;
    for (size_t i = 1; i < arr.size(); ++i) {
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
Requests::Requests(int size) : size_(size) {
    requests_ = std::make_unique<MPI_Request[]>(size_);
    std::fill_n(requests_.get(), size_, MPI_REQUEST_NULL);
}

void Requests::wait() const {
    for (int r = 0; r < size_; ++r) {
        if (requests_[r] != MPI_REQUEST_NULL) {
            MPI_Wait(&requests_[r], MPI_STATUS_IGNORE);
        }
    }
}

void RequestsList::reserve(int size) {
    requests_.reserve(size);
}

void RequestsList::operator+=(Requests&& requests) {
    requests_.emplace_back(std::move(requests));
}

void RequestsList::wait() const {
    for (auto& req: requests_) {
        req.wait();
    }
}

Router::Router() {
    size_ = mpi::size();
    send_count_  = std::vector<index_t>(size_, 0);
    send_offset_ = std::vector<index_t>(size_, 0);
    recv_count_  = std::vector<index_t>(size_, 0);
    recv_offset_ = std::vector<index_t>(size_, 0);
}

void Router::set_send_count(const std::vector<index_t> &send_count) {
    assert(size_ == send_count.size());

    send_count_  = send_count;
    send_offset_ = accumulate(send_count);
}

void Router::set_recv_count(const std::vector<index_t> &recv_count) {
    assert(size_ == recv_count.size());

    recv_count_  = recv_count;
    recv_offset_ = accumulate(recv_count);

    send_recv_.clear();
}

void Router::fill_partial() {
    // На получение
    mpi::all_to_all(send_count_, recv_count_);

    // Посчитать смещения
    recv_offset_ = accumulate(recv_count_);
}

void Router::fill_complete() {
    // Полные обмены числами (все со всеми)
    send_recv_.resize(size_ * size_);

    // Полные обмены числами (все со всеми)
    MPI_Allgather(send_count_.data(), size_, mpi::type<index_t>(),
                  send_recv_.data(), size_, mpi::type<index_t>(),
                  mpi::comm());

    // Соберем массив recv_count
    recv_count_.resize(size_);
    for (int r = 0; r < size_; ++r) {
        recv_count_[r] = get(r, mpi::rank());
    }

    // Посчитать смещения
    recv_offset_ = accumulate(recv_count_);
}

index_t Router::get(int i, int j) const {
    return send_recv_[size_ * i + j];
}

index_t Router::operator()(int i, int j) const {
    return get(i, j);
}

index_t Router::send_buffer_size() const {
    return send_offset_.back() + send_count_.back();
}

index_t Router::recv_buffer_size() const {
    return recv_offset_.back() + recv_count_.back();
}

void Router::print() const {
    if (complete()) {
        print_complete();
    } else {
        print_partial();
    }
}

void Router::print_partial() const {
    std::cout << "Rank " << mpi::rank() << ". send count: " << send_count_ << "\n";
    std::cout << "        recv count: " << recv_count_ << "\n";
}

void Router::print_complete() const {
    int n = 7;
    std::cout << "from \\ to |";
    for (int i = 0; i < size_; ++i) {
        std::cout << std::setw(n) << i << " |";
    }
    std::cout << "\n";

    for (int i = 0; i < size_; ++i) {
        std::cout << "   " << i << "      |";
        for (int j = 0; j < size_; ++j) {
            std::cout << std::setw(n) << get(i, j) << " |";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

Requests Router::isend(const utils::Buffer& src, MpiTag tag) const {
    Requests send_req(size_);
    for (int r = 0; r < size_; ++r) {
        if (send_count_[r] > 0) {
            MPI_Isend(src.get_ptr(send_offset_[r]), send_count_[r],
                      src.dtype(), r, int(tag), utils::mpi::comm(), &send_req[r]);
        }
    }
    return send_req;
}

Requests Router::irecv(utils::Buffer& dst, MpiTag tag) const {
    Requests recv_req(size_);
    for (int r = 0; r < size_; ++r) {
        if (recv_count_[r] > 0) {
            MPI_Irecv(dst.get_ptr(recv_offset_[r]), recv_count_[r],
                      dst.dtype(), r, int(tag), utils::mpi::comm(), &recv_req[r]);
        }
    }
    return recv_req;
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh