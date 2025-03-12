#pragma once
#include <vector>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

using zephyr::utils::mpi;

/// @brief Матрица пересылок
struct Router {

    int& operator()(int i, int j) {
        return m_data1[0];
    };

    void sync() {

    }

    // Какой лучше? Двумерный массив из boost не хочу
    std::vector<int> m_data1;
    std::vector<std::vector<int>> m_data2;

};


class MigrationService {
public:

    Router m_router;

    std::vector<int> m_send_counts;
    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_counts;
    std::vector<int> m_recv_offsets;

    AmrStorage m_migrants;

    //
    std::vector<int> m_i_sum;
    std::vector<int> m_i;
    std::vector<int> m_sum;
    std::vector<int> m;

    MigrationService() : 
        m(mpi::size(), 0),
        m_sum(mpi::size(), 0),
        m_i(mpi::size(), 0),
        m_i_sum(mpi::size(), 0),
        m_send_counts(mpi::size(), 0),
        m_send_offsets(mpi::size(), 0),
        m_recv_counts(mpi::size(), 0),
        m_recv_offsets(mpi::size(), 0)
    {}

    // Должно вызываться в начале migrate
    void reset() {
        std::fill(m_i.begin(), m_i.end(), 0);
        std::fill(m_sum.begin(), m_sum.end(), 0);
        m_i_sum[0] = 0;

        m_recv_offsets[0] = 0;
        m_send_offsets[0] = 0;
    }

#ifdef ZEPHYR_MPI
    std::vector<MPI_Status> m_status;
#endif
};

} // namespace zephyr::mesh