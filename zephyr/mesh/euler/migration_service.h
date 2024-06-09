#pragma once

namespace zephyr::mesh {

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

/// @class Просто забавное название)
class MigrationService {
public:

    Router m_router;

    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

    AmrStorage m_migrants;

    std::vector<MPI_Status> m_status;
};

} // namespace zephyr::mesh