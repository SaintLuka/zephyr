#pragma once

namespace zephyr::mesh {

/// @class Просто забавное название)
class TouristAgency {
public:
    std::vector<int> m_count_to_send;
    std::vector<int> m_count_to_recv;

    std::vector<int> m_send_offsets;
    std::vector<int> m_recv_offsets;

    std::vector<std::vector<int>> m_border_indices;

    AmrStorage m_border;

#ifdef ZEPHYR_ENABLE_MPI
    std::vector<MPI_Status> m_status;
#endif
};

} // namespace zephyr::mesh