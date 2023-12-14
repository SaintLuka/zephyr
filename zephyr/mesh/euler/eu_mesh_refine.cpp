#include <chrono>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/mesh/amr/apply.h>
#include <zephyr/mesh/amr/balancing.h>
#include <zephyr/mesh/amr/balancing_fast.h>


namespace zephyr { namespace mesh {

using zephyr::utils::mpi;
using zephyr::utils::Stopwatch;


void EuMesh::init_amr() {
    m_max_level = 0;

    if (m_locals.empty()) {
        return;
    }

    int rank = mpi::rank();
    int size = mpi::size();

    auto cells_nums = mpi::all_gather((int)m_locals.size());
    std::vector<int> offset(size, 0);
    for (int r = 0; r < size - 1; ++r) {
        offset[r + 1] = offset[r] + cells_nums[r];
    }

    for (int ic = 0; ic < m_locals.size(); ++ic) {
        auto& cell = m_locals[ic];

        cell.b_idx = offset[rank] + ic;
        cell.z_idx = 0;
        cell.level = 0;
        cell.flag = 0;
    }
}

bool EuMesh::is_adaptive() const {
    return m_max_level > 0;
}

int EuMesh::max_level() const {
    return m_max_level;
}

void EuMesh::set_max_level(int max_level) {
    m_max_level = std::max(0, std::min(max_level, 15));
}

void EuMesh::set_distributor(Distributor distr) {
    distributor = std::move(distr);
}

void EuMesh::refine() {
    break_nodes();

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::is_single() && m_locals.empty()) {
        throw std::runtime_error("EuMesh::refine() error: Empty mesh");
    }

    full.resume();

    balance.resume();
    if (mpi::is_single()) {
#if FAST_BALANCING
        amr::balance_flags_fast(m_locals, m_max_level);
#else
        amr::balance_flags_slow(m_locals, m_aliens, m_max_level);
#endif
    }
    else {
        //amr::balance_flags_slow(decomposition, max_level);
    }
    balance.stop();

    apply.resume();
    if (mpi::is_single()) {
        amr::apply(m_locals, distributor);
    }
    else {
        //amr::apply(decomposition, distributor);
    }
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    std::cout << "  Balance steps elapsed: " << std::setw(10) << balance.milliseconds() << " ms\n";
    std::cout << "  Apply steps elapsed:   " << std::setw(10) << apply.milliseconds() << " ms\n";
    std::cout << "Refine steps elapsed:    " << std::setw(10) << full.milliseconds() << " ms\n";
#endif
}

} // mesh
} // zephyr