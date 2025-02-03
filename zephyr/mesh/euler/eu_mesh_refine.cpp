#include <chrono>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/mesh/amr/apply.h>
#include <zephyr/mesh/amr/balancing.h>


namespace zephyr::mesh {

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

void EuMesh::balance_flags() {
#if SCRUTINY
    static bool first_time = true;
    if (first_time) {
        int res = check_base();
        if (res < 0) {
            throw std::runtime_error("Check base failed");
        }
        first_time = false;
    }
#endif
    if (mpi::single()) {
        amr::balance_flags(m_locals, m_max_level);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::balance_flags(m_locals, m_aliens, m_max_level, *this);
    }
#endif
}

void EuMesh::apply_flags() {
    if (mpi::single()) {
        amr::apply(m_locals, distributor);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::apply(m_locals, m_aliens, distributor, *this);
    }
#endif

#if SCRUTINY
    int res = check_refined();
    if (res < 0) {
        throw std::runtime_error("Check refined failed");
    }
#endif
}

void EuMesh::refine() {
    if (!is_adaptive()) { return; }

    break_nodes();

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::single() && m_locals.empty()) {
        throw std::runtime_error("EuMesh::refine() error: Empty mesh");
    }

    full.resume();

    balance.resume();
    balance_flags();
    balance.stop();

    apply.resume();
    apply_flags();
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "  Balance steps elapsed: " << std::setw(10) << balance.milliseconds() << " ms\n";
        std::cout << "  Apply steps elapsed:   " << std::setw(10) << apply.milliseconds() << " ms\n";
        std::cout << "Refine steps elapsed:    " << std::setw(10) << full.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

} // namespace zephyr::mesh