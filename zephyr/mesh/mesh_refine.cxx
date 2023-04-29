#include <chrono>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/mesh/mesh.h>

#include <zephyr/mesh/amr/apply.h>
#include <zephyr/mesh/amr/balancing.h>
#include <zephyr/mesh/amr/balancing_fast.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr { namespace mesh {

using zephyr::utils::mpi;
using zephyr::utils::Stopwatch;


void Mesh::init_amr() {
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
        auto cell = m_locals[ic];

        cell.geom().b_idx = offset[rank] + ic;
        cell.geom().z_idx = 0;
        cell.geom().level = 0;
        cell.geom().flag = 0;
    }
}

int Mesh::max_level() const {
    return m_max_level;
}

void Mesh::set_max_level(int max_level) {
    m_max_level = std::max(0, std::min(max_level, 15));
}

void Mesh::set_distributor(Distributor distr) {
    distributor = std::move(distr);
}

void Mesh::refine() {
    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::is_single() && m_locals.empty()) {
        throw std::runtime_error("Mesh::refine() error: Empty mesh");
    }

    full.resume();

    balance.resume();
    if (mpi::is_single()) {
#if FAST_BALANCING
        amr::balance_flags_fast(m_locals, max_level);
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
    std::cout << "  Balance steps elapsed: " << balance.seconds() << " sec\n";
    std::cout << "  Apply steps elapsed:   " << apply.seconds() << " sec\n";
    std::cout << "Refine steps elapsed: " << full.seconds() << " sec\n";
#endif
}

} // mesh
} // zephyr