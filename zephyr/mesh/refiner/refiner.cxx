#include <chrono>

#include <zephyr/mesh/refiner/refiner.h>
#include <zephyr/mesh/refiner/impl/apply.h>
#include <zephyr/mesh/refiner/impl/balancing.h>
#include <zephyr/mesh/refiner/impl/balancing_fast.h>
#include <zephyr/mesh/refiner/impl/check.h>

namespace zephyr { namespace mesh {

using namespace zephyr::data;
using amrData;

void Refiner::init(Storage& cells) {
#ifdef ZEPHYR_ENABLE_MPI
    if (network::mpi::main().size() > 1) {
        std::cerr << "Warning: MPI is enabled, but single process version Refiner::init() is used\n";
    }
#endif
    if (cells.empty()) return;

    impl::test_storage(cells);

    for (size_t ic = 0; ic < cells.size(); ++ic) {
        auto cell = cells[ic];

        cell[amrData].base_id = ic;
        cell[amrData].z = 0;
        cell[amrData].level = 0;
        cell[amrData].flag = 0;
    }
}

#ifdef ZEPHYR_ENABLE_MPI

std::vector<int> gather(Network& network, int value) {
    // TODO: Вынести в Network
    std::vector<int> values(network.size());
    MPI_Allgather(&value, 1, MPI_INT, values.data(), 1, MPI_INT, network.comm());
    return values;
}

void Refiner::init(Network& network, Storage& cells) {
    if (cells.empty()) return;

    impl::test_storage(cells);

    int size = network.size();
    int rank = network.rank();

    auto cells_nums = gather(network, (int)cells.size());
    std::vector<size_t> offset(network.size(), 0);
    for (int r = 0; r < size - 1; ++r) {
        offset[r + 1] = offset[r] + cells_nums[r];
    }

    for (size_t ic = 0; ic < cells.size(); ++ic) {
        auto cell = cells[ic];

        cell[amrData].base_id = offset[rank] + ic;
        cell[amrData].z = 0;
        cell[amrData].level = 0;
        cell[amrData].flag = 0;
    }
}
#endif

void print_cell_info(Storage::iterator cell) {
    impl::print_cell_info(cell);
}

void print_cell_info(Storage& locals, Storage& aliens, size_t ic) {
    impl::print_cell_info(locals, aliens, ic);
}

void visualize_cell(Storage::iterator cell) {
    impl::visualize_cell(cell);
}

int Refiner::check_base_mesh(Storage& cells) {
    Storage aliens;
    return impl::check_base_mesh(cells, aliens, 0);
}

#ifdef ZEPHYR_ENABLE_MPI
int Refiner::check_base_mesh(Decomposition& decomposition) {
    auto& network = decomposition.network();
    int res = 0;
    // Последовательно на каждом процессе для удобства отладки
    for (unsigned int r = 0; r < network.size(); ++r) {
        if (network.rank() == r) {
            res = impl::check_base_mesh(
                    decomposition.inner_elements(),
                    decomposition.outer_elements(),
                    decomposition.network().rank()
            );
        }
        MPI_Barrier(network.comm());
    }
    return res;
}
#endif

int Refiner::check_refined_mesh(Storage& cells) {
    Storage aliens;
    return impl::check_refined_mesh(cells, aliens, 0);
}

#ifdef ZEPHYR_ENABLE_MPI
int Refiner::check_refined_mesh(Decomposition& decomposition) {
    auto& network = decomposition.network();
    int res = 0;
    // Последовательно на каждом процессе для удобства отладки
    for (unsigned int r = 0; r < network.size(); ++r) {
        if (network.rank() == r) {
            res = impl::check_refined_mesh(
                    decomposition.inner_elements(),
                    decomposition.outer_elements(),
                    decomposition.network().rank()
            );
        }
        MPI_Barrier(network.comm());
    }
    return res;
}
#endif

void Refiner::refine(
        Storage& cells,
        unsigned int max_level,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads))
{
    using zephyr::performance::timer::Stopwatch;

    if (cells.empty())
        return;

    auto dim = cells[0][element].dimension;

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    full.resume();

    balance.resume();
#if FAST_BALANCING
    impl::balance_flags_fast(cells, max_level if_multithreading(, threads));
#else
    impl::balance_flags_slow(cells, max_level if_multithreading(, threads));
#endif
    balance.stop();

    apply.resume();
    impl::apply(cells, op if_multithreading(, threads));
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    std::cout << "  Balance steps elapsed: " << balance.times().wall() << "\n";
    std::cout << "  Apply steps elapsed:   " << apply.times().wall() << "\n";
    std::cout << "Refine steps elapsed: " << full.times().wall() << "\n";
#endif
}

#ifdef ZEPHYR_ENABLE_MPI
void Refiner::refine(
        Decomposition& decomposition,
        unsigned int max_level,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads))
{
    using zephyr::performance::timer::Stopwatch;

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    full.resume();

    balance.resume();
    impl::balance_flags_slow(decomposition, max_level if_multithreading(, threads));
    balance.stop();

    apply.resume();
    impl::apply(decomposition, op if_multithreading(, threads));
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    std::cout << "  Balance steps elapsed: " << balance.times().wall() << "\n";
    std::cout << "  Apply steps elapsed:   " << apply.times().wall() << "\n";
    std::cout << "Refine steps elapsed: " << full.times().wall() << "\n";
#endif
}
#endif

} // namespace zephyr
} // namespace mesh
