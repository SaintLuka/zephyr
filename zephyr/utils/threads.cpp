#include <memory>
#include <thread>
#include <cmath>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>

namespace zephyr::utils {

int threads::n_threads = 1;

#ifdef ZEPHYR_TBB
tbb::global_control threads::m_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, 1);
#endif

#ifdef ZEPHYR_STD_THREADS
std::unique_ptr<ThreadPool> threads::pool = nullptr;
#endif

int threads::recommended() {
    int HC = int(std::thread::hardware_concurrency());
    int n_tasks = zephyr::utils::mpi::n_tasks();

    // Проблема не решена, округляю в большую сторону,
    // как вариант разделить между процессами, чтобы все
    // точно складывалось. Но как тогда балансировать?
    int res = std::ceil(HC / double(n_tasks));

    return res;
}

void threads::on() {
    on(recommended());
}

void threads::on(int count) {
    n_threads = std::max(1, std::min(count, recommended()));

#ifdef ZEPHYR_TBB
    m_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, n_threads);
#endif

#ifdef ZEPHYR_OPENMP
    omp_set_num_threads(n_threads);
#endif

#ifdef ZEPHYR_STD_THREADS
    if (n_threads < 2) {
        pool = nullptr;
    } else {
        pool = std::make_unique<ThreadPool>(n_threads);
    }
#endif
}

void threads::info() {
    if (mpi::single()) {
        std::cout << "Threads count: " << threads::count() << "\n\n";
    }
    else {
        mpi::cout << "MPI processes: " << mpi::size() << "\n";
        mpi::for_each([]() {
            std::cout << "  " << "Threads count: " << threads::count() << "\n";
            std::cout.flush();
        });
        mpi::cout << "\n";
    }
}

void threads::off() {
    n_threads = 1;
#ifdef ZEPHYR_TBB
    m_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, 1);
#endif

#ifdef ZEPHYR_OPENMP
    omp_set_num_threads(1);
#endif

#ifdef ZEPHYR_STD_THREADS
    pool = nullptr;
#endif
}

bool threads::active() {
    return n_threads > 1;
}

int threads::count() {
    return n_threads;
}

} // namespace zephyr::utils