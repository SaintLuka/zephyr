#include <memory>
#include <cmath>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>

namespace zephyr::utils {

int threads::n_threads = 1;
std::unique_ptr<ThreadPool> threads::pool = nullptr;

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

    if (n_threads < 2) {
        pool = nullptr;
    } else {
        pool = std::make_unique<ThreadPool>(n_threads);
    }
}

void threads::off() {
    n_threads = 1;
    pool = nullptr;
}

bool threads::active() {
    return n_threads > 1;
}

int threads::count() {
    return n_threads;
}

} // namespace zephyr::utils