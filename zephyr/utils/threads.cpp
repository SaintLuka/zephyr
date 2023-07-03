#include <zephyr/utils/threads.h>

namespace zephyr { namespace utils {

int threads::n_threads = 1;
std::unique_ptr<ThreadPool> threads::pool = nullptr;

void threads::on() {
    n_threads = int(std::thread::hardware_concurrency());
    pool = std::unique_ptr<ThreadPool>(new ThreadPool(n_threads));
}

void threads::on(int count) {
    if (count < 1) {
        n_threads = 1;
    }
    else {
        int HC = int(std::thread::hardware_concurrency());
        n_threads = count < HC ? count : HC;
    }
    pool = std::unique_ptr<ThreadPool>(new ThreadPool(n_threads));
}

void threads::off() {
    n_threads = 1;
    pool = nullptr;
}

bool threads::is_on() {
    return n_threads > 1;
}

bool threads::is_off() {
    return n_threads < 2;
}

int threads::count() {
    return n_threads;
}

} // utils
} // zephyr