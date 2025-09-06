#include <memory>
#include <thread>
#include <cmath>
#include <boost/program_options.hpp>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>

namespace po = boost::program_options;

namespace zephyr::utils {

int threads::n_threads = 1;

#ifdef ZEPHYR_TBB
std::unique_ptr<tbb::global_control> threads::m_control = nullptr;
#else
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
    if (n_threads < 2) {
        m_control = nullptr;
    } else {
        m_control = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, n_threads);
    }
#else
    if (n_threads < 2) {
        pool = nullptr;
    } else {
        pool = std::make_unique<ThreadPool>(n_threads);
    }
#endif
}

namespace {
po::variables_map parse_options(int argc, char **argv) {
    po::options_description description;

    description.add_options()("threads,t", po::value<std::string>()->default_value("on"));

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).allow_unregistered().run(), vm);
    po::notify(vm);

    return vm;
}
}

void threads::init(int argc, char** argv) {
    auto opts = parse_options(argc, argv);
    std::string arg = opts["threads"].as<std::string>();

    if (arg == "on") {
        threads::on();
        return;
    }
    if (arg == "off") {
        threads::off();
        return;
    }

    int count;
    try {
        count = std::stoi(arg);
    }
    catch (std::exception &e) {
        throw std::runtime_error("Wrong command line argument threads=" + arg);
    }

    if (count < 0) {
        std::cerr << "Incorrect command line argument threads=" << count << "; threads::on();\n";
    }
    else {
        threads::on(count);
    }
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
    m_control = nullptr;
#else
    pool = nullptr;
#endif
}


} // namespace zephyr::utils