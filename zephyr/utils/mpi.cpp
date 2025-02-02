#include <zephyr/utils/mpi.h>

namespace zephyr::utils {

#ifdef ZEPHYR_MPI

static int g_size  = 1;
static int g_rank  = 0;
static int g_tasks = 1;

// Собирает уникальные имена процессоров/узлов
std::vector<std::string> proc_names();

void mpi::init() {
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

    auto names = proc_names();
    std::string my_name = names[g_rank];

    g_tasks = 0;
    for (int r = 0; r < g_size; ++r) {
        if (names[r] == my_name) {
            g_tasks += 1;
        }
    }
}

void mpi::finalize() {
    MPI_Finalize();
}

int mpi::rank() {
    return g_rank;
}

std::string mpi::srank() {
    return std::to_string(g_rank);
}

int mpi::size() {
    return g_size;
}

int mpi::n_tasks() {
    return g_tasks;
}

bool mpi::master() {
    return g_rank < 1;
}

bool mpi::single() {
    return g_size < 2;
}

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

// Собирает уникальные имена процессоров/узлов
std::vector<std::string> proc_names() {
    const int max_length = MPI_MAX_PROCESSOR_NAME;

    int length = 0;
    char name[max_length];

    MPI_Get_processor_name(name, &length);

    std::vector<int> lengths = mpi::all_gather(length);
    std::vector<int> offsets(lengths.size());
    offsets[0] = 0;
    for (int i = 1; i < g_size; ++i) {
        offsets[i] = offsets[i - 1] + lengths[i - 1];
    }

    int buff_size = offsets.back() + lengths.back();
    char *all_names = new char[buff_size + 1];

    MPI_Allgatherv(name, length, MPI_CHAR, all_names,
                   lengths.data(), offsets.data(), MPI_CHAR, MPI_COMM_WORLD);

    all_names[buff_size] = '\0';

    std::vector<std::string> names(g_size);
    for (int r = 0; r < g_size; ++r) {
        names[r] = std::string(all_names + offsets[r], lengths[r]);
    }

    //mpi::cout << "Process names:\n";
    //for (int r = 0; r < g_size; ++r) {
    //    mpi::cout << "  '" << names[r] << "'\n";
    //}

    return names;
}

mpi::master_stream mpi::cout = {};

#else

std::ostream& mpi::cout = std::cout;

#endif

} // namespace zephyr::utils
