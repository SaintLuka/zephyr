#include <zephyr/utils/mpi.h>

namespace zephyr { namespace utils {

#ifdef ZEPHYR_ENABLE_MPI
static int g_size = -1;
static int g_rank = -1;
#else
static int g_size = 1;
static int g_rank = 0;
#endif

void mpi::init() {
#ifdef ZEPHYR_ENABLE_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
#else
    g_size = 1;
    g_rank = 0;
#endif
}

void mpi::finalize() {
#ifdef ZEPHYR_ENABLE_MPI
    MPI_Finalize();
#endif
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

bool mpi::is_master() {
    return g_rank < 1;
}

bool mpi::is_single() {
    return g_size < 2;
}

void mpi::barrier() {
#ifdef ZEPHYR_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

} // utils
} // zephyr
