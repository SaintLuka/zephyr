#include <zephyr/utils/mpi.h>

namespace zephyr::utils {

#ifdef ZEPHYR_ENABLE_MPI


static int g_size = -1;
static int g_rank = -1;

void mpi::init() {
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
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

bool mpi::master() {
    return g_rank < 1;
}

bool mpi::single() {
    return g_size < 2;
}

void mpi::barrier() {
    MPI_Barrier(MPI_COMM_WORLD);
}

#else

std::ostream& mpi::cout = std::cout;

#endif

} // namespace zephyr::utils
