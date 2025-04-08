#include <iostream>
#include <zephyr/utils/mpi.h>

using zephyr::utils::mpi;

int main() {
    mpi::handler init;

    mpi::cout << "Hello from master process\n";

    mpi::for_each([]() {
        std::cout << "Hello from " << mpi::rank() << " / " << mpi::size() << "\n";
    });

    return 0;
}