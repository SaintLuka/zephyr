#include <iostream>
#include <zephyr/utils/mpi.h>

#include <zephyr/mesh/euler/eu_mesh.h>

using zephyr::utils::mpi;

int main() {
    mpi::init();

    mpi::cout << "Hello from master process\n";

    mpi::for_each([]() {
        std::cout << "Hello from " << mpi::rank() << " / " << mpi::size() << "\n";
    });

    mpi::finalize();
    return 0;
}