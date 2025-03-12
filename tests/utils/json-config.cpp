/// @brief В проекте реализована поддержка JSON, некоторые классы
/// могут быть инициализированы из JSON строки/файла.

#include <zephyr/utils/json.h>
#include <zephyr/utils/mpi.h>
#include <zephyr/io/pvd_file.h>

using zephyr::utils::mpi;
using zephyr::utils::Json;
using zephyr::io::PvdFile;
using zephyr::mesh::EuMesh;

struct U {
    double x;
};

// TODO: Написать тест нормально
int main(int argc, char** argv) {
    mpi::init();

    Json config = Json::load(argc, argv, "json-config.json");

    EuMesh mesh(config["mesh"], 7);

    PvdFile pvd(config["io"]);

    pvd.save(mesh, 0.0);

    mpi::finalize();
    return 0;
}
