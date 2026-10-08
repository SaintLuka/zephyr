// В проекте реализована поддержка JSON, некоторые классы можно
// инициализовать из JSON строки/файла.

#include <zephyr/mesh/mesh.h>
#include <zephyr/utils/json.h>
#include <zephyr/utils/mpi.h>
#include <zephyr/io/pvd_file.h>


using zephyr::utils::mpi;
using zephyr::utils::Json;
using zephyr::io::PvdFile;
using zephyr::mesh::Mesh;

struct U {
    double x;
};

// TODO: Написать тест нормально
int main(int argc, char** argv) {
    mpi::handler init;

    Json config = Json::load(argc, argv, "json-config.json");

    Mesh mesh(config["mesh"]);

    PvdFile pvd(config["io"]);

    pvd.save(mesh, 0.0);

    return 0;
}
