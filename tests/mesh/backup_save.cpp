/// @brief Тест с полным сохранением/восстановлением сетки

#include <iostream>
#include <iomanip>
#include <random>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using namespace zephyr::geom;

using zephyr::utils::mpi;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Rectangle;
using zephyr::utils::Stopwatch;

int main() {
    mpi::handler init;

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(50);

    EuMesh mesh(rect);
    mesh.set_decomposition("XY");
    mesh.set_max_level(3);

    // Add various fields to mesh
    constexpr int n_norms = 4;
    constexpr int n_fracs = 5;
    auto density = mesh.add<double>("density");
    auto velocity = mesh.add<Vector3d>("velocity");
    auto normals = mesh.add<Vector3d[]>("normals", n_norms);
    auto fractions = mesh.add<double[n_fracs]>("fractions");

    auto init_cells = [&]() {
        for (auto cell: mesh) {
            double x = cell.x();
            double y = cell.y();
            cell[density] = std::sin(5*x + 8*y);
            cell[velocity] = Vector3d{
                std::cos(1.0 / (x * x + y * y + 0.1)),
                std::sin(1.0 / (x * x + y * y + 0.01)),
                std::tan(1.0 / (x * x + y * y + 0.01)),
            };
            for (int i = 0; i < n_norms; ++i) {
                cell[normals][i] = {x + i, y + i, 0.0 + i};
            }
            for (int i = 0; i < n_fracs; ++i) {
                cell[fractions][i] = i * x * y;
            }
        }
    };

    auto set_flags = [&]() {
        for (auto cell: mesh) {
            if (cell[velocity].x() > 0.6) {
                cell.set_flag(1);
            }
            else {
                cell.set_flag(-1);
            }
        }
    };

    // Initialize cells data
    for (int lvl = 0; mesh.adaptive() ? lvl <= mesh.max_level() : 0; ++lvl) {
        init_cells();
        set_flags();
        mesh.refine();
    }

    // Find better decomposition
    mesh.prebalancing(50);
    init_cells();

    // Save mesh for ParaView
    PvdFile pvd("backup_save", "output");
    pvd.variables = {"level"};
    pvd.variables.append("density", density);
    pvd.variables.append("velocity", velocity);
    for (int i = 0; i < 4; ++i) {
        pvd.variables.append<Vector3d>(
            "normals[" + std::to_string(i) + "]",
            [normals, i](EuCell& cell) -> Vector3d { return cell[normals][i]; }
        );
    }
    for (int i = 0; i < 5; ++i) {
        pvd.variables.append<double>(
            "fractions[" + std::to_string(i) + "]",
            [fractions, i](EuCell& cell) -> double { return cell[fractions][i]; }
        );
    }
    pvd.save(mesh, 0.0);

    // Save full mesh
    mesh.backup("output/save", {"density", "velocity", "normals", "fractions"});

    return 0;
}