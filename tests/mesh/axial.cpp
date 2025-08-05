// @brief Создание сетки с осевой симметрией. Наличие осевой симметрии можно
// увидеть по расположению центров ячеек и граней, а также объемам/длинам граней.

#include <filesystem>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/lagrange/la_node.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/io/vtu_file.h>

using zephyr::geom::Vector3d;
using zephyr::mesh::EuMesh;
using zephyr::mesh::AmrStorage;
using zephyr::mesh::NodeStorage;
using zephyr::geom::generator::Sector;
using zephyr::geom::generator::Rectangle;

using zephyr::io::VtuFile;
using zephyr::io::Variables;

namespace fs = std::filesystem;

struct _U_ {
    double value;
};

_U_ U;

double get_volume(AmrStorage::Item& cell) { return cell.volume; }
double get_vol_as(AmrStorage::Item& cell) { return cell.volume_alt; }

// Массив центров ячеек и граней
NodeStorage centers(EuMesh& mesh) {
    size_t count = 0;
    for (auto& cell: mesh) {
        ++count;
        for (auto& face: cell.faces()) {
            ++count;
        }
    }

    NodeStorage nodes_plain(count, 0);

    count = 0;
    for (auto& cell: mesh) {
        nodes_plain[count].coords = cell.center();
        ++count;
        for (auto& face: cell.faces()) {
            nodes_plain[count].coords = face.center();
            ++count;
        }
    }
    return nodes_plain;
}

double calc_volume(EuMesh& mesh, bool axial) {
    double volume = 0.0;
    for (auto &cell: mesh) {
        volume += cell.volume(axial);
    }
    return axial ? 2.0 * M_PI * volume : volume;
}

int main() {
    double R = 2.0;

    // Два варианта теста: для прямоугольника/цилиндра и диска/сферы
#if 1
    double H = 1.5 * R;

    Rectangle gen(0.0, H, 0.0, R);
    gen.set_nx(30);

    double volume_plain = R * H;
    double volume_axial = M_PI * std::pow(R, 2) * H;
#else
    Sector gen(R, R / 3.0, M_PI, false);
    gen.set_n_phi(60);

    double volume_plain = M_PI_2 * std::pow(R, 2);
    double volume_axial = 4.0 / 3.0 * M_PI * std::pow(R, 3);
#endif

    gen.set_axial(false);
    EuMesh mesh_plain(gen, U);

    gen.set_axial(true);
    EuMesh mesh_axial(gen, U);

    std::cout << "Area   (plain): " << calc_volume(mesh_plain, false) << " / " << volume_plain << "\n";
    std::cout << "Volume (axial): " << calc_volume(mesh_axial, true ) << " / " << volume_axial << "\n";

    Variables vars;
    vars += {"volume", get_volume};
    vars += {"vol_as", get_vol_as};

    if (!fs::exists("out") || (fs::exists("out") && fs::is_directory("out"))) {
        fs::create_directory("out");
    }

    NodeStorage centers_pl = centers(mesh_plain);
    NodeStorage centers_ax = centers(mesh_axial);

    VtuFile::save("out/mesh_plain.vtu", mesh_plain, vars, true);
    VtuFile::save("out/mesh_axial.vtu", mesh_axial, vars, true);

    VtuFile::save("out/centers_pl.vtu", centers_pl);
    VtuFile::save("out/centers_ax.vtu", centers_ax);

    return 0;
}
