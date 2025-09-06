// Создание сетки с осевой симметрией. Наличие осевой симметрии можно увидеть
// по расположению центров ячеек и граней, а также объемам/длинам граней.

#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/mesh/euler/eu_mesh.h>

using zephyr::geom::Vector3d;
using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
using zephyr::geom::generator::Sector;
using zephyr::geom::generator::Rectangle;

using zephyr::io::VtuFile;
using zephyr::io::Variables;

double get_volume(EuCell& cell) { return cell.volume(); }
double get_vol_as(EuCell& cell) { return cell.volume(true); }

// Массив центров ячеек и граней
EuMesh centers(EuMesh& mesh) {
    EuMesh markers(2, false);

    // Добавляем точки в виде маркеров
    for (auto& cell: mesh) {
        double size = 0.1 * cell.linear_size();
        markers.add_marker(cell.center(), size);
        for (auto& face: cell.faces()) {
            markers.add_marker(face.center(), 0.1 * face.area());
        }
    }
    return markers;
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
    EuMesh mesh_plain(gen);

    gen.set_axial(true);
    EuMesh mesh_axial(gen);

    std::cout << "Area   (plain): " << calc_volume(mesh_plain, false) << " / " << volume_plain << "\n";
    std::cout << "Volume (axial): " << calc_volume(mesh_axial, true ) << " / " << volume_axial << "\n";

    Variables vars;
    vars += {"volume", get_volume};
    vars += {"vol_as", get_vol_as};
    EuMesh centers_pl = centers(mesh_plain);
    EuMesh centers_ax = centers(mesh_axial);

    VtuFile::save("out/mesh_plain.vtu", mesh_plain, vars, true);
    VtuFile::save("out/mesh_axial.vtu", mesh_axial, vars, true);

    VtuFile::save("out/centers_pl.vtu", centers_pl);
    VtuFile::save("out/centers_ax.vtu", centers_ax);

    return 0;
}
