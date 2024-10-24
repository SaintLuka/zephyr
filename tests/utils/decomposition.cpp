/// @brief Декомпозиция сетки, уравнение переноса

#include <iostream>
#include <iomanip>

#include <zephyr/utils/mpi.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>


using namespace zephyr::mesh;
using zephyr::utils::mpi;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;

using zephyr::mesh::decomp::ORB;

struct _U_ {
    double u1, u2;
};

_U_ U;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  {
    return cell(U).u1;
}

double get_vx(AmrStorage::Item& cell) {
    return velocity(cell.center).x();
}

double get_vy(AmrStorage::Item& cell) {
    return velocity(cell.center).y();
}


int main() {
    mpi::init();

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables += {"u",  get_u};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_nx(256);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Создаем сетку
    EuMesh mesh(U);
    if (mpi::master()) {
        mesh = EuMesh(U, gen);
    }

    // Добавляем декомпозицию
    ORB orb(domain, "XY", mpi::size());
    mesh.add_decomposition(orb, false);

    // Распределить ячейки))
    mesh.redistribute();

    // Дальше простая схема, ничего интересного
    // Начальные данные
    Vector3d vc = domain.center();
    double D = 0.1 * domain.diameter();
    for (auto& cell: mesh) {
        //if(mpi::rank()==0)
        //    printf("i_m: %d\n", cell.geom().index);
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time <= 1.0) {
        mesh.exchange();
        
        if (curr_time >= next_write) {
            std::cout << "\tРанг: " << mpi::rank() << ";"
                      << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh.locals(), curr_time);
            next_write += 0.02;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto& cell: mesh) {
            double max_area = 0.0;
            for (auto &face: cell.faces()) {
                max_area = std::max(max_area, face.area());
            }
            double dx = cell.volume() / max_area;
            dt = std::min(dt, dx / velocity(cell.center()).norm());
        }
        dt *= CFL;


        // Расчет по схеме upwind
        for (auto& cell: mesh) {
            auto& zc = cell(U);

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                const auto& zn = face.neib(U);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc.u1 + a_m * zn.u1) * face.area();
            }

            zc.u2 = zc.u1 - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto& cell: mesh) {
            cell(U).u1 = cell(U).u2;
            cell(U).u2 = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    mpi::finalize();
    return 0;
}