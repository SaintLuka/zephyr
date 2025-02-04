/// @brief Декомпозиция сетки, уравнение переноса

#include <iostream>
#include <iomanip>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>


using namespace zephyr::mesh;
using zephyr::utils::mpi;
using zephyr::utils::Stopwatch;
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

void set_initials(EuMesh& mesh, Box domain) {
    Vector3d vc = domain.vmin + 0.2 * domain.size();
    double D = 0.1 * domain.diameter();
    for (auto& cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }
}

void set_flags(EuMesh& mesh) {
    for (auto cell: mesh) {
        cell.set_flag(-1);

        double u1 = cell(U).u1;
        for (auto face: cell.faces()) {
            double u2 = face.neib(U).u1;

            if (std::abs(u1 - u2) > 0.01) {
                cell.set_flag(1);
                break;
            }
        }
    }
}

// Если поставить обычный распределитель, то возникают
// странные дефекты, я уж думал ошибка адаптации.
Distributor conv_distributor() {
    Distributor distr = Distributor::simple();
    distr.merge = [](Children &children, AmrStorage::Item &parent) {
        double sum = 0.0;
        for (auto &child: children) {
            sum += child(U).u1 * child.volume();
        }
        parent(U).u1 = sum / parent.volume();
    };
    return distr;
}

int main() {
    mpi::init();

    threads::off();

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank", "index", "flag", "faces"};
    pvd.variables += {"u",  get_u};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_nx(123);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Создаем сетку
    EuMesh mesh(U, gen);

    mesh.set_max_level(3);
    mesh.set_distributor(conv_distributor());

    // Добавляем декомпозицию
    mesh.set_decomposition("XY");

    // Начальные данные
    for (int i = 0; i < mesh.max_level(); ++i) {
        set_initials(mesh, domain);
        mesh.exchange();
        set_flags(mesh);
        mesh.refine();
    }
    set_initials(mesh, domain);
    mesh.exchange();

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time <= 1.0 && n_step < 2000000000) {
        // Балансировка декомпозиции
        //if (n_step % 10 == 0) {
            mesh.balancing(double(mesh.locals().size()));
            mesh.redistribute();
            mesh.exchange();
        //}

        if (curr_time >= next_write) {
            mpi::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
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
        dt = mpi::min(CFL * dt);

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

        mesh.exchange();
        set_flags(mesh);
        mesh.exchange();
        mesh.refine();
        mesh.exchange();

        n_step += 1;
        curr_time += dt;
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    mpi::finalize();
    return 0;
}