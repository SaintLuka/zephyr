/// @brief Декомпозиция сетки, уравнение переноса

#include <iostream>
#include <iomanip>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>
#include <zephyr/mesh/decomp/VD3.h>
#include <zephyr/mesh/decomp/rwalk.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using zephyr::utils::mpi;
using zephyr::utils::Stopwatch;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;

using zephyr::mesh::decomp::ORB;
using zephyr::mesh::decomp::VD3;
using zephyr::mesh::decomp::RWalk;

// Решаем два волновых уравнения по простой схеме
struct _U_ {
    double u;      ///< Некоторая гладкая функция
    double wflag;  ///< Равномерно распределенная величина [0, 1]
};

_U_ U;

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)     { return cell(U).u; }
double get_wflag(AmrStorage::Item& cell) { return cell(U).wflag; }

class Membrane {
public:
    explicit Membrane(int seed = 0) {
        count = 125;
        A.resize(count);

        n.resize(count);
        m.resize(count);
        q.resize(count);

        x0.resize(count);
        y0.resize(count);
        t0.resize(count);

        const int N = 15; // Влияет на величину формаций

        std::mt19937_64 gen(seed);
        std::uniform_int_distribution<int> uni_i(0, N);
        std::uniform_real_distribution<double> uni_d(0.0, 1.0);

        for (int i = 0; i < count; ++i) {

            x0[i] = uni_d(gen);
            y0[i] = uni_d(gen);
            t0[i] = uni_d(gen);

            n[i] = M_PI * uni_i(gen);
            m[i] = M_PI * uni_i(gen);
            q[i] = std::sqrt(n[i] * n[i] + m[i] * m[i]);

            A[i] = (2.0 * uni_d(gen) - 1.0) / std::pow((n[i] + 3) * (m[i] + 3), 0.25);
        }
    }

    double operator()(const Vector3d& v, double t) const {
        double x = v.x();
        double y = v.y();

        double res = 0;
        for (int i = 0; i < count; ++i) {
            res += A[i] *
                   std::cos(n[i] * (x - x0[i])) *
                   std::cos(m[i] * (y - y0[i])) *
                   std::cos(q[i] * (t - t0[i]));
        }
        return res;
    }

protected:
    int count;
    std::vector<double> A, n, m, q;
    std::vector<double> x0, y0, t0;
};

void setup_values(EuMesh& mesh, double t = 0.0) {
    static Membrane func_u(1);

    double u_min = +1.0e300;
    double u_max = -1.0e300;
    for (auto &cell: mesh) {
        double u = func_u(cell.center(), t);
        cell(U).u = u;
        u_min = std::min(u, u_min);
        u_max = std::max(u, u_max);
    }

    u_min = mpi::min(u_min);
    u_max = mpi::max(u_max);

    // Нормируем от 0.0 до 1.0
    for (auto &cell: mesh) {
        cell(U).u = (cell(U).u - u_min) / (u_max - u_min);
    }

    double count = 100.0;

    double volume = 0.0;
    std::vector<double> vols_u(count);
    for (auto &cell: mesh) {
        // Индекс от 0 до count
        double idx_u = std::floor(cell(U).u * count);

        double V = cell.volume();
        vols_u[idx_u] += V;
        volume += V;
    }

    volume = mpi::sum(volume);
    vols_u = mpi::sum(vols_u);

    // Кумулятивная сумма
    vols_u[0] /= volume;
    for (int i = 1; i < count; ++i) {
        vols_u[i] += vols_u[i - 1] / volume;
    }

    for (auto &cell: mesh) {
        int idx_u = int(std::floor(cell(U).u * count));
        cell(U).wflag = vols_u[idx_u];
    }
}

void set_flags(EuMesh& mesh) {
    for (auto cell: mesh) {
        if (cell(U).wflag < 0.5) {
            cell.set_flag(-1);
        } else if (cell(U).wflag < 0.9) {
            cell.set_flag(0);
        } else {
            cell.set_flag(1);
        }
    }
}


int main() {
    mpi::init();

    threads::off();

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank", "index", "level", "faces"};
    pvd.variables += {"u", get_u};
    pvd.variables += {"wflag", get_wflag};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_nx(50);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Создаем сетку
    EuMesh mesh(U, gen);

    mesh.set_max_level(1);
    mesh.set_distributor(Distributor::simple());

    //auto decmp = ORB::create(domain, "XY", mpi::size());
    auto decmp = RWalk::create(domain, mpi::size());

    // Добавляем декомпозицию
    mesh.set_decomposition(decmp);

    // Начальные данные
    for (int i = 0; i < mesh.max_level(); ++i) {
        setup_values(mesh);
        mesh.exchange();
        set_flags(mesh);
        mesh.refine();
    }
    setup_values(mesh);
    mesh.exchange();

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time <= 1.0 && n_step < 20000000) {
        // Балансировка декомпозиции
        //if (n_step % 10 == 0) {
            mesh.balancing(mesh.n_cells());
            mesh.redistribute();
            mesh.exchange();
        //}

        if (curr_time >= next_write) {
            mpi::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += 0.0;
        }

        setup_values(mesh, curr_time);

        mesh.exchange();
        set_flags(mesh);
        mesh.exchange();
        mesh.refine();
        mesh.exchange();

        n_step += 1;
        curr_time += 0.005;
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    mpi::finalize();
    return 0;
}