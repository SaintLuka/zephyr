// @brief Тяжелый тест на динамическую декомпозицию и адаптацию.
// Процессы хаотично перемещаются, ячейки хаотично адаптируются.

#include <iostream>
#include <iomanip>
#include <random>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/rwalk.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using namespace zephyr::geom;

using zephyr::utils::mpi;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;
using zephyr::utils::Stopwatch;

using zephyr::mesh::decomp::RWalk;

// Некоторая странная функция (случайный тригонометрический ряд
// для генерации случайных возмущений)
struct Membrane {
    int count;
    std::vector<double> A, n, m, q;
    std::vector<double> x0, y0, t0;

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
};

// Поле cell(u) задается из Membrane, функция
// Поле cell(wflag) задается по процентилям от cell(u)
void setup_values(EuMesh& mesh, double t,
        Storable<double> val, Storable<double> wflag) {
    static Membrane func_u(1);

    double u_min = +1.0e300;
    double u_max = -1.0e300;
    for (auto &cell: mesh) {
        double u = func_u(cell.center(), t);
        cell(val) = u;
        u_min = std::min(u, u_min);
        u_max = std::max(u, u_max);
    }

    u_min = mpi::min(u_min);
    u_max = mpi::max(u_max);

    // Нормируем от 0.0 до 1.0
    for (auto &cell: mesh) {
        cell(val) = (cell(val) - u_min) / (u_max - u_min);
    }

    double count = 100.0;

    double volume = 0.0;
    std::vector<double> vols_u(count);
    for (auto &cell: mesh) {
        // Индекс от 0 до count
        double idx_u = std::floor(cell(val) * count);

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
        int idx_u = int(std::floor(cell(val) * count));
        cell(wflag) = vols_u[idx_u];
    }
}

// Флаги выставляются по процентилям в wflag
void set_flags(EuMesh& mesh, Storable<double> wflag) {
    for (auto cell: mesh) {
        if (cell(wflag) < 0.5) {
            cell.set_flag(-1);
        } else if (cell(wflag) < 0.9) {
            cell.set_flag(0);
        } else {
            cell.set_flag(1);
        }
    }
}

int main() {
    mpi::handler init;

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_nx(50);

    // Создаем сетку
    EuMesh mesh(gen);
    mesh.set_max_level(3);

    // Добавить переменные на сетку
    auto [u, wflag] = mesh.add<double>("u", "wflag");

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank", "index", "faces2D"};
    pvd.variables.append("u", u);
    pvd.variables.append("wflag", wflag);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Сложная случайная декомпозиция
    auto decmp = RWalk::create(domain, mpi::size());

    // Добавляем декомпозицию
    mesh.set_decomposition(decmp);

    //if (mesh.check_base() < 0) {
    //    std::cout << "Bad init mesh\n";
    //    return 0;
    //}

    // Начальные данные
    for (int i = 0; i < mesh.max_level(); ++i) {
        setup_values(mesh, 0.0, u, wflag);
        set_flags(mesh, wflag);
        mesh.refine();

        //if (mesh.check_refined() < 0) {
        //    std::cout << "Bad init refinement\n";
        //    return 0;
        //}
    }
    setup_values(mesh, 0.0, u, wflag);
    mesh.sync(u, wflag);

    Stopwatch elapsed(true);
    for (int n_step = 0; n_step <= 1000; ++n_step) {
        // redistribute на каждом шаге
        mesh.balancing(mesh.locals().size());
        mesh.redistribute();

        setup_values(mesh, 0.01 * n_step, u, wflag);

        if (n_step % 10 == 0) {
            mpi::cout << "\tШаг: " << std::setw(6) << n_step << "\n";
            pvd.save(mesh, n_step);
        }

        set_flags(mesh, wflag);
        mesh.refine();

        //if (mesh.check_refined() < 0) {
        //    throw std::runtime_error("Bad refined mesh");
        //}
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}