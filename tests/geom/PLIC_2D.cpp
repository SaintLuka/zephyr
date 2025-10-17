// PLIC реконструкция границы на двумерных декартовых сетках

#include <boost/type_traits/type_with_alignment.hpp>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/sections.h>

#include <zephyr/geom/geom.h>

#include "zephyr/math/funcs.h"

using namespace zephyr;
using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Rectangle;

Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

/// @brief Характеристическая функция (функция-индикатор)
using InFunction = std::function<bool(const Vector3d &)>;

/// @brief Пространственная функция
using SpFunction = std::function<double(const Vector3d &)>;

// 1. Гладкая замкнутая кривая
InFunction smooth_func = [](const Vector3d &v) -> bool {
    const double R1 = 0.3;
    const double e1 = 0.2;
    const double R2 = 0.7;
    const double e2 = 0.3;
    const double k = 5.0;

    double phi = std::atan2(v.y(), v.x());
    double r = v.norm();
    double r1 = R1 * (1.0 + e1 * std::sin(k * phi));
    double r2 = R2 * (1.0 + e2 * std::sin(k * phi));
    return r1 < r && r < r2;
};

// 2. Некоторая угловатая функция
InFunction angle_func = [](const Vector3d &v) -> bool {
    static bool initialized = false;
    static Polygon poly(10);
    if (!initialized) {
        double R = 0.8;
        double r = 0.382 * R;
        for (int i = 0; i < 10; i += 2) {
            Vector3d v1 = to_cartesian(R, M_PI*(i + 0)/5.0);
            Vector3d v2 = to_cartesian(r, M_PI*(i + 1)/5.0);
            poly.set(i, v1);
            poly.set(i + 1, v2);
        }
        initialized = true;
    }
    return poly.inside(v);
};

// 3. Функция с диффузной границей
SpFunction diffuse_func = [](const Vector3d &v) -> double {
    const double R1 = 0.3;
    const double e1 = 0.2;
    const double R2 = 0.7;
    const double e2 = 0.3;
    const double k = 5.0;

    double phi = std::atan2(v.y(), v.x());
    double r = v.norm();
    double r1 = R1 * (1.0 + e1 * std::sin(k * phi));
    double r2 = R2 * (1.0 + e2 * std::sin(k * phi));
    double h = 0.12;
    return 0.25 * (1 + math::sign_p(r - r1, h)) * (1 + math::sign_p(r2 - r, h));
};

// 4. Зашумленная функция (случайный шум)
InFunction noise_func = [](const Vector3d &v) -> bool {
    return false;
};

// 5. Функция с мелкими деталями (маленькие круги, тонкие стержни)
// Соприкасающиеся окружности.
InFunction detail_func = [](const Vector3d &v) -> bool {
    return false;
};

static Storable<double> a;
static Storable<double> a2;
// Простые производные
static Storable<Vector3d> n1;
static Storable<double> p1;
static Storable<double> e1;
// Моя формула
static Storable<Vector3d> n2;
static Storable<double> p2;
static Storable<double> e2;

inline bool mixed(EuCell& cell) {
    return 0.0 < cell[a] && cell[a] < 1.0;
}

void make_interface(EuMesh& mesh) {
    mesh.for_each([](EuCell& cell) {
        if (!mixed(cell)) {
            cell[n1] = Vector3d::Zero();
            cell[n2] = Vector3d::Zero();
            cell[p1] = 0.0;
            cell[p2] = 0.0;
            return;
        }

        Vector3d grad1 = Vector3d::Zero();
        Vector3d grad2 = Vector3d::Zero();
        for (auto& face: cell.faces()) {
            double a_f1 = 0.5 * (cell[a] + face.neib(a));
            Vector3d S = face.area() * face.normal();
            grad1 -= a_f1 * S;

            double a_f2 = face_fraction(cell[a], face.neib(a));
            grad2 -= a_f2 * S;
        }
        cell[n1] = grad1.normalized();
        cell[n2] = grad2.normalized();

        double h = cell.linear_size();
        cell[p1] = quad_find_section(cell[a], cell[n1]) * h;
        cell[p2] = quad_find_section(cell[a], cell[n2]) * h;
    });
}

EuMesh body(EuMesh& mesh, Storable<double> p, Storable<Vector3d> n) {
    EuMesh clipped(2, false);
    for (auto cell: mesh) {
        if (cell[a] <= 0.0 || (cell[a] < 0.5 && cell[n].isZero())) {
            continue;
        }
        if (cell[a] >= 1.0 || (cell[a] > 0.5 && cell[n].isZero())) {
            clipped.push_back(cell.polygon());
            continue;
        }

        Vector3d point = cell.center() + cell[p] * cell[n];
        auto poly = cell.polygon().clip(point, cell[n]);
        clipped.push_back(poly);
    }
    return clipped;
}

std::tuple<double, double> calc_errors(EuMesh& mesh, InFunction func, int nx) {
    double err1 = 0.0;
    double err2 = 0.0;

    mesh.for_each([func, nx, &err1, &err2](EuCell& cell) {
        if (!mixed(cell)) {
            cell[e1] = 0.0;
            cell[e2] = 0.0;
            return;
        }

        Box box = cell.bbox();

        double pos1 = cell[p1];
        double pos2 = cell[p2];
        Vector3d norm1 = cell[n1];
        Vector3d norm2 = cell[n2];

        int counter1 = 0;
        int counter2 = 0;
        for (int i = 0; i < nx; ++i) {
            double x = box.vmin.x() + (box.vmax.x() - box.vmin.x()) * (i + 0.5) / nx;
            for (int j = 0; j < nx; ++j) {
                double y = box.vmin.y() + (box.vmax.y() - box.vmin.y()) * (j + 0.5) / nx;
                Vector3d r = {x, y, 0.0};
                bool in0 = func(r);
                bool in1 = (r - cell.center()).dot(norm1) < pos1;
                bool in2 = (r - cell.center()).dot(norm2) < pos2;
                if (in0 != in1) { ++counter1; }
                if (in0 != in2) { ++counter2; }
            }
        }
        if (counter1 == nx * nx) {
            // Выяснить, что здесь такое вообще
            cell[e1] = 0.0;
            cell[e2] = 0.0;
        }
        else {
            cell[e1] = double(counter1) / (nx * nx);// * cell.volume();
            cell[e2] = double(counter2) / (nx * nx);// * cell.volume();

            err1 += cell[e1];
            err2 += cell[e2];
        }
    });

    int mix_count = mesh.sum([](EuCell& cell) -> int {
        return mixed(cell) ? 1 : 0;
    }, 0);
    err1 /= mix_count;
    err2 /= mix_count;

    return {err1, err2};
}

void save_mesh(EuMesh& mesh) {
    Variables vars;
    vars.append("a", a);
    vars.append("a2", a2);
    vars.append("p1", p1);
    vars.append("p2", p2);
    vars.append("n1", n1);
    vars.append("n2", n2);
    vars.append("e1", e1);
    vars.append("e2", e2);

    VtuFile::save("output/mesh.vtu", mesh, vars);

    auto body1 = body(mesh, p1, n1);
    VtuFile::save("output/body1.vtu", body1, {}, false, true);

    //auto body2 = body(mesh, p2, n2);
    //VtuFile::save("output/body2.vtu", body2, {}, false, true);
}

// Для обычной характеристической функции
void show_classic(EuMesh& mesh, InFunction func) {
    mesh.for_each([func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 10000);
    });

    make_interface(mesh);

    auto [error1, error2] = calc_errors(mesh, func, 200);

    std::cout << "Error1: " << error1 << "\n";
    std::cout << "Error2: " << error2 << "\n";

    save_mesh(mesh);
}

void show_diffuse(EuMesh& mesh) {
    mesh.for_each([](EuCell& cell) {
        cell[a] = cell.integrate_low(diffuse_func, 10) / cell.volume();
        cell[a] = math::between(cell[a], 0.0, 10);
    });

    make_interface(mesh);

    save_mesh(mesh);
}

void show_noise(EuMesh& mesh, InFunction func) {
    mesh.for_each([func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 100000);
    });

    // Добавить шум
    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            double avg_a = 0.0;
            for (int ii = i - 2; ii <= i + 2; ++ii) {
                for (int jj = j - 2; jj <= j + 2; ++jj) {
                    avg_a += mesh(ii, jj)[a];
                }
            }
            avg_a /= 25;
            mesh(i, j)[a2] = 0.5 - std::abs(avg_a - 0.5);
        }
    }

    mesh.for_each([func](EuCell& cell) {
        double r = 2.0 * rand() / double(RAND_MAX) - 1.0;
        double alp = cell[a] + 0.05 * r * cell[a2];
        cell[a] = math::between(alp, 0.0, 1.0);
    });

    make_interface(mesh);

    auto [error1, error2] = calc_errors(mesh, func, 200);

    std::cout << "Error1: " << error1 << "\n";
    std::cout << "Error2: " << error2 << "\n";

    save_mesh(mesh);
}

// В норме L_inf: ц.р. не сходятся, мой метод с 1-ым порядком
// В норме L_1: ц.р. имеют 1-ый порядок, мой метод 2-ой
// Средняя от L_inf: ц.р. имеют 1-ый порядок, мой метод 2-ой
// Во всех случаях обычные разности дают погрешности на порядок больше
int main() {
    utils::threads::on();

    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_nx(100);

    EuMesh mesh(gen);

    a = mesh.add<double>("a");
    a2 = mesh.add<double>("a2");
    p1 = mesh.add<double>("p1");
    n1 = mesh.add<Vector3d>("n1");
    e1 = mesh.add<double>("e1");
    p2 = mesh.add<double>("p2");
    n2 = mesh.add<Vector3d>("n2");
    e2 = mesh.add<double>("e2");

    int test = 4;

    switch (test) {
        case 1: show_classic(mesh, smooth_func); return 0;
        case 2: show_classic(mesh, angle_func); return 0;
        case 3: show_diffuse(mesh); return 0;
        case 4: show_noise(mesh, smooth_func); return 0;
        default: return 0;
    }
}