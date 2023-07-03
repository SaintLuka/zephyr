/// @file Тестирование AMR на задачах с движением больших областей с высоким
/// уровнем адаптации. Под большими областями понимаются области, которые имеют
/// ту же меру размерности, что и сама область.

#include <zephyr/mesh/mesh.h>
#include <zephyr/mesh/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;


struct _U_ {
    int bit;
};

_U_ U;

double get_bit(Storage::Item cell) {
    return cell(U).bit;
}

inline double sqr(double x) {
    return x * x;
}

/// @brief Пятиконечная звезда
struct Star {
    std::array<Vector3d, 10> vs;
    double R;

    /// @param R Радиус описаной окружности
    explicit Star(double R) : R(R) {
        std::fill(vs.begin(), vs.end(), Vector3d::Zero());
        update(0.0);
    }

    /// @brief Точка внутри полигона?
    bool is_inside(const Vector3d& v) const {
        int counter = 0;

        int n_points = int(vs.size());
        for (int i = 0; i < n_points; ++i) {
            int j = (i + 1) % n_points;

            auto v1 = vs[i];
            auto v2 = vs[j];
            v1.y() -= v.y();
            v2.y() -= v.y();
            if (v1.y() * v2.y() < 0.0) {
                double x = v1.x() - v1.y() * (v2.x() - v1.x()) / (v2.y() - v1.y());
                if (x > v.x()) {
                    ++counter;
                }
            }
        }
        // Нечетное число пересечений
        return counter % 2 == 1;
    }

    void update(double t) {
        double xi = 2.0 * M_PI * t;
        double phi = 0.00123 + xi;

        double r2 = R;
        double r1 = 0.4 * R;
        double r = r1 + (r2 - r1) * sqr(std::sin(0.5 * xi));

        double a = 0.9;
        double xc = a * std::sqrt(2.0) * std::cos(xi) / (1.0 + sqr(std::sin(xi)));
        double yc = xc * std::sin(xi);

        for (int i = 0; i < 10; i += 2) {
            vs[i].x() = xc + R * std::cos(M_PI * i / 5.0 + phi);
            vs[i].y() = yc + R * std::sin(M_PI * i / 5.0 + phi);
            vs[i + 1].x() = xc + r * std::cos(M_PI * (i + 1) / 5.0 + phi);
            vs[i + 1].y() = yc + r * std::sin(M_PI * (i + 1) / 5.0 + phi);
        }
    }
};

int solution_step(Mesh& mesh, Star& star) {
    for (auto cell: mesh) {
        cell(U).bit = star.is_inside(cell.center());
    }

    for (auto& cell: mesh) {
        if (cell(U).bit > 0) {
            cell.set_flag(1);
        }
        else {
            cell.set_flag(-1);
        }
    }

    mesh.refine();

    return mesh.check_refined();
}

int main() {
    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};
    pvd.variables += { "bit", get_bit };

    Rectangle rect(-2.0, 2.0, -1.0, 1.0);
    rect.set_nx(20);

    Mesh mesh(U, &rect);

    mesh.set_max_level(5);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    // Фигура для адаптации
    Star star(0.5);

    // Начальная адаптация
    std::cout << "Начальная адапция\n";
    for (int lvl = 0; lvl < mesh.max_level() + 4; ++lvl) {
        std::cout << "  Уровень " << lvl << "\n";
        solution_step(mesh, star);
    }

    std::cout << "\nРасчет\n";
    for (int step = 0; step < 1000; ++step) {
        if (step % 10 == 0) {
            std::cout << "  Шаг " << step << " / 1000\n";
            pvd.save(mesh, step);
        }

        star.update(step / 1000.0);

        solution_step(mesh, star);
    }

    return 0;
}