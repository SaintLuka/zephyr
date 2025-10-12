// Тестирование класса Polygon. Проверка функций интегрирования по полигону,
// отсечения прямой и кругом.

#include <iostream>

#include <zephyr/math/funcs.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/geom/sections.h>
#include <zephyr/utils/numpy.h>
#include <zephyr/utils/pyplot.h>

using namespace zephyr;
using namespace zephyr::geom;

// Пересечение грани между парой ячеек
struct Inter {
    double inter; // Доля пересечения
    double cos;   // Косинус нормали к грани
};

// Найти сечение двух ячеек в заданных объемных долях в лоб с помощью
// метода дихотомии.
Inter face_fraction_direct(Polygon &cell, Polygon& neib, double a1, double a2) {
    auto func = [](double a1, double a2, double cos) -> double {
        Vector3d n = {cos, std::sqrt(1.0 - cos * cos), 0.0};
        double p1 = quad_find_section(n, a1);
        double p2 = quad_find_section(n, a2);
        return p2 + cos - p1;
    };

    double cos_min = -1.0;
    double cos_max = +1.0;

    double f_min = func(a1, a2, cos_min);
    double f_max = func(a1, a2, cos_max);

    assert(f_min * f_max <= 0.0);

    if (std::abs(f_min) < 1.0e-8) {
        cos_max = cos_min;
    }
    if (std::abs(f_max) < 1.0e-8) {
        cos_min = cos_max;
    }

    while (cos_max - cos_min > 1.0e-8) {
        double cos = 0.5 * (cos_min + cos_max);
        double f_avg = func(a1, a2, cos);

        if (f_min * f_avg <= 0.0) {
            cos_max = cos;
            f_max = f_avg;
        } else {
            cos_min = cos;
            f_min = f_avg;
        }
    }

    double cos = 0.5 * math::sign(a1 - a2) * std::abs(cos_min + cos_max);

    auto[a_min, a_max] = math::sorted(a1, a2);

    double inter;
    if (a_min == 0.0 && a_max == 1.0) {
        inter = 0.5;
    }
    else if (a_min == 0.0) {
        inter = 0.0;
    }
    else if (a_max == 1.0) {
        inter = 1.0;
    }
    else {
        Vector3d n = {cos, std::sqrt(1 - cos * cos), 0.0};

        Vector3d p1 = cell.find_section(n, a1);
        Vector3d p2 = neib.find_section(n, a2);

        inter = math::between(p1.y() + (1 - p1.x()) * (p2.y() - p1.y()) / (p2.x() - p1.x()), 0.0, 1.0);
    }

    return {inter, cos};
}

// Интерфейс через пару ячеек
void test_interface(double a1, double a2) {
    // Тестовый многоугольник
    Polygon cell = {
            Vector3d{0.0, 0.0, 0.0},
            Vector3d{1.0, 0.0, 0.0},
            Vector3d{1.0, 1.0, 0.0},
            Vector3d{0.0, 1.0, 0.0}
    };
    cell.sort();
    Polygon neib = {
            Vector3d{1.0, 0.0, 0.0},
            Vector3d{2.0, 0.0, 0.0},
            Vector3d{2.0, 1.0, 0.0},
            Vector3d{1.0, 1.0, 0.0}
    };
    neib.sort();

    // Визуализация для пары параметров (a1, a2)
    if (true) {
        auto [inter1, cos1] = face_fraction_direct(cell, neib, a1, a2);

        Vector3d n = {cos1, std::sqrt(1 - cos1 * cos1), 0.0};

        Vector3d p1 = cell.find_section(n, a1);
        Vector3d p2 = neib.find_section(n, a2);

        auto clip1 = cell.clip(p1, n);
        auto clip2 = neib.clip(p2, n);

        double cos2 = face_fraction_cos(a1, a2);
        double inter2 = face_fraction(a1, a2);

        std::cout << "  Line cos:     " << cos1 << " " << cos2 << ";\terror: " << std::abs(cos1 - cos2) << "\n";
        std::cout << "  Intersection: " << inter1 << " " << inter2 << ";\terror: " << std::abs(inter1 - inter2) << "\n";

        // Строим многоугольник и сечение
        utils::pyplot plt;

        plt.figure({.figsize={10.0, 5.0}});
        plt.title("Интерфейс через пару ячеек");
        plt.set_aspect_equal();

        plt.plot(cell.xs(), cell.ys(), {.color="black"});
        plt.plot(neib.xs(), neib.ys(), {.color="black"});

        plt.marker(1.0, inter1, {.color="black",  .marker="o"});
        plt.marker(1.0, inter2, {.color="yellow", .marker="."});

        plt.fill(clip1.xs(), clip1.ys(), {.color="#0000ff3f"});
        plt.fill(clip2.xs(), clip2.ys(), {.color="#00ff003f"});

        plt.tight_layout();
        plt.show();
    }

    auto[V1, V2] = np::meshgrid(
            np::linspace(0.0, 1.0, 100),
            np::linspace(0.0, 1.0, 100)
    );

    auto asig1 = np::zeros_like(V1);
    auto asig2 = np::zeros_like(V1);
    auto error = np::zeros_like(V1);
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            asig1[i][j] = face_fraction_direct(cell, neib, V1[i][j], V2[i][j]).inter;
            asig2[i][j] = face_fraction(V1[i][j], V2[i][j]);

            //asig1[i][j] = face_fraction_direct(cell, neib, V1[i][j], V2[i][j]).cos;
            //asig2[i][j] = face_fraction_cos(V1[i][j], V2[i][j]);

            error[i][j] = std::abs(asig1[i][j] - asig2[i][j]);
        }
    }

    utils::pyplot plt;

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(V1, V2, asig1, {.cmap="jet"});
    plt.title("Итерационное вычисление");
    plt.tight_layout();

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(V1, V2, asig2, {.cmap="jet"});
    plt.title("Точная формула");
    plt.tight_layout();

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(V1, V2, error, {.cmap="jet"});
    plt.title("Погрешность");
    plt.tight_layout();

    plt.show();

}

// Попытки аппроксимации потока? C -- локальное число Куранта?
void test_flux(double C) {
    using namespace zephyr::geom;

    // Тестовый многоугольник
    Polygon cell = {
            Vector3d{0.0, 0.0, 0.0},
            Vector3d{1.0, 0.0, 0.0},
            Vector3d{1.0, 1.0, 0.0},
            Vector3d{0.0, 1.0, 0.0}
    };
    cell.sort();

    Polygon part = {
            Vector3d{1.0 - C, 0.0, 0.0},
            Vector3d{1.0, 0.0, 0.0},
            Vector3d{1.0, 1.0, 0.0},
            Vector3d{1.0 - C, 1.0, 0.0}
    };
    part.sort();

    // Просто визуализация того, что мы здесь будем считать
    if (true) {
        double cosn = std::cos(0.1 * M_PI);
        double alpha = 0.734;

        Vector3d n = {cosn, std::sqrt(1.0 - cosn * cosn), 0.0};
        Vector3d p = cell.find_section(n, alpha);

        double res1 = part.clip_area(p, n);
        double res2 = C * average_flux(alpha, cosn, C);

        // чисто для визуализации
        auto clip1 = cell.clip(p, n);
        auto clip2 = part.clip(p, n);

        // Строим многоугольник и сечение
        utils::pyplot plt;

        plt.figure({.figsize={6.0, 6.0}, .dpi=150});
        plt.title("Поток через правую грань");
        plt.set_aspect_equal();
        plt.text(0.1, 0.8, "Объемная доля: " + std::to_string(alpha));
        plt.text(0.1, 0.65, "Поток: " + std::to_string(res1));
        plt.text(0.1, 0.50, "Поток: " + std::to_string(res2));
        plt.plot(cell.xs(), cell.ys(), {.color="black"});
        plt.arrow(p.x(), p.y(), 0.1 * n.x(), 0.1 * n.y(), {.face_color="black", .head_length=0.05, .head_width=0.03});
        plt.plot(part.xs(), part.ys(), {.linestyle="dashed", .color="black"});
        plt.fill(clip1.xs(), clip1.ys(), {.color="#0000ff3f"});
        plt.fill(clip2.xs(), clip2.ys(), {.color="#00ff003f"});
        plt.tight_layout();
        plt.show();
    }

    int nc = 200;
    int nv = 100;
    auto [cosn, vols] = np::meshgrid(
            np::linspace(-1.0, 1.0, nc),
            np::linspace(1.0e-5, 1.0 - 1.0e-5, nv));

    auto flux1 = np::zeros(nc, nv);
    auto flux2 = np::zeros(nc, nv);
    auto error = np::zeros(nc, nv);

    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nv; ++j) {
            double cos = cosn[i][j];
            double vol = vols[i][j];

            Vector3d n = {cos, std::sqrt(1.0 - cos * cos), 0.0};
            Vector3d p = cell.find_section(n, vol);

            flux1[i][j] = part.clip_area(p, n);

            flux2[i][j] = C * average_flux(vol, cos, C);

            error[i][j] = flux2[i][j] - flux1[i][j];
        }
    }

    utils::pyplot plt;

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(cosn, vols, flux1, {.cmap="jet"});
    plt.tight_layout();

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(cosn, vols, flux2, {.cmap="jet"});
    plt.tight_layout();

    plt.figure({.figsize={6.0, 6.0}});
    plt.plot_surface(cosn, vols, error, {.cmap="jet"});
    plt.tight_layout();

    plt.show();
}

int main() {

    test_interface(0.23, 0.45);

    test_flux(0.3);

    return 0;
}