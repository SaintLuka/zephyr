/// @file Тестирование WENO реконструкции

#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/mesh/mesh.h>

#include <zephyr/math/calc/weno.h>

#include <zephyr/utils/matplotlib.h>

using namespace zephyr::geom;
using generator::Strip;
using generator::Rectangle;

using namespace zephyr::mesh;

using zephyr::math::WENO5;

namespace plt = zephyr::utils::matplotlib;

struct _D1_ {
    double u;
    double du_dx;

    _D1_() = default;
};

class TestFunc {
public:

    TestFunc() = default;

    static TestFunc Step();

    static TestFunc Parabola();

    static TestFunc Smooth();

    static TestFunc Hardcore();


    std::function<double(double)> func;
    std::function<double(double)> integral;
};

void test_1D(TestFunc test) {
    _D1_ U{};

    Strip gen(0.0, 1.0, Strip::Type::UNIFORM);
    gen.set_nx(54);

    EuMesh cells(U, gen);

    // Установить начальные данные (средние величины в ячейках)
    for (int i = 0; i < cells.nx(); ++i) {
        auto cell = cells[i];
        double x1 = cell.face(Side::L).x();
        double x2 = cell.face(Side::R).x();
        cell(U).u = (test.integral(x2) - test.integral(x1)) / (x2 - x1);
    }

    // Посчитаем численную производную
    for (int i = 0; i < cells.nx(); ++i) {
        auto cell = cells[i];
        auto neib_r = cell.face(Side::R).neib();
        auto neib_l = cell.face(Side::L).neib();
        cell(U).du_dx = (neib_r(U).u - neib_l(U).u) / (neib_r.center().x() - neib_l.center().x());
    }


    // Строим все полученные реконструкции

    plt::figure_size(13, 8, 150);

    plt::xlim(0.0, 1.0);
    plt::ylim(-0.1, 1.1);

    for (int i = 0; i < cells.nx(); ++i) {
        auto cell = cells[i];
        auto neib_r = cell.face(Side::R).neib();
        auto neib_l = cell.face(Side::L).neib();

        double x1 = cell.face(Side::L).x();
        double x2 = cell.face(Side::R).x();

        // Для внутренней реконструкции
        int M = 20;
        auto xs = linspace(x1, x2, M);

        // Границы ячеек
        plt::plot({x1, x1}, {-0.1, 1.1}, {{"linestyle", "solid"}, {"color", "lightgray"}, {"linewidth", "0.5"}});

        // Точная функция
        std::vector<double> exact(xs);
        for (int k = 0; k < M; ++k) {
            exact[k] = test.func(xs[k]);
        }

        plt::plot(xs, exact, {{"linestyle", "solid"}, {"color", "blue"}, {"linewidth", "0.5"}});

        // Средние значения
        plt::plot({x1, x2}, {cell(U).u, cell(U).u}, {{"linestyle", "dashed"}, {"color", "black"}, {"linewidth", "1.0"}});

        // Наклоны
        std::vector<double> slopes(xs);
        for (int k = 0; k < M; ++k) {
            slopes[k] = cell(U).u + (xs[k] - cell.x()) * cell(U).du_dx;;
        }

        plt::plot(xs, slopes, {{"linestyle", "solid"}, {"color", "green"}, {"linewidth", "1.0"}});

        // WENO
        WENO5 weno {
            cells(i - 2, 0)(U).u,
            cells(i - 1, 0)(U).u,
            cells(i - 0, 0)(U).u,
            cells(i + 1, 0)(U).u,
            cells(i + 2, 0)(U).u
        };

        plt::plot({x1, x2}, {weno.m(), weno.p()}, {{"linestyle", "solid"}, {"color", "orange"}, {"linewidth", "1.0"}});
    }

    plt::tight_layout();
    plt::show();
}



int main() {
    //TestFunc test = TestFunc::Smooth();
    //TestFunc test = TestFunc::Parabola();
    //TestFunc test = TestFunc::Step();
    TestFunc test = TestFunc::Hardcore();

    test_1D(test);

    return 0;
}

inline double sqr(double x) {
    return x * x;
}

inline double cube(double x) {
    return x * x * x;
}

TestFunc TestFunc::Step() {
    double x1{0.2}, x2{0.4}, x3{0.6}, x4{0.8};
    TestFunc res;
    res.func = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return 1.0;
        }
        if (x < x3) {
            return 0.5;
        }
        if (x < x4) {
            return 1.0;
        }
        return 0.0;
    };
    res.integral = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return x - x1;
        }
        if (x < x3) {
            return x2 - x1 + 0.5 * (x - x2);
        }
        if (x < x4) {
            return x2 - x1 + 0.5 * (x3 - x2) + x - x3;
        }
        return x2 - x1 + 0.5 * (x3 - x2) + (x4 - x3);
    };
    return res;
}

TestFunc TestFunc::Parabola() {
    double x1{0.2}, x2{0.5}, x3{0.8};
    TestFunc res;
    res.func = [=](double x) -> double {
        if (x < x1 ) {
            return 0.0;
        }
        if (x < x2) {
            return 4.0 / sqr(x2 - x1) * (x - x1) * (x2 - x);
        }
        if (x < x3) {
            return 4.0 / sqr(x3 - x2) * (x - x2) * (x3 - x);
        }
        return 0.0;
    };
    res.integral = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return -2.0 / (3.0 * sqr(x2 - x1)) * (sqr(x - x1) * (2 * x + x1 - 3 * x2));
        }
        if (x < x3) {
            return 2.0 / 3.0 * (x2 - x1) - 2.0 / (3.0 * sqr(x3 - x2)) * (sqr(x - x2) * (2 * x + x2 - 3 * x3));;
        }
        return 2.0 / 3.0 * (x2 - x1 + x3 - x2);
    };
    return res;
}

TestFunc TestFunc::Smooth() {
    double x1{0.2}, x2{0.8};
    double xc = 0.5 * (x2 + x1);
    double w = 2.0 * M_PI / (x2 - x1);

    TestFunc res;
    res.func = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return std::pow(std::sin(w * (x - xc)), 4);
        }
        return 0.0;
    };
    res.integral = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return (12.0 * w * (x - x1) - 16.0 * std::sin(w * (x - x1)) * std::cos(w * (x - x2))
                    + 2.0 * std::sin(2.0 * w * (x - x1)) * std::cos(2.0 * w * (x - x2))) / (32.0 * w);
        }
        return (12.0 * w * (x2 - x1) - 16.0 * std::sin(w * (x2 - x1)) +
                2.0 * std::sin(2.0 * w * (x2 - x1))) / (32.0 * w);
    };
    return res;
}

TestFunc TestFunc::Hardcore() {
    double x1{0.15}, x2{0.35}, x3{0.6}, x4{0.9};
    double x0 = 0.52;
    double xc = x2 - 2.0 * x0;

    TestFunc res;
    res.func = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return 30.0 * std::pow(x - x1, 3);
        }
        if (x < x3) {
            return 1.0 + 15.0 * (x - x2) * (x + xc);
        }
        if (x < x4) {
            return 350.0 * (x - x3) * std::pow(x4 - x, 3);
        }
        return 0.0;
    };
    res.integral = [=](double x) -> double {
        if (x < x1) {
            return 0.0;
        }
        if (x < x2) {
            return 7.5 * std::pow(x - x1, 4);
        }
        double Int = 7.5 * std::pow(x2 - x1, 4);
        if (x < x3) {
            return Int + 0.5 * ((x - x2) * (2 + 10 * sqr(x) + 15 * x * xc - 5 * x * x2 - 15 * xc * x2 - 5 * sqr(x2)));
        }
        Int += 0.5 * (x2 - x3) * (-2 + 5 * sqr(x2) + 15 * xc * (x2 - x3) + 5 * x2 * x3 - 10 * sqr(x3));
        if (x < x4) {
            return Int - 17.5 * sqr(x - x3) *
                         (4 * cube(x) + 3 * (x3 - 5 * x4) * sqr(x) + 2 * (sqr(x3) - 5 * x3 * x4 + 10 * sqr(x4)) * x +
                          cube(x3) - 10 * cube(x4) + 10 * x3 * sqr(x4) - 5 * sqr(x3) * x4);
        }
        Int += 17.5 * std::pow(x4 - x3, 5);
        return Int;
    };
    return res;
}