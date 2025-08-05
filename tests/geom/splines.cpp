// @brief В geom/curves определено несколько классов параметрических кривых
// и сплайнов, протестируем их здесь.
#include <zephyr/geom/curves/linear_spline.h>
#include <zephyr/geom/curves/cubic_spline.h>
#include <zephyr/geom/curves/lagrange.h>

#include <zephyr/utils/matplotlib.h>


namespace plt = zephyr::utils::matplotlib;

using namespace zephyr::geom;
using namespace zephyr::geom::curves;

using Dictionary = std::map<std::string, std::string>;

Dictionary node_style = {
        {"color",      "black"},
        {"linestyle",  "none"},
        {"marker",     "o"},
        {"markersize", "4"}
};

// Обычные интерполянты: y(x)
void ordinary_interpolants() {
    // Оригинальная функция
    double L = 5.0;
    auto func = [L](double x) -> double {
        return std::sin(2.0 * M_PI * x / L);
    };
    std::vector<double> x_all = linspace(0.0, L, 1000);
    std::vector<double> y_all(x_all.size());
    for (size_t i = 0; i < x_all.size(); ++i) {
        y_all[i] = func(x_all[i]);
    }

    // Узлы сплайна
    std::vector<double> xs = {0.0, 0.12, 0.3, 0.46, 0.6, 0.82, 0.94, 1.0};
    std::vector<double> ys(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        xs[i] *= L;
        ys[i] = func(xs[i]);
    }
    ys.back() = ys.front();

    LinearSpline line1(xs, ys, SplineBound::Free, SplineBound::Crop);
    LinearSpline line2(xs, ys, SplineBound::Free, SplineBound::Periodic);
    Lagrange     poly1(xs, ys, SplineBound::Free, SplineBound::Crop);
    Lagrange     poly2(xs, ys, SplineBound::Free, SplineBound::Periodic);
    CubicSpline cubic1(xs, ys, SplineBound::Crop, SplineBound::Free);
    CubicSpline cubic2(xs, ys, SplineBound::Crop, SplineBound::Periodic);

    Dictionary orig_style = {
            {"label",     "Original"},
            {"color",     "black"},
            {"linestyle", "solid"},
            {"linewidth", "0.5"}
    };

    plt::figure_size(9.0, 5.0, 150);

    // Crop/Free на границах
    plt::subplot(2, 1, 1);
    plt::grid(true);

    plt::plot(line1.xs(2000, -0.1 * L, 1.1 * L), line1.ys(2000, -0.1 * L, 1.1 * L),
              {{"label",     "Linear [Free, Crop]"},
               {"linestyle", "dashed"},
               {"color",     "blue"}});

    plt::plot(cubic1.xs(200, -0.1 * L, 1.1 * L), cubic1.ys(200, -0.1 * L, 1.1 * L),
              {{"label",     "Cubic  [Crop, Free]"},
               {"linestyle", "solid"},
               {"color",     "orange"}});

    plt::plot(poly1.xs(200, -0.1 * L, 1.1 * L), poly1.ys(200, -0.1 * L, 1.1 * L),
              {{"label",     "L poly [Free, Crop]"},
               {"linestyle", "dashed"},
               {"color",     "red"}});

    plt::plot(x_all, y_all, orig_style);

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Периодические граничные условия
    plt::subplot(2, 1, 2);
    plt::grid(true);

    plt::plot(cubic2.xs(200, -L, 1.5 * L), cubic2.ys(200, -L, 1.5 * L),
              {{"label",     "Cubic  [Periodic]"},
               {"linestyle", "solid"},
               {"color",     "orange"}});

    plt::plot(line2.xs(2000, -L, 1.5 * L), line2.ys(2000, -L, 1.5 * L),
              {{"label",     "Linear [Periodic]"},
               {"linestyle", "dashed"},
               {"color",     "blue"}});

    plt::plot(poly2.xs(200, -L, 1.5 * L), poly2.ys(200, - L, 1.5 * L),
              {{"label",     "L poly [Periodic]"},
               {"linestyle", "dashed"},
               {"color",     "red"}});

    plt::plot(x_all, y_all, orig_style);

    plt::plot(xs, ys, node_style);
    plt::legend();

    plt::suptitle("Обычные сплайны: $y(x)$");
    plt::tight_layout();
    plt::show();
}

// Параметрические интерполянты: x(t), y(t)
void parametric_interpolants() {
    // Узлы сплайна
    std::vector<double> xs = {0.0, 4.0, 10.0, 7.0, 8.0, 1.0};
    std::vector<double> ys = {1.0, 5.0,  6.0, 3.0, 0.0, 0.0};

    PLinearSpline line1(xs, ys, SplineBound::Free, SplineBound::Crop);
    PLinearSpline line2(xs, ys, SplineBound::Free, SplineBound::Periodic);
    PLagrange     poly1(xs, ys, SplineBound::Free, SplineBound::Crop);
    PLagrange     poly2(xs, ys, SplineBound::Free, SplineBound::Periodic);
    PCubicSpline cubic1(xs, ys, SplineBound::Crop, SplineBound::Free);
    PCubicSpline cubic2(xs, ys, SplineBound::Crop, SplineBound::Periodic);

    plt::figure_size(9.0, 3.8, 150);

    // Crop/Free на границах
    plt::subplot(1, 2, 1);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(line1.xs(2000, -0.1, 1.1), line1.ys(2000, -0.1, 1.1),
              {{"label",     "Linear [Free, Crop]"},
               {"linestyle", "dashed"},
               {"color",     "blue"}});

    plt::plot(cubic1.xs(200, -0.1, 1.1), cubic1.ys(200, -0.1, 1.1),
              {{"label",     "Cubic  [Crop, Free]"},
               {"linestyle", "solid"},
               {"color",     "orange"}});

    plt::plot(poly1.xs(200, -0.1, 1.1), poly1.ys(200, -0.1, 1.1),
              {{"label",     "L poly [Free, Crop]"},
               {"linestyle", "solid"},
               {"color",     "red"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Периодические граничные условия
    plt::subplot(1, 2, 2);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(cubic2.xs(200, 0.0, 1.0), cubic2.ys(200, 0.0, 1.0),
              {{"label",     "Cubic  [Periodic]"},
               {"linestyle", "solid"},
               {"color",     "orange"}});

    plt::plot(line2.xs(2000, 0.0, 1.0), line2.ys(2000, 0.0, 1.0),
              {{"label",     "Linear [Periodic]"},
               {"linestyle", "dashed"},
               {"color",     "blue"}});

    plt::plot(poly2.xs(2000, 0.0, 1.0), poly2.ys(2000, 0.0, 1.0),
              {{"label",     "L poly [Periodic]"},
               {"linestyle", "solid"},
               {"color",     "red"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    plt::suptitle("Параметрические иннтерполянты: $x(t)$, $y(t)$");
    plt::tight_layout();
    plt::show();
}

// Параметризация кубических сплайнов
void cubic_parametrization() {
    // Узлы сплайна
    std::vector<double> xs = {0.0, 4.0, 10.0, 7.0, 8.0, 1.0};
    std::vector<double> ys = {1.0, 5.0, 6.0, 3.0, 0.0, 0.0};

    PCubicSpline cubic1(xs, ys, SplineBound::Crop, SplineBound::Periodic, Parametrization::Uniform);
    PCubicSpline cubic2(xs, ys, SplineBound::Crop, SplineBound::Periodic, Parametrization::Chord);
    PCubicSpline cubic3(xs, ys, SplineBound::Crop, SplineBound::Periodic, Parametrization::Chebyshev);

    // Число тестовых точек
    int n = 50;
    int N = 200;

    plt::figure_size(9.0, 7, 150);

    // Все виды параметризации
    plt::subplot(2, 2, 1);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(cubic1.xs(N), cubic1.ys(N),
              {{"label", "Uniform"},
               {"color", "blue"}});
    plt::plot(cubic2.xs(N), cubic2.ys(N),
              {{"label", "Chord Length"},
               {"color", "green"}});
    plt::plot(cubic3.xs(N), cubic3.ys(N),
              {{"label", "Chebyshev"},
               {"color", "orange"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Равномерная параметризация
    plt::subplot(2, 2, 2);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(cubic1.xs(N), cubic1.ys(N),
              {{"label", "Uniform"},
               {"color", "blue"}});

    plt::plot(cubic1.xs(n), cubic1.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "blue"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Параметризация по длине хорд
    plt::subplot(2, 2, 3);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(cubic2.xs(N), cubic2.ys(N),
              {{"label", "Chord Length"},
               {"color", "green"}});

    plt::plot(cubic2.xs(n), cubic2.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "green"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Параметризация по корням полинома Чебышёва
    plt::subplot(2, 2, 4);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(cubic3.xs(N), cubic3.ys(N),
              {{"label", "Chebyshev"},
               {"color", "orange"}});

    plt::plot(cubic3.xs(n), cubic3.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "orange"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    plt::suptitle("Параметризация кубических сплайнов");
    plt::tight_layout();
    plt::show();
}

// Параметризация полиномов Лагранжа
void lagrange_parametrization() {
    // Узлы сплайна
    std::vector<double> xs = {0.0, 4.0, 10.0, 7.0, 8.0, 1.0};
    std::vector<double> ys = {1.0, 5.0, 6.0, 3.0, 0.0, 0.0};

    PLagrange poly1(xs, ys, SplineBound::Crop, SplineBound::Crop, Parametrization::Uniform);
    PLagrange poly2(xs, ys, SplineBound::Crop, SplineBound::Crop, Parametrization::Chord);
    PLagrange poly3(xs, ys, SplineBound::Crop, SplineBound::Crop, Parametrization::Chebyshev);

    // Число тестовых точек
    int n = 50;
    int N = 200;

    plt::figure_size(9.0, 7, 150);

    // Все виды параметризации
    plt::subplot(2, 2, 1);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(poly1.xs(N), poly1.ys(N),
              {{"label", "Uniform"},
               {"color", "blue"}});
    plt::plot(poly2.xs(N), poly2.ys(N),
              {{"label", "Chord Length"},
               {"color", "green"}});
    plt::plot(poly3.xs(N), poly3.ys(N),
              {{"label", "Chebyshev"},
               {"color", "orange"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Равномерная параметризация
    plt::subplot(2, 2, 2);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(poly1.xs(N), poly1.ys(N),
              {{"label", "Uniform"},
               {"color", "blue"}});

    plt::plot(poly1.xs(n), poly1.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "blue"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Параметризация по длине хорд
    plt::subplot(2, 2, 3);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(poly2.xs(N), poly2.ys(N),
              {{"label", "Chord Length"},
               {"color", "green"}});

    plt::plot(poly2.xs(n), poly2.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "green"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    // Параметризация по корням полинома Чебышёва
    plt::subplot(2, 2, 4);
    plt::grid(true);
    plt::set_aspect_equal();
    plt::xlim(-2.0, 11.0);
    plt::ylim(-2.0, 7.0);

    plt::plot(poly3.xs(N), poly3.ys(N),
              {{"label", "Chebyshev"},
               {"color", "orange"}});

    plt::plot(poly3.xs(n), poly3.ys(n),
              {{"linestyle", "none"},
               {"marker",    "."},
               {"color",     "orange"}});

    plt::plot(xs, ys, node_style);
    plt::legend();

    plt::suptitle("Параметризация полиномов Лагранжа");
    plt::tight_layout();
    plt::show();
}

int main() {
    // Обычные интерполянты: y(x)
    ordinary_interpolants();

    // Параметрические интерполянты: x(t), y(t)
    parametric_interpolants();

    // Параметризация кубических сплайнов
    cubic_parametrization();

    // Параметризация полиномов Лагранжа
    lagrange_parametrization();

    return 0;
}