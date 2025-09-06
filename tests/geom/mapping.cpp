// Проверка геометрии линейных отображений (Line, Quad, Cube) и их
// квадратичных аналогов (SqLine, SqQuad, SqCube), сравнение результатов
// вычислений с результатами численного интегрирования

#include <iostream>
#include <iomanip>
#include <functional>

#include <zephyr/geom/geom.h>
#include <zephyr/utils/numpy.h>

using namespace zephyr;
using namespace zephyr::geom;

// Одномерный интеграл
template <class V>
V simpson(const std::vector<double>& t, const std::function<V(double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(t.size()) - 1; ++i) {
        res += (t[i + 1] - t[i]) * (f(t[i]) + 4.0 * f(0.5 * (t[i] + t[i + 1])) + f(t[i + 1]));
    }
    return res / 6.0;
}

// Двумерный интеграл
template <class V>
V simpson(
        const std::vector<double>& xi,
        const std::vector<double>& eta,
        const std::function<V(double, double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(xi.size()) - 1; ++i) {
        double xi1 = xi[i];
        double xi2 = xi[i + 1];
        double xic = 0.5*(xi1 + xi2);
        for (int j = 0; j < int(eta.size()) - 1; ++j) {
            double eta1 = eta[j];
            double eta2 = eta[j + 1];
            double etac = 0.5 * (eta1 + eta2);

            res += (xi2 - xi1) * (eta2 - eta1) * (
                    f(xi1, eta1) + 4 * f(xi1, etac) + f(xi1, eta2) +
                    4 * f(xic, eta1) + 16 * f(xic, etac) + 4 * f(xic, eta2) +
                    f(xi2, eta1) + 4 * f(xi2, etac) + f(xi2, eta2)
            );
        }
    }
    return res / 36.0;
}

// Трехмерный интеграл
template <class V>
V simpson(
        const std::vector<double>& xi,
        const std::vector<double>& eta,
        const std::vector<double>& chi,
        const std::function<V(double, double, double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(xi.size()) - 1; ++i) {
        double xi1 = xi[i];
        double xi2 = xi[i + 1];
        double xic = 0.5 * (xi1 + xi2);
        for (int j = 0; j < int(eta.size()) - 1; ++j) {
            double eta1 = eta[j];
            double eta2 = eta[j + 1];
            double etac = 0.5 * (eta1 + eta2);

            for (int k = 0; k < int(chi.size()) - 1; ++k) {
                double chi1 = chi[k];
                double chi2 = chi[k + 1];
                double chic = 0.5 * (chi1 + chi2);


                res += (xi2 - xi1) * (eta2 - eta1) * (chi2 - chi1) * (
                        64.0 * f(xic, etac, chic) +
                        16.0 * (
                                f(xic, etac, chi1) + f(xic, etac, chi2) +
                                f(xic, eta1, chic) + f(xic, eta2, chic) +
                                f(xi1, etac, chic) + f(xi2, etac, chic)
                        ) +
                        4.0 * (
                                f(xic, eta1, chi1) + f(xic, eta1, chi2) +
                                f(xic, eta2, chi1) + f(xic, eta2, chi2) +
                                f(xi1, etac, chi1) + f(xi1, etac, chi2) +
                                f(xi2, etac, chi1) + f(xi2, etac, chi2) +
                                f(xi1, eta1, chic) + f(xi1, eta2, chic) +
                                f(xi2, eta1, chic) + f(xi2, eta2, chic)
                        ) +
                        f(xi1, eta1, chi1) + f(xi1, eta1, chi2) +
                        f(xi1, eta2, chi1) + f(xi1, eta2, chi2) +
                        f(xi2, eta1, chi1) + f(xi2, eta1, chi2) +
                        f(xi2, eta2, chi1) + f(xi2, eta2, chi2)

                );
            }
        }
    }
    return res / (6.0 * 6.0 * 6.0);
}

int main() {
    // Основной тестовый двумерный набор
    SqQuad test_sqquad = {
            Vector3d(0.1, 1.0, 0.0), Vector3d(0.7, 1.1, 0.0), Vector3d(1.3, 1.1, 0.0),
            Vector3d(0.2, 1.5, 0.0), Vector3d(0.8, 1.7, 0.3), Vector3d(1.4, 1.8, 0.0),
            Vector3d(0.0, 2.0, 0.0), Vector3d(0.5, 2.4, 0.0), Vector3d(1.2, 2.6, 0.0)
    };

    // Основной тестовый трехмерный набор
    Cube test_cube = {
            Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.2, 0.1),
            Vector3d(0.0, 1.0, 0.3), Vector3d(1.0, 1.2, 0.2),
            Vector3d(0.0, 0.1, 1.0), Vector3d(1.5, 0.5, 1.0),
            Vector3d(0.4, 1.0, 1.0), Vector3d(1.6, 1.5, 1.3),
    };

    // Характерная длина, площадь, объем
    double L1 = (test_sqquad.vs<+1, +1>() - test_sqquad.vs<-1, -1>()).norm();
    double L2 = L1 * L1;
    double L3 = L1 * L1 * L1;

    // Центр ячейки
    Vector3d C2 = test_sqquad.center();
    Vector3d C3 = test_cube.center();

    std::cout << std::setprecision(6);

    {
        std::cout << "Длина и нормаль простого отрезка\n";

        Line line = {test_sqquad[0], test_sqquad[2]};

        // int normal(x) * dl(x) = int normal(x) |J(x)| dx = length() * normal()
        Vector3d res1 = line.length() * line.normal(C2);

        std::function<Vector3d(double)> func =
                [&line, &C2](double x) -> Vector3d {
                    return line.Jacobian() * line.normal(C2);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res2 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Длина и нормаль криволинейного отрезка\n";

        SqLine line = {test_sqquad[0], test_sqquad[1], test_sqquad[2]};

        // int normal(x) * dl(x) = int normal(x) |J(x)| dx = length() * normal()
        Vector3d res1 = line.length() * line.normal(C2);

        std::function<Vector3d(double)> func =
                [&line, &C2](double x) -> Vector3d {
                    return line.Jacobian(x) * line.normal(x, C2);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Площадь простой двумерной ячейки\n";

        Quad quad = test_sqquad.reduce();

        double res1 = quad.area();

        std::function<double(double, double)> func =
                [&quad](double x, double y) -> double {
                    return quad.Jacobian(x, y);
                };

        double res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                np::linspace(-1.0, 1.0, 51),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Площадь и нормаль простой изогнутой двумерной ячейки\n";

        Quad quad = test_sqquad.reduce();
        quad(+1, -1).z() += 0.1;
        quad(+1, +1).z() += 0.2;

        Vector3d res1 = quad.area() * quad.normal(C3);

        std::function<Vector3d(double, double)> func =
                [&quad, C3](double x, double y) -> Vector3d {
                    return quad.Jacobian(x, y) * quad.normal(x, y, C3);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                np::linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L2 << "\n\n";
    }

    {
        std::cout << "Площадь криволинейной двумерной ячейки\n";

        SqQuad& quad = test_sqquad;
        for (int i = 0; i < 9; ++i) {
            quad[i].z() = 0.0;
        }

        double res1 = quad.area();

        std::function<double(double, double)> func =
                [&quad](double x, double y) -> double {
                    return quad.Jacobian(x, y);
                };

        double res2 = simpson(
                np::linspace(-1.0, 1.0, 101),
                np::linspace(-1.0, 1.0, 101),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем обычной двумерной ячейки в осесимметричной постновке\n";

        Quad quad = test_sqquad.reduce();

        double res1 = quad.volume_as();

        std::function<double(double, double)> func =
                [&quad](double x, double y) -> double {
                    return quad(x, y).y() * quad.Jacobian(x, y);
                };

        double res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                np::linspace(-1.0, 1.0, 51),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем криволинейной ячейки в осесимметричной постновке\n";

        SqQuad& quad = test_sqquad;
        for (int i = 0; i < 9; ++i) {
            quad[i].z() = 0.0;
        }

        double res1 = quad.volume_as();

        std::function<double(double, double)> func =
                [&quad](double x, double y) -> double {
                    return quad(x, y).y() * quad.Jacobian(x, y);
                };

        double res2 = simpson(
                np::linspace(-1.0, 1.0, 501),
                np::linspace(-1.0, 1.0, 501),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем трехмерной ячейки\n";

        Cube& cube = test_cube;

        double res1 = cube.volume();

        std::function<double(double, double, double)> func =
                [&cube](double x, double y, double z) -> double {
                    return cube.Jacobian(x, y, z);
                };

        double res2 = simpson(
                np::linspace(-1.0, 1.0, 41),
                np::linspace(-1.0, 1.0, 41),
                np::linspace(-1.0, 1.0, 41),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L3 << "\n\n";
    }

    {
        std::cout << "Барицентр простой двумерной ячейки\n";

        Quad quad = test_sqquad.reduce();

        double S = quad.area();

        Vector3d res1 = quad.centroid(S);

        std::function<Vector3d(double, double)> func =
                [&quad](double x, double y) -> Vector3d {
                    return quad(x, y) * quad.Jacobian(x, y);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 51),
                np::linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0}) / S;

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Барицентр криволинейной двумерной ячейки\n";

        SqQuad& quad = test_sqquad;

        double S = quad.area();

        Vector3d res1 = quad.centroid(S);

        std::function<Vector3d(double, double)> func =
                [&quad](double x, double y) -> Vector3d {
                    return quad(x, y) * quad.Jacobian(x, y);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 501),
                np::linspace(-1.0, 1.0, 501),
                func, {0.0, 0.0, 0.0}) / S;

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Барицентр трехмерной ячейки\n";

        Cube& cube = test_cube;

        double V = cube.volume();
        Vector3d res1 = cube.centroid(V);

        std::function<Vector3d(double, double, double)> func =
                [&cube](double x, double y, double z) -> Vector3d {
                    return cube(x, y, z) * cube.Jacobian(x, y, z);
                };

        Vector3d res2 = simpson(
                np::linspace(-1.0, 1.0, 11),
                np::linspace(-1.0, 1.0, 11),
                np::linspace(-1.0, 1.0, 11),
                func, {0.0, 0.0, 0.0}) / V;

        std::cout << "    res (formula): " << res1.transpose() << "\n";
        std::cout << "    res (numeric): " << res2.transpose() << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }
}
