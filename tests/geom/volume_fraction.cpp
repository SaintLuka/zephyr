#include <iostream>
#include <iomanip>

#include <zephyr/geom/primitives/cube.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/geom/primitives/polygon.h>

using namespace zephyr::geom;

template <typename T>
inline T sqr(const T& x) {
    return x * x;
}

// Гладкая функция
double test_func1(const Vector3d& v) {
    return std::sin(10.0 * v.y() * v.y()) + std::sin(10.0 * v.x() * v.x()) + std::cos(5.0 * v.x() + 5 * v.y());
}

// Интеграл test_func1 по квадрату [0, 1]^2.
double test_integral1() {
    return 0.498542242881323557872;
}

inline double acc_ord(int n1, double err1, int n2, double err2) {
    return std::log(err2 / err1) / std::log(double(n1) / n2);
}

void print_n(int n1, int n2, int n3) {
    std::cout << "                n: " << std::setw(10) << n1 << std::setw(12) << n2 << std::setw(12) << n3 << "\n";
}

void print_errors(std::string title, int n1, int n2, int n3, double err1, double err2, double err3) {
    std::cout << title << ":   ";
    std::cout << std::scientific << std::setprecision(2);
    std::cout << err1 << ",   " << err2 << ",   " << err3 << ";\t\tacc. ord.: ";
    double acc1 = acc_ord(n1, err1, n2, err2);
    double acc2 = acc_ord(n1, err1, n3, err3);
    double acc3 = acc_ord(n2, err2, n3, err3);
    double acc_min = std::min(acc1, std::min(acc2, acc3));
    double acc_max = std::max(acc1, std::min(acc2, acc3));
    std::cout << std::fixed;
    std::cout << acc_min << " -- " << acc_max << "\n";
}

// Интегрирование гладкой функции
void test1() {
    int n1, n2, n3;

    if (true) {
        std::cout << "Интегрирование гладкой функции по квадратному элементу:\n";

        Vector3d v11 = {0.0, 0.0, 0.0};
        Vector3d v21 = {1.0, 0.0, 0.0};
        Vector3d v12 = {0.0, 1.0, 0.0};
        Vector3d v22 = {1.0, 1.0, 0.0};

        Quad quad(v11, v21, v12, v22);

        double I = test_integral1();
        std::cout << "  Интеграл: " << std::fixed << std::setprecision(6) << I << "\n";

        n1 = 10;
        n2 = 100;
        n3 = 500;
        print_n(n1, n2, n3);
        print_errors("  Err(low   acc.)", n1, n2, n3,
                     std::abs(quad.integrate_low(test_func1, n1) - I),
                     std::abs(quad.integrate_low(test_func1, n2) - I),
                     std::abs(quad.integrate_low(test_func1, n3) - I));
        print_errors("  Err(mid   acc.)", n1, n2, n3,
                     std::abs(quad.integrate_mid(test_func1, n1) - I),
                     std::abs(quad.integrate_mid(test_func1, n2) - I),
                     std::abs(quad.integrate_mid(test_func1, n3) - I));

        n1 = 10;
        n2 = 20;
        n3 = 50;
        print_n(n1, n2, n3);
        print_errors("  Err(high  acc.)", n1, n2, n3,
                     std::abs(quad.integrate_high(test_func1, n1) - I),
                     std::abs(quad.integrate_high(test_func1, n2) - I),
                     std::abs(quad.integrate_high(test_func1, n3) - I));

        n1 = 2;
        n2 = 4;
        n3 = 6;
        print_n(n1, n2, n3);
        print_errors("  Err(extra acc.)", n1, n2, n3,
                     std::abs(quad.integrate_extra(test_func1, n1) - I),
                     std::abs(quad.integrate_extra(test_func1, n2) - I),
                     std::abs(quad.integrate_extra(test_func1, n3) - I));

        std::cout << "\n";
    }

    if (true) {
        std::cout << "Интегрирование гладкой функции по треугольному элементу:\n";

        Vector3d v1 = {0.0, 0.0, 0.0};
        Vector3d v2 = {1.0, 0.0, 0.0};
        Vector3d v3 = {1.0, 1.0, 0.0};

        Triangle tri(v1, v2, v3);

        // интеграл по треугольнику (функция симметрична)
        double I = 0.5 * test_integral1();
        std::cout << "  Интеграл: " << std::fixed << std::setprecision(6) << I << "\n";

        n1 = 10;
        n2 = 100;
        n3 = 500;
        print_n(n1, n2, n3);
        print_errors("  Err(low   acc.)", n1, n2, n3,
                     std::abs(tri.integrate_low(test_func1, n1) - I),
                     std::abs(tri.integrate_low(test_func1, n2) - I),
                     std::abs(tri.integrate_low(test_func1, n3) - I));
        print_errors("  Err(mid   acc.)", n1, n2, n3,
                     std::abs(tri.integrate_mid(test_func1, n1) - I),
                     std::abs(tri.integrate_mid(test_func1, n2) - I),
                     std::abs(tri.integrate_mid(test_func1, n3) - I));

        n1 = 10;
        n2 = 20;
        n3 = 50;
        print_n(n1, n2, n3);
        print_errors("  Err(high  acc.)", n1, n2, n3,
                     std::abs(tri.integrate_high(test_func1, n1) - I),
                     std::abs(tri.integrate_high(test_func1, n2) - I),
                     std::abs(tri.integrate_high(test_func1, n3) - I));

        n1 = 2;
        n2 = 4;
        n3 = 6;
        print_n(n1, n2, n3);
        print_errors("  Err(extra acc.)", n1, n2, n3,
                     std::abs(tri.integrate_extra(test_func1, n1) - I),
                     std::abs(tri.integrate_extra(test_func1, n2) - I),
                     std::abs(tri.integrate_extra(test_func1, n3) - I));

        std::cout << "\n";
    }

    std::cout << "Интегрирование гладкой функции по скошенному четырехугольнику:\n";
    if (true) {

        Vector3d v1 = {-0.2, -0.4, 0.0};
        Vector3d v2 = {1.2, 0.0, 0.0};
        Vector3d v3 = {0.0, 1.0, 0.0};
        Vector3d v4 = {0.8, 1.1, 0.0};

        Quad quad(v1, v2, v3, v4);
        Triangle tri1(v1, v2, v4);
        Triangle tri2(v1, v4, v3);
        Triangle tri3(v1, v2, v3);
        Triangle tri4(v2, v4, v3);

        // Точный интеграл неизвестен, считаем приближенные тремя способами
        double I1 = quad.integrate_high(test_func1, 90); // Погр ~ 1.1e-16
        double I2 = tri1.integrate_high(test_func1, 80) + // Погр ~ 6.1e-16
                    tri2.integrate_high(test_func1, 80);
        double I3 = tri3.integrate_high(test_func1, 80) + // Погр ~ 6.1e-16
                    tri4.integrate_high(test_func1, 80);

        std::cout << "  Интеграл: " << std::fixed << std::setprecision(6) << I1 << " " << I2 << " " << I3 << "; \t";
        std::cout << "Отклонение: " << std::scientific << std::setprecision(2)
                  << std::max(I1, std::max(I2, I3)) - std::min(I1, std::min(I2, I3)) << "\n";
        double I = (I1 + I2 + I3) / 3.0;

        n1 = 10;
        n2 = 100;
        n3 = 500;
        print_n(n1, n2, n3);
        print_errors("  Err(low   acc.)", n1, n2, n3,
                     std::abs(quad.integrate_low(test_func1, n1) - I),
                     std::abs(quad.integrate_low(test_func1, n2) - I),
                     std::abs(quad.integrate_low(test_func1, n3) - I));
        print_errors("  Err(mid   acc.)", n1, n2, n3,
                     std::abs(quad.integrate_mid(test_func1, n1) - I),
                     std::abs(quad.integrate_mid(test_func1, n2) - I),
                     std::abs(quad.integrate_mid(test_func1, n3) - I));

        n1 = 10;
        n2 = 20;
        n3 = 50;
        print_n(n1, n2, n3);
        print_errors("  Err(high  acc.)", n1, n2, n3,
                     std::abs(quad.integrate_high(test_func1, n1) - I),
                     std::abs(quad.integrate_high(test_func1, n2) - I),
                     std::abs(quad.integrate_high(test_func1, n3) - I));

        n1 = 2;
        n2 = 4;
        n3 = 6;
        print_n(n1, n2, n3);
        print_errors("  Err(extra acc.)", n1, n2, n3,
                     std::abs(quad.integrate_extra(test_func1, n1) - I),
                     std::abs(quad.integrate_extra(test_func1, n2) - I),
                     std::abs(quad.integrate_extra(test_func1, n3) - I));

        std::cout << "\n";
    }
    
}

// elem -- ячейка
// test_func -- Характ функция
// inside -- та же функция, но со значениями true/false
// I -- значение интеграла по ячейке
template<class Elem>
void volume_test(const Elem& elem,
              const std::function<double(const Vector3d&)>& test_func,
              const std::function<bool(const Vector3d&)>& inside, double I) {
    int n1 = 10;
    int n2 = 100;
    int n3 = 1000;

    std::cout << "  Интеграл: " << std::fixed << std::setprecision(6) << I << "\n";

    print_n(n1, n2, n3);
    print_errors("  Err(low   acc.)", n1, n2, n3,
                 std::abs(elem.integrate_low(test_func, n1) - I),
                 std::abs(elem.integrate_low(test_func, n2) - I),
                 std::abs(elem.integrate_low(test_func, n3) - I));
    print_errors("  Err(mid   acc.)", n1, n2, n3,
                 std::abs(elem.integrate_mid(test_func, n1) - I),
                 std::abs(elem.integrate_mid(test_func, n2) - I),
                 std::abs(elem.integrate_mid(test_func, n3) - I));
    print_errors("  Err(high  acc.)", n1, n2, n3,
                 std::abs(elem.integrate_high(test_func, n1) - I),
                 std::abs(elem.integrate_high(test_func, n2) - I),
                 std::abs(elem.integrate_high(test_func, n3) - I));
    print_errors("  Err(extra acc.)", n1, n2, n3,
                 std::abs(elem.integrate_extra(test_func, n1) - I),
                 std::abs(elem.integrate_extra(test_func, n2) - I),
                 std::abs(elem.integrate_extra(test_func, n3) - I));
    print_errors("  Vol. frac. func", n1, n2, n3,
                 std::abs(elem.volume_fraction(inside, n1 * n1) * elem.area() - I),
                 std::abs(elem.volume_fraction(inside, n2 * n2) * elem.area() - I),
                 std::abs(elem.volume_fraction(inside, n3 * n3) * elem.area() - I));

    std::cout << "\n";
}

// Отсечение от элемента части круга
void test2() {
    const double R = 0.5;
    const double R2 = R * R;

    auto inside = [&R](const Vector3d &v) -> bool {
        return v.norm() < R;
    };

    auto test_func = [&R](const Vector3d &v) -> double {
        return v.norm() < R ? 1.0 : 0.0;
    };

    // квадрат
    Vector3d v11 = {0.0, 0.0, 0.0};
    Vector3d v21 = {1.0, 0.0, 0.0};
    Vector3d v12 = {0.0, 1.0, 0.0};
    Vector3d v22 = {1.0, 1.0, 0.0};

    Quad quad(v11, v21, v12, v22);

    // Четверть круга
    double I1 = 0.25 * M_PI * R2;

    // скошеный четырехугольник
    Vector3d p11 = {0.0, 0.0, 0.0};
    Vector3d p21 = {1.0, 0.2, 0.0};
    Vector3d p12 = {0.0, 1.0, 0.0};
    Vector3d p22 = {1.2, 1.4, 0.0};

    Quad skewed(p11, p21, p12, p22);

    Triangle triangle(p11, p21, p12);

    // угол при p11
    Vector3d a = p12 - p11;
    Vector3d b = p21 - p11;
    double alpha = std::acos(a.dot(b) / (a.norm() * b.norm()));
    double I2 = 0.5 * alpha * R2;

    std::cout << "Объемная доля части круга на квадратном элементе:\n";
    volume_test(quad, test_func, inside, I1);

    std::cout << "Объемная доля части круга на треугольном элементе:\n";
    volume_test(triangle, test_func, inside, I2);

    std::cout << "Объемная доля части круга на скошеном квадратном элементе:\n";
    volume_test(skewed, test_func, inside, I2);
}

// Отсечение от элемента полуплоскости
void test3() {
    // Внешняя нормаль плоскости
    Vector3d n = (Vector3d{0.7, 0.3, 0.0}).normalized();

    // Точка прямой
    Vector3d p = Vector3d{0.2523, 0.3863, 0.0};

    auto inside = [&p, &n](const Vector3d &v) -> bool {
        return (v - p).dot(n) < 0.0;
    };

    auto test_func = [&p, &n](const Vector3d &v) -> double {
        return (v - p).dot(n) < 0.0 ? 1.0 : 0.0;
    };

    // квадрат
    Vector3d v11 = {0.0, 0.0, 0.0};
    Vector3d v21 = {1.0, 0.0, 0.0};
    Vector3d v12 = {0.0, 1.0, 0.0};
    Vector3d v22 = {1.0, 1.0, 0.0};

    Quad quad(v11, v21, v12, v22);

    // Отсечение прямой
    double I1 = Polygon({v11, v21, v12, v22}, true).clip_area(p, n);

    // скошеный четырехугольник
    Vector3d p11 = {0.0, 0.0, 0.0};
    Vector3d p21 = {1.0, 0.2, 0.0};
    Vector3d p12 = {0.0, 1.0, 0.0};
    Vector3d p22 = {1.2, 1.4, 0.0};

    Quad skewed(p11, p21, p12, p22);

    // Отсечение прямой
    double I2 = Polygon({p11, p21, p12, p22}, true).clip_area(p, n);

    Triangle triangle(p11, p21, p12);

    // Отсечение прямой
    double I3 = Polygon({p11, p21, p12}, true).clip_area(p, n);

    std::cout << "Объемная доля полуплоскости на квадратном элементе:\n";
    volume_test(quad, test_func, inside, I1);

    std::cout << "Объемная доля полуплоскости на треугольном элементе:\n";
    volume_test(triangle, test_func, inside, I2);

    std::cout << "Объемная доля полуплоскости на скошеном квадратном элементе:\n";
    volume_test(skewed, test_func, inside, I2);
}

int main() {
    test1();
    test2();
    test3();

    return 0;
}