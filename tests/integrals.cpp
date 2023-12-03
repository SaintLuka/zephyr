#include <iostream>

#include <zephyr/geom/maps.h>

using namespace zephyr::geom;

double test_func(const Vector3d& v) {
    return std::sin(10.0 * v.x() * v.x()) + std::cos(10.0 * v.x() + 5 * v.y());
}

int main() {
    Vector3d v11 = {0.0, 0.0, 0.0};
    Vector3d v21 = {1.0, 0.0, 0.0};
    Vector3d v12 = {0.0, 1.0, 0.0};
    Vector3d v22 = {1.0, 1.0, 0.0};

    Quad q(v11, v21, v12, v22);

    double I = 0.225228774825515581034;

    int n = 100;
    std::cout << "Integrate(r): " << std::abs(q.integrate_r(test_func, n) - I) << "\n";
    std::cout << "Integrate(t): " << std::abs(q.integrate_t(test_func, n) - I) << "\n";
    std::cout << "Integrate(g): " << std::abs(q.integrate_g(test_func, n) - I) << "\n";
    std::cout << "Integrate(s): " << std::abs(q.integrate_s(test_func, n) - I) << "\n";
    std::cout << "Integrate(m): " << std::abs(q.integrate_m(test_func, n) - I) << "\n";

    return 0;
}