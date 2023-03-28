#include <algorithm>

#include <zephyr/mesh/generator/vertex.h>
#include <zephyr/mesh/generator/curve/cubic.h>

namespace zephyr { namespace mesh { namespace generator {


inline double sqr(double x) {
    return x * x;
}

Cubic::Cubic(const Vector3d &v1, const Vector3d &v2,
      const Vector3d &v3, const Vector3d &v4)
      : v1(v1), v2(v2), v3(v3), v4(v4) {

    size_t N = 10;
    for (size_t i = 0; i <= N; ++i) {
        double t = (2.0*i)/N - 1.0;
        m_points.emplace_back(get_V(t));
    }
}

Curve::Ptr Cubic::create(const Vector3d &v1, const Vector3d &v2,
                         const Vector3d &v3, const Vector3d &v4) {
    return std::make_shared<Cubic>(v1, v2, v3, v4);
}

Curve::Ptr Cubic::create(BaseVertex::Ref v1, BaseVertex::Ref v2,
                         BaseVertex::Ref v3, BaseVertex::Ref v4) {
    return std::make_shared<Cubic>(v1->v(), v2->v(), v3->v(), v4->v());
}

double Cubic::get_X(double t) const {
    double A = -1.0 / 16.0 * (t - 1) * (3 * t - 1) * (3 * t + 1);
    double B = +9.0 / 16.0 * (t - 1) * (3 * t - 1) * (t + 1);
    double C = -9.0 / 16.0 * (t - 1) * (3 * t + 1) * (t + 1);
    double D = +1.0 / 16.0 * (3 * t - 1) * (3 * t + 1) * (t + 1);
    return v1.x() * A + v2.x() * B + v3.x() * C + v4.x() * D;
}

double Cubic::get_Y(double t) const {
    double A = -1.0 / 16.0 * (t - 1) * (3 * t - 1) * (3 * t + 1);
    double B = +9.0 / 16.0 * (t - 1) * (3 * t - 1) * (t + 1);
    double C = -9.0 / 16.0 * (t - 1) * (3 * t + 1) * (t + 1);
    double D = +1.0 / 16.0 * (3 * t - 1) * (3 * t + 1) * (t + 1);
    return v1.y() * A + v2.y() * B + v3.y() * C + v4.y() * D;
}

Vector3d Cubic::get_V(double t) const {
    double A = -1.0 / 16.0 * (t - 1) * (3 * t - 1) * (3 * t + 1);
    double B = +9.0 / 16.0 * (t - 1) * (3 * t - 1) * (t + 1);
    double C = -9.0 / 16.0 * (t - 1) * (3 * t + 1) * (t + 1);
    double D = +1.0 / 16.0 * (3 * t - 1) * (3 * t + 1) * (t + 1);
    return v1 * A + v2 * B + v3 * C + v4 * D;
}

Vector3d Cubic::get_T(double t) const {
    return Curve::get_T(t);
}

Vector3d Cubic::get_N(double t) const {
    return Curve::get_N(t);
}

Vector3d Cubic::projection(const Vector3d &v) const {
    std::vector<double> dists(m_points.size());
    for (size_t i = 0; i < m_points.size(); ++i) {
        dists[i] = (v - m_points[i]).squaredNorm();
    }

    size_t idx = std::distance(dists.begin(), std::min_element(dists.begin(), dists.end()));

    double t = 2.0 * idx / (m_points.size() - 1.0) - 1.0;

    double eps = 1.0;
    int counter = 0;
    while (eps > 1.0e-8 && counter < 20) {
        double dt = -obj_func(v, t) / obj_func_deriv(v, t);
        t += dt;

        if (std::abs(t) >= 1.0) {
            return t < 0.0 ? v1 : v4;
        }

        eps = std::abs(dt);
        ++counter;
    }
    return get_V(t);
}

double Cubic::dX(double t) const {
    double A = (-27 * sqr(t) + 18 * t + 1) / 16.0;
    double B = (+81 * sqr(t) - 18 * t - 27) / 16.0;
    double C = (-81 * sqr(t) - 18 * t + 27) / 16.0;
    double D = (+27 * sqr(t) + 18 * t - 1) / 16.0;

    return v1.x() * A + v2.x() * B + v3.x() * C + v4.x() * D;
}

double Cubic::dY(double t) const {
    double A = (-27 * sqr(t) + 18 * t + 1) / 16.0;
    double B = (+81 * sqr(t) - 18 * t - 27) / 16.0;
    double C = (-81 * sqr(t) - 18 * t + 27) / 16.0;
    double D = (+27 * sqr(t) + 18 * t - 1) / 16.0;

    return v1.y() * A + v2.y() * B + v3.y() * C + v4.y() * D;
}

double Cubic::d2X(double t) const {
    double A = (-27 * t + 9) / 8.0;
    double B = (+81 * t - 9) / 8.0;
    double C = (-81 * t - 9) / 8.0;
    double D = (+27 * t + 9) / 8.0;

    return v1.x() * A + v2.x() * B + v3.x() * C + v4.x() * D;
}

double Cubic::d2Y(double t) const {
    double A = (-27 * t + 9) / 8.0;
    double B = (+81 * t - 9) / 8.0;
    double C = (-81 * t - 9) / 8.0;
    double D = (+27 * t + 9) / 8.0;

    return v1.y() * A + v2.y() * B + v3.y() * C + v4.y() * D;
}

double Cubic::obj_func(const Vector3d& v, double t) const {
    return
        (v.x() - get_X(t)) * dX(t) +
        (v.y() - get_Y(t)) * dY(t);
}

double Cubic::obj_func_deriv(const Vector3d& v, double t) const {
    return
        - sqr(dX(t)) - sqr(dY(t)) +
        (v.x() - get_X(t)) * d2X(t) +
        (v.y() - get_Y(t)) * d2Y(t);
}

} // namespace generator
} // namespace mesh
} // namespace zephyr
