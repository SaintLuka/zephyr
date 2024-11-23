#include <iostream>

#include <zephyr/geom/curves/linear_spline.h>

namespace zephyr::geom::curves {

// Линейная интерполяция на сегменте
inline double segment(double x1, double x2, double y1, double y2, double x) {
    return x1 != x2 ? y1 + (y2 - y1) / (x2 - x1) * (x - x1) : 0.5 * (y1 + y2);
}

// Линейная интерполяция на сегменте
inline double segment(const std::vector<double>& xs, const std::vector<double>& ys, int n, double x) {
    return segment(xs[n], xs[n + 1], ys[n], ys[n + 1], x);
}

static double basic_get(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        SplineBound L, SplineBound R,
        double x) {

    if (L == SplineBound::Periodic ||
        R == SplineBound::Periodic) {
        fit_into_period(xs.front(), xs.back(), x);
    }

    if (x < xs.front()) {
        switch (L) {
            case SplineBound::Crop:
                return ys.front();
            case SplineBound::Free:
                return segment(xs, ys, 0, x);
            case SplineBound::Warning: {
                std::cerr << "LinearSpline: parameter out of bounds\n";
                return ys.front();
            }
            default:
                throw std::runtime_error("LinearSpline error: bad case");
        }
    } else if (x > xs.back()) {
        switch (R) {
            case SplineBound::Crop:
                return ys.back();
            case SplineBound::Free:
                return segment(xs, ys, int(xs.size()) - 2, x);
            case SplineBound::Warning: {
                std::cerr << "LinearSpline: parameter out of bounds\n";
                return ys.back();
            }
            default:
                throw std::runtime_error("LinearSpline error: bad case");
        }
    } else {
        // Нормальная ситуация
        int n = find_segment(xs, x);
        return segment(xs, ys, n, x);
    }
}

LinearSpline::LinearSpline(
        const std::vector<double> &x,
        const std::vector<double> &y,
        SplineBound left, SplineBound right) 
{
    m_left  = left;
    m_right = right;
    m_xs    = x;
    m_ys    = y;
    
    if (m_xs.size() != m_ys.size()) {
        throw std::runtime_error("LinearSpline error: different sizes (ts, xs)");
    }
    if (m_xs.size() < 2) {
        throw std::runtime_error("LinearSpline error: need at least two points");
    }

    if (m_left == SplineBound::Periodic || m_right == SplineBound::Periodic) {
        m_left = m_right = SplineBound::Periodic;

        if (m_ys.front() != m_ys.back()) {
            throw std::runtime_error("Periodic linear spline error: "
                                     "first and last values are different");
        }
    }

    if (monotonic(m_xs) < 1) {
        throw std::runtime_error("LinearSpline error: not monotonic arguments");
    }
}

double LinearSpline::get(double x) const {
    return basic_get(m_xs, m_ys, m_left, m_right, x);
}


PLinearSpline::PLinearSpline(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        SplineBound left, SplineBound right) {
    std::vector<double> zs(xs.size(), 0.0);
    build(xs, ys, zs, left, right);
}
PLinearSpline::PLinearSpline(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right) {
    build(xs, ys, zs, left, right);
}

PLinearSpline::PLinearSpline(
        const std::vector<Vector3d> &vs,
        SplineBound left, SplineBound right) {
    std::vector<double> xs(vs.size());
    std::vector<double> ys(vs.size());
    std::vector<double> zs(vs.size());
    for (size_t i = 0; i < vs.size(); ++i) {
        xs[i] = vs[i].x();
        ys[i] = vs[i].y();
        zs[i] = vs[i].z();
    }
    build(xs, ys, zs, left, right);
}

void PLinearSpline::build(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right) {

    if (xs.size() < 2) {
        throw std::runtime_error("LinearSpline error: need at least two points");
    }

    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        throw std::runtime_error("PLinearSpline error: sizes mismatch");
    }

    int n = xs.size();
    int N = n;

    m_left  = left;
    m_right = right;

    if (left == SplineBound::Periodic || right == SplineBound::Periodic) {
        m_left = m_right = SplineBound::Periodic;

        // Добавим один для замыкания
        if (xs.front() != xs.back() ||
            ys.front() != zs.back() ||
            zs.front() != zs.back()) {
            N += 1;
        }
    }

    m_ts.resize(N);
    m_xs.resize(N);
    m_ys.resize(N);
    m_zs.resize(N);
    for (size_t i = 0; i < n; ++i) {
        m_xs[i] = xs[i];
        m_ys[i] = ys[i];
        m_zs[i] = zs[i];
    }

    m_ts[0] = 0.0;
    for (size_t i = 1; i < n; ++i) {
        Vector3d p1 = {xs[i], ys[i], zs[i]};
        Vector3d p2 = {xs[i - 1], ys[i - 1], zs[i - 1]};
        m_ts[i] = m_ts[i - 1] + (p2 - p1).norm();
    }

    // Замкнем периодический случай
    if (xs.size() < N) {
        m_xs.back() = xs.front();
        m_ys.back() = ys.front();
        m_zs.back() = zs.front();

        Vector3d p1 = {xs.front(), ys.front(), zs.front()};
        Vector3d p2 = {xs.back(), ys.back(), zs.back()};
        m_ts.back() = m_ts[n - 1] + (p2 - p1).norm();
    }

    double L = m_ts.back();
    for (size_t i = 1; i < m_ts.size(); ++i) {
        m_ts[i] /= L;
    }
}

double PLinearSpline::x(double t) const {
    return basic_get(m_ts, m_xs, m_left, m_right, t);
}

double PLinearSpline::y(double t) const {
    return basic_get(m_ts, m_ys, m_left, m_right, t);
}

double PLinearSpline::z(double t) const {
    return basic_get(m_ts, m_zs, m_left, m_right, t);
}

} // namespace zephyr::geom::curves