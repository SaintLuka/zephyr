#include <iostream>

#include <zephyr/geom/curves/lagrange.h>

namespace zephyr::geom::curves {

static double lagrange_get(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
    int n = xs.size();
    double y = 0.0;
    for (int i = 0; i < n; ++i) {
        double L = 1.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                L *= (x - xs[j]) / (xs[i] - xs[j]);
            }
        }
        y += ys[i] * L;
    }
    return y;
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
                return lagrange_get(xs, ys, x);
            case SplineBound::Warning: {
                std::cerr << "Lagrange interpolant: parameter out of bounds\n";
                return ys.front();
            }
            default:
                throw std::runtime_error("Lagrange interpolant error: bad case");
        }
    } else if (x > xs.back()) {
        switch (R) {
            case SplineBound::Crop:
                return ys.back();
            case SplineBound::Free:
                return lagrange_get(xs, ys, x);
            case SplineBound::Warning: {
                std::cerr << "Lagrange interpolant: parameter out of bounds\n";
                return ys.back();
            }
            default:
                throw std::runtime_error("Lagrange interpolant error: bad case");
        }
    } else {
        // Нормальная ситуация
        return lagrange_get(xs, ys, x);
    }
}

Lagrange::Lagrange(
        const std::vector<double> &x,
        const std::vector<double> &y,
        SplineBound left, SplineBound right)
{
    m_left  = left;
    m_right = right;
    m_xs    = x;
    m_ys    = y;

    if (m_xs.size() != m_ys.size()) {
        throw std::runtime_error("Lagrage interpolant error: different sizes (ts, xs)");
    }
    if (m_xs.size() < 2) {
        throw std::runtime_error("Lagrage interpolant error: need at least two points");
    }

    if (m_left == SplineBound::Periodic || m_right == SplineBound::Periodic) {
        m_left = m_right = SplineBound::Periodic;

        if (m_ys.front() != m_ys.back()) {
            throw std::runtime_error("Periodic Lagrage interpolant error: "
                                     "first and last values are different");
        }
    }

    if (monotonic(m_xs) < 1) {
        throw std::runtime_error("Lagrange interpolant error: not monotonic arguments");
    }
}

double Lagrange::get(double x) const {
    return basic_get(m_xs, m_ys, m_left, m_right, x);
}




PLagrange::PLagrange(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        SplineBound left, SplineBound right,
        Parametrization param) {
    std::vector<double> zs(xs.size(), 0.0);
    build(xs, ys, zs, left, right, param);
}

PLagrange::PLagrange(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right,
        Parametrization param) {
    build(xs, ys, zs, left, right, param);
}

PLagrange::PLagrange(
        const std::vector<Vector3d> &vs,
        SplineBound left, SplineBound right,
        Parametrization param) {
    std::vector<double> xs(vs.size());
    std::vector<double> ys(vs.size());
    std::vector<double> zs(vs.size());
    for (size_t i = 0; i < vs.size(); ++i) {
        xs[i] = vs[i].x();
        ys[i] = vs[i].y();
        zs[i] = vs[i].z();
    }
    build(xs, ys, zs, left, right, param);
}

void PLagrange::build(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right,
        Parametrization param) {

    if (xs.size() < 2) {
        throw std::runtime_error("Lagrange interpolant error: need at least two points");
    }

    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        throw std::runtime_error("Lagrange interpolant error: sizes mismatch");
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

    // Замкнем периодический случай
    if (xs.size() < N) {
        m_xs.back() = xs.front();
        m_ys.back() = ys.front();
        m_zs.back() = zs.front();
    }

    // Параметризация
    switch (param) {
        case Parametrization::Uniform:
            m_ts = linspace(0.0, 1.0, m_xs.size());
            break;
        case Parametrization::Chord:
            m_ts = chord_parametrization(m_xs, m_ys, m_zs);
            break;
        case Parametrization::Chebyshev:
            m_ts = chebyshev_parametrization(m_xs.size());
            break;
        default: {
            std::string message = "Lagrange interpolant error: Unknown parametrization type";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }
    }
}

double PLagrange::x(double t) const {
    return basic_get(m_ts, m_xs, m_left, m_right, t);
}

double PLagrange::y(double t) const {
    return basic_get(m_ts, m_ys, m_left, m_right, t);
}

double PLagrange::z(double t) const {
    return basic_get(m_ts, m_zs, m_left, m_right, t);
}

} // namespace zephyr::geom::curves