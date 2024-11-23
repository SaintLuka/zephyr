#include <cmath>

#include <zephyr/geom/curves/interpolant.h>

namespace zephyr::geom::curves {

void fit_into_period(double x_min, double x_max, double& x) {
    double period = x_max - x_min;
    x -= std::floor((x - x_min) / period) * period;
}

int find_segment(const std::vector<double>& xs, double x) {
    int n = std::distance(xs.begin(), std::lower_bound(xs.begin(), xs.end(), x));
    return std::max(0, n - 1);
}

int monotonic(const std::vector<double>& arr) {
    if (arr.size() < 2) {
        return 1;
    }
    double sign = arr[1] - arr[0];
    if (sign == 0.0) {
        return 0;
    }
    for (size_t i = 2; i < arr.size(); ++i) {
        if ((arr[i] - arr[i - 1]) * sign <= 0.0) {
            return 0;
        }
    }
    return sign > 0.0 ? +1 : -1;
}

std::vector<double> chord_parametrization(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs)
{
    std::vector<double> ts(xs.size());
    ts[0] = 0.0;
    for (size_t i = 1; i < xs.size(); ++i) {
        Vector3d p1 = {xs[i], ys[i], zs[i]};
        Vector3d p2 = {xs[i - 1], ys[i - 1], zs[i - 1]};
        ts[i] = ts[i - 1] + (p2 - p1).norm();
    }

    double L = ts.back();
    for (size_t i = 1; i < ts.size(); ++i) {
        ts[i] /= L;
    }
    return ts;
}

std::vector<double> chebyshev_parametrization(int n) {
    std::vector<double> ts(n);
    for (int i = 0; i < n; ++i) {
        ts[i] = 0.5 * (1.0 - std::cos((2 * i + 1) * M_PI / (2 * n)) / std::cos(M_PI / (2 * n)));
    }
    return ts;
}

// ============================================================================
//                               INTERPOLANT
// ============================================================================

double Interpolant::x_min() const {
    return m_xs.front();
}

double Interpolant::x_max() const {
    return m_xs.back();
}

const std::vector<double>& Interpolant::xs() const {
    return m_xs;
}

const std::vector<double>& Interpolant::ys() const {
    return m_ys;
}

std::vector<double> Interpolant::xs(int N) const {
    return xs(N, x_min(), x_max());
}

std::vector<double> Interpolant::ys(int N) const {
    return ys(N, x_min(), x_max());
}

std::vector<double> Interpolant::xs(int N, double x1, double x2) const {
    return geom::linspace(x1, x2, N);
}

std::vector<double> Interpolant::ys(int N, double x1, double x2) const {
    std::vector<double> x = xs(N, x1, x2);
    std::vector<double> y(x.size());
    for (int i = 0; i < int(x.size()); ++i) {
        y[i] = get(x[i]);
    }
    return y;
}

// ============================================================================
//                             PARAMTERIC INTERPOLANT
// ============================================================================

double PInterpolant::t_min() const {
    return m_ts.front();
}

double PInterpolant::t_max() const {
    return m_ts.back();
}

Vector3d PInterpolant::get(double t) const {
    return {x(t), y(t), z(t)};
}

std::vector<Vector3d> PInterpolant::vs() const {
    return zip(m_xs, m_ys, m_zs);
}

const std::vector<double> &PInterpolant::ts() const {
    return m_ts;
}

const std::vector<double> &PInterpolant::xs() const {
    return m_xs;
}

const std::vector<double> &PInterpolant::ys() const {
    return m_ys;
}

const std::vector<double> &PInterpolant::zs() const {
    return m_zs;
}

std::vector<Vector3d> PInterpolant::vs(int N) const {
    return vs(N, t_min(), t_max());
}

std::vector<double> PInterpolant::ts(int N) const {
    return ts(N, t_min(), t_max());
}

std::vector<double> PInterpolant::xs(int N) const {
    return xs(N, t_min(), t_max());
}

std::vector<double> PInterpolant::ys(int N) const {
    return ys(N, t_min(), t_max());
}

std::vector<double> PInterpolant::zs(int N) const {
    return zs(N, t_min(), t_max());
}

std::vector<Vector3d> PInterpolant::vs(int N, double t1, double t2) const {
    std::vector<double> t = ts(N, t1, t2);
    std::vector<Vector3d> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = v(t[i]);
    }
    return res;
}

std::vector<double> PInterpolant::ts(int N, double t1, double t2) const {
    return geom::linspace(t1, t2, N);
}

std::vector<double> PInterpolant::xs(int N, double t1, double t2) const {
    std::vector<double> t = ts(N, t1, t2);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = x(t[i]);
    }
    return res;
}

std::vector<double> PInterpolant::ys(int N, double t1, double t2) const {
    std::vector<double> t = ts(N, t1, t2);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = y(t[i]);
    }
    return res;
}

std::vector<double> PInterpolant::zs(int N, double t1, double t2) const {
    std::vector<double> t = ts(N, t1, t2);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = z(t[i]);
    }
    return res;
}

} // namespace zephyr::geom::curves