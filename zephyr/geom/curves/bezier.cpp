#include <zephyr/geom/curves/bezier.h>

namespace zephyr::geom::curves {

Bezier::Bezier(const std::vector<Vector3d>& vs)
    : m_xs(vs.size()), m_ys(vs.size()), m_zs(vs.size()) {

    for (size_t i = 0; i < vs.size(); ++i) {
        m_xs[i] = vs[i].x();
        m_ys[i] = vs[i].y();
        m_zs[i] = vs[i].z();
    }
}

Bezier::Bezier(std::vector<double> xs, std::vector<double> ys)
        : m_xs(std::move(xs)), m_ys(std::move(ys)) {

    if (m_xs.size() != m_ys.size()) {
        throw std::runtime_error("curves::Beizer error: different size");
    }

    m_zs = std::vector<double>(m_xs.size(), 0.0);
}

Bezier::Bezier(std::vector<double> xs, std::vector<double> ys, std::vector<double> zs)
        : m_xs(std::move(xs)), m_ys(std::move(ys)), m_zs(std::move(zs)) {

    if (m_xs.size() != m_ys.size() || m_xs.size() != m_ys.size()) {
        throw std::runtime_error("curves::Beizer error: different size");
    }
}

double binomial(int k, int n) {
    double res = 1.0;
    for (int i = 1; i <= std::min(k, n - k); ++i) {
        res *= (n - i + 1.0) / i;
    }
    return res;
}

Vector3d Bezier::get(double t) const {
    double x{0.0}, y{0.0}, z{0.0};

    int n = int(m_xs.size()) - 1;
    for (int k = 0; k <= n; ++k) {
        double bc = binomial(k, n) * std::pow(t, k) * std::pow(1.0 - t, n - k);

        x += m_xs[k] * bc;
        y += m_ys[k] * bc;
        z += m_zs[k] * bc;
    }
    return {x, y, z};
}

Vector3d Bezier::tangent(double t) const {
    return {NAN, NAN, NAN};
}

Vector3d Bezier::normal(double t, Vector3d p) const {
    return {NAN, NAN, NAN};
}

const std::vector<double>& Bezier::xs() const { return m_xs; }

const std::vector<double>& Bezier::ys() const { return m_ys; }

const std::vector<double>& Bezier::zs() const { return m_zs; }

std::vector<Vector3d> Bezier::vs() const {
    return zip(m_xs, m_ys, m_zs);
}

std::vector<double> Bezier::xs(int N) const {
    return xs(N, 0.0, 1.0);
}

std::vector<double> Bezier::ys(int N) const {
    return ys(N, 0.0, 1.0);
}

std::vector<double> Bezier::zs(int N) const {
    return zs(N, 0.0, 1.0);
}

std::vector<Vector3d> Bezier::vs(int N) const {
    return vs(N, 0.0, 1.0);
}

std::vector<double> Bezier::xs(int N, double t1, double t2) const {
    double h = (t2 - t1) / (N - 1.0);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = get(t1 + i * h).x();
    }
    return res;
}

std::vector<double> Bezier::ys(int N, double t1, double t2) const {
    double h = (t2 - t1) / (N - 1.0);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = get(t1 + i * h).y();
    }
    return res;
}

std::vector<double> Bezier::zs(int N, double t1, double t2) const {
    double h = (t2 - t1) / (N - 1.0);
    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = get(t1 + i * h).z();
    }
    return res;
}

std::vector<Vector3d> Bezier::vs(int N, double t1, double t2) const {
    double h = (t2 - t1) / (N - 1.0);
    std::vector<Vector3d> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = get(t1 + i * h);
    }
    return res;
}

} // namespace zephyr::geom::curves