#include <iostream>

#include <zephyr/math/sle/tridiagonal.h>
#include <zephyr/math/calc/derivatives.h>
#include <zephyr/geom/curves/cubic_spline.h>

namespace zephyr::geom::curves {

// return a + b t + c t^2 + d t^3
// a + t (b + t (c + d t)
inline double cubic_segment(double a, double b, double c, double d, double t) {
    return a + t * (b + t * (c + t * d));
}

inline double cubic_segment(const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d,
        int n, double x) {
    return cubic_segment(ys[n], b[n], c[n], d[n], x - xs[n]);
}

std::vector<double> get_c(
        const std::vector<double>& xs,
        const std::vector<double>& ys)
{
    using namespace zephyr::math;
    using tridiagonal::array;

    // Количество частей сплайна,
    // на единицу меньше числа узлов
    int n = int(xs.size()) - 1;

    // Находим коэффициенты сплайнов,
    // пропорциональные вторым производным.
    std::vector<double> c(n + 1);

    array A(n + 1), B(n + 1), C(n + 1);

    A[0] = A[n] = 0.0;
    B[0] = B[n] = 1.0;
    C[0] = C[n] = 0.0;

    c[0] = c[n] = 0.0;
    for (int i = 1; i < n; ++i) {
        // Коэффициенты матрицы для прогонки
        int ip = i + 1;
        int im = i - 1;

        double hp = xs[ip] - xs[i];
        double hm = xs[i] - xs[im];

        A[i] = hm;
        B[i] = 2.0 * (hm + hp);
        C[i] = hp;

        // Правые части (потом решение)
        c[i] = 3.0 * ((ys[ip] - ys[i]) / hp - (ys[i] - ys[im]) / hm);
    }

    // Решение СЛАУ записывается в последний массив
    tridiagonal::solve_l(A, B, C, c);
    return c;
}

std::vector<double> get_c_periodic(
        const std::vector<double>& xs,
        const std::vector<double>& ys)
{
    using namespace zephyr::math;
    using tridiagonal::array;

    // Количество частей сплайна (на единицу меньше числа узлов)
    int n = int(xs.size()) - 1;

    double period = xs.back() - xs.front();

    // Находим коэффициенты сплайнов,
    // пропорциональные вторым производным.
    std::vector<double> c(n + 1);

    c.resize(n);

    array A(n), B(n), C(n);

    for (int i = 0; i < n; ++i) {
        // Коэффициенты матрицы для прогонки
        int ip = (i + 1) % n;
        int im = (i - 1 + n) % n;

        double hp = xs[ip] - xs[i];
        double hm = xs[i] - xs[im];

        hp = hp > 0.0 ? hp : hp + period;
        hm = hm > 0.0 ? hm : hm + period;

        A[i] = hm;
        B[i] = 2.0 * (hm + hp);
        C[i] = hp;

        // Правые части (потом решение)
        c[i] = 3.0 * ((ys[ip] - ys[i]) / hp - (ys[i] - ys[im]) / hm);
    }

    // Решение СЛАУ записывается в последний массив
    tridiagonal::solve_cyclic_l(A, B, C, c);

    c.push_back(c.front());
    return c;
}

// Построить кубический сплайн. b, c, d -- выходные параметры
void basic_build(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        SplineBound L, SplineBound R,
        std::vector<double>& b,
        std::vector<double>& c,
        std::vector<double>& d) {

    using namespace zephyr::math;

    bool periodic = L == SplineBound::Periodic ||
                    R == SplineBound::Periodic;

    // Рассматриваем тривиальные случаи
    if (xs.size() < 2) {
        std::string message = "Build CubicSpline error: can't build spline on "
                              + std::to_string(xs.size()) + " points.";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    else if (xs.size() == 2) {
        b = { xs[1] - xs[0] };
        c = { 0.0 };
        d = { 0.0 };
        return;
    }

    // Проверим периодику, если не замкнута
    if (periodic) {
        if (ys.front() != ys.back()) {
            std::cout << xs.front() << " " << xs.back() << "\n";
            std::cout << ys.front() << " " << ys.back() << "\n";
            std::string message = "Periodic cubic spline error: "
                                  "first and last values are different";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }
    }

    // Находим коэффициенты сплайнов, пропорциональные вторым производным
    if (!periodic) {
        c = get_c(xs, ys);
    }
    else {
        c = get_c_periodic(xs, ys);
    }

    // Количество частей сплайна (на единицу меньше числа узлов)
    int n = int(xs.size()) - 1;

    // Находим коэффициенты сплайна, пропорциональные третьим производным
    d.resize(n);

    for (int i = 0; i < n; ++i) {
        d[i] = (c[i + 1] - c[i]) / (3.0 * (xs[i + 1] - xs[i]));
    }

    // Находим коэффициенты сплайна, пропорциональные первым производным
    b.resize(n);

    for (int i = 0; i < n; ++i) {
        b[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) -
               (xs[i + 1] - xs[i]) * (c[i + 1] + 2 * c[i]) / 3.0;
    }

    // Удаляем последний вспомогательный элемент
    c.pop_back();
}

static double basic_get(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d,
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
                return cubic_segment(ys[0], b[0], 0.0, 0.0, x - xs[0]);
            case SplineBound::Warning: {
                std::cerr << "CubicSpline: parameter out of bounds\n";
                return ys.front();
            }
            default:
                throw std::runtime_error("CubicSpline error: bad case");
        }
    } else if (x > xs.back()) {
        switch (R) {
            case SplineBound::Crop:
                return ys.back();
            case SplineBound::Free: {
                int n = int(xs.size()) - 1;
                double dx = xs[n] - xs[n - 1];
                double B = b[n - 1] + 2 * c[n - 1] * dx + 3 * d[n - 1] * dx * dx;
                return cubic_segment(ys[n], B, 0.0, 0.0, x - xs[n]);
            }
            case SplineBound::Warning: {
                std::cerr << "CubicSpline: parameter out of bounds\n";
                return ys.back();
            }
            default:
                throw std::runtime_error("CubicSpline error: bad case");
        }
    } else {
        // Нормальная ситуация
        int n = find_segment(xs, x);
        return cubic_segment(xs, ys, b, c, d, n, x);
    }
}


CubicSpline::CubicSpline(
        const std::vector<double>& x,
        const std::vector<double>& y,
        SplineBound left,
        SplineBound right) {
    m_left  = left;
    m_right = right;
    m_xs    = x;
    m_ys    = y;

    if (m_xs.size() != m_ys.size()) {
        throw std::runtime_error("CubicSpline error: different sizes (ts, xs)");
    }
    if (m_xs.size() < 2) {
        throw std::runtime_error("CubicSpline error: need at least two points");
    }

    if (m_left == SplineBound::Periodic || m_right == SplineBound::Periodic) {
        m_left = m_right = SplineBound::Periodic;

        if (m_ys.front() != m_ys.back()) {
            throw std::runtime_error("Periodic cubic spline error: "
                                     "first and last values are different");
        }
    }

    if (!monotonic(m_xs)) {
        throw std::runtime_error("CubicSpline error: not monotonic arguments");
    }

    basic_build(m_xs, m_ys, m_left, m_right, m_b, m_c, m_d);
}

double CubicSpline::get(double x) const {
    return basic_get(m_xs, m_ys, m_b, m_c, m_d, m_left, m_right, x);
}

PCubicSpline::PCubicSpline(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        SplineBound left, SplineBound right,
        Parametrization param) {
    std::vector<double> zs(xs.size(), 0.0);
    build(xs, ys, zs, left, right, param);
}

PCubicSpline::PCubicSpline(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right,
        Parametrization param) {
    build(xs, ys, zs, left, right, param);
}

PCubicSpline::PCubicSpline(
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

inline double maxmod(double x, double y, double z) {
    return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
}

void PCubicSpline::build(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs,
        SplineBound left, SplineBound right,
        Parametrization param) {

    if (xs.size() < 2) {
        throw std::runtime_error("PCubicSpline error: need at least two points");
    }

    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        throw std::runtime_error("PCubicSpline error: sizes mismatch");
    }

    int n = xs.size();
    int N = n;

    bool periodic = left == SplineBound::Periodic ||
                    right == SplineBound::Periodic;

    m_left  = periodic ? SplineBound::Periodic : left;
    m_right = periodic ? SplineBound::Periodic : right;

    if (periodic) {
        double L = 0.0;
        for (int i = 1; i < n; ++i) {
            double dist = maxmod(xs[i] - xs[i - 1],
                                 ys[i] - ys[i - 1],
                                 zs[i] - zs[i - 1]);
            L = std::max(dist, L);
        }

        // Добавим один для замыкания
        double dist = maxmod(xs.back() - xs.front(),
                             ys.back() - ys.front(),
                             zs.back() - zs.front());
        if (dist >= 1.0e-8 * L) {
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
    if (periodic) {
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
            std::string message = "CubicSpline error: Unknown parametrization type";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }
    }

    basic_build(m_ts, m_xs, m_left, m_right, m_xb, m_xc, m_xd);
    basic_build(m_ts, m_ys, m_left, m_right, m_yb, m_yc, m_yd);
    basic_build(m_ts, m_zs, m_left, m_right, m_zb, m_zc, m_zd);
}

double PCubicSpline::x(double t) const {
    return basic_get(m_ts, m_xs, m_xb, m_xc, m_xd, m_left, m_right, t);
}

double PCubicSpline::y(double t) const {
    return basic_get(m_ts, m_ys, m_yb, m_yc, m_yd, m_left, m_right, t);
}

double PCubicSpline::z(double t) const {
    return basic_get(m_ts, m_zs, m_zb, m_zc, m_zd, m_left, m_right, t);
}

} // namespace zephyr::geom::curves