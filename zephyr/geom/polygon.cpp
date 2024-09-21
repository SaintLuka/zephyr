#include <iomanip>
#include <algorithm>
#include <iostream>

#include <zephyr/math/calc/integrals.h>

#include <zephyr/geom/polygon.h>
#include <zephyr/geom/triangle.h>
#include <zephyr/geom/intersection.h>

namespace zephyr::geom {

inline bool bad_normal(const Vector3d& n) {
    return n.hasNaN() || n.squaredNorm() == 0.0;
}

Polygon::Polygon()
    : vs({}) , m_center(Vector3d::Zero()) {
}

Polygon::Polygon(int size)
    : vs(size, Vector3d::Zero()), m_center(Vector3d::Zero()) {

}

Polygon::Polygon(const Vector3d* buff, int size, bool sort) {
    vs.resize(size);
    std::memcpy((void *) vs.data(), (void *) buff, size * sizeof(Vector3d));
    setup_center();
    if (sort) {
        Polygon::sort();
    }
}

Polygon::Polygon(const std::vector<Vector3d>& vertices, bool sort)
        : vs(vertices) {
    setup_center();
    if (sort) {
        Polygon::sort();
    }
}

Polygon::Polygon(std::vector<Vector3d>&& vertices, bool sort)
    :vs(std::move(vertices)) {
    setup_center();
    if (sort) {
        Polygon::sort();
    }
}

Polygon::Polygon(std::initializer_list<Vector3d>&& list, bool sort)
    : vs(std::move(list)) {
    setup_center();
    if (sort) {
        Polygon::sort();
    }
}

void Polygon::setup_center() {
    m_center = {0.0, 0.0, 0.0};
    for (auto &v: vs) {
        m_center += v;
    }
    m_center /= vs.size();
}

void Polygon::reserve(int size) {
    vs.reserve(size);
}

void Polygon::set(int idx, const Vector3d& p) {
    m_center += (p - vs[idx]) / vs.size();
    vs[idx] = p;
}

void Polygon::operator+=(const Vector3d& p) {
    m_center = (vs.size() * m_center + p) / (vs.size() + 1);
    vs.emplace_back(p);
}

std::vector<double> Polygon::xs() const {
    if (empty()) {
        return {};
    }
    std::vector<double> res(size() + 1);
    for (int i = 0; i < size(); ++i) {
        res[i] = vs[i].x();
    }
    res[size()] = vs[0].x();
    return res;
}

std::vector<double> Polygon::ys() const {
    if (empty()) {
        return {};
    }
    std::vector<double> res(size() + 1);
    for (int i = 0; i < size(); ++i) {
        res[i] = vs[i].y();
    }
    res[size()] = vs[0].y();
    return res;
}
/*
void Polygon::print_numpy() const {
    std::cout << "x = np.array([";
    for (int i = 0; i < size(); ++i) {
        std::cout << vs[i].x() << ", ";
    }
    if (!empty()) {
        std::cout << vs[0].x();
    }
    std::cout << "])\n";

    std::cout << "y = np.array([";
    for (int i = 0; i < size(); ++i) {
        std::cout << vs[i].y() << ", ";
    }
    if (!empty()) {
        std::cout << vs[0].y();
    }
    std::cout << "])\n";
}
*/
Box Polygon::bbox() const {
    if (empty()) {
        return {Vector3d::Zero(), Vector3d::Zero()};
    }

    const double max = +std::numeric_limits<double>::infinity();
    const double min = -std::numeric_limits<double>::infinity();

    Vector3d vmin = {max, max, min};
    Vector3d vmax = {min, min, max};
    for (int i = 0; i < size(); ++i) {
        if (vs[i].x() < vmin.x()) {
            vmin.x() = vs[i].x();
        }
        else if (vs[i].x() > vmax.x()) {
            vmax.x() = vs[i].x();
        }

        if (vs[i].y() < vmin.y()) {
            vmin.y() = vs[i].y();
        }
        else if (vs[i].y() > vmax.y()) {
            vmax.y() = vs[i].y();
        }
    }

    return {vmin, vmax};
}

Vector3d Polygon::center() const {
    return m_center;
}

void Polygon::sort() {
    std::sort(vs.begin(), vs.end(),
              [this](const Vector3d &a, const Vector3d &b) -> bool {
                  double phi_a = std::atan2(a.y() - m_center.y(), a.x() - m_center.x());
                  double phi_b = std::atan2(b.y() - m_center.y(), b.x() - m_center.x());
                  return phi_a < phi_b;
              });
}

// vector product [v1, v2].z
double cross(const Vector3d& v1, const Vector3d& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

// vector product [v1 - c, v2 - c].z
double cross(const Vector3d& c, const Vector3d& v1, const Vector3d& v2) {
    return (v1.x() - c.x()) * (v2.y() - c.y()) -
           (v1.y() - c.y()) * (v2.x() - c.x());
}

bool Polygon::inside(const Vector3d& p) const {
    obj::ray ray = {p, Vector3d{0.123, 0.346433, 0.0}};

    // even-odd rule
    int counter = 0;
    for (int i = 0; i < size(); ++i) {
        obj::segment seg = {vs[i], vs[(i + 1) % size()]};
        if (intersection2D::exist(ray, seg)) {
            ++counter;
        }
    }
    return counter % 2 > 0;
}

double Polygon::area() const {
    double S = 0.0;
    for (int i = 0; i < size(); ++i) {
        const auto &v1 = vs[i];
        const auto &v2 = vs[(i + 1) % size()];

        S += cross(m_center, v1, v2);
    }
    return 0.5 * S;
}

Vector3d Polygon::centroid(double S) const {
    if (S == 0.0) {
        S = area();
    }

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        auto &v1 = vs[i];
        auto &v2 = vs[(i + 1) % size()];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * S);

    return C;
}

double Polygon::volume_as() const {
    return NAN;
}

double Polygon::clip_area(const std::function<bool(const Vector3d&)>& inside, int N) const {
    int n = std::max(1, N / size());

    double area = 0.0;
    for (int i = 0; i < size(); ++i) {
        const auto &v1 = vs[i];
        const auto &v2 = vs[(i + 1) % size()];

        Triangle tri(m_center, v1, v2);
        area += tri.clip_area(inside, n);
    }

    return area;
}

Polygon Polygon::clip(const Vector3d& p, const Vector3d& n) const {
    if (empty()) {
        return {};
    }

    std::vector<Vector3d> part;
    part.reserve(size() + 1);
    for (int i = 0; i < size(); ++i) {
        const auto &v1 = vs[i];
        const auto &v2 = vs[(i + 1) % size()];

        double dot1 = (v1 - p).dot(n);
        double dot2 = (v2 - p).dot(n);
        double Det = (v2 - v1).dot(n);

        // Проверяем вершину (если внутри)
        if (dot1 <= 1.0e-14) {
            part.emplace_back(v1);
        }

        // Проверяем пересечение
        if (dot1 * dot2 < 0.0) {
            part.emplace_back((dot2 * v1 - dot1 * v2) / Det);
        }
    }
    return Polygon(std::move(part));
}

struct AnS {
    double area;   ///< Площадь отсечения
    Vector3d p1;   ///< Точка сечения
    Vector3d p2;   ///< Точка сечения
    bool tri_in;   ///< Отсекается треугольник
    bool tri_out;  ///< Остается треугольник

    // Длина сечения
    double slice() const {
        return (p2 - p1).norm();
    }
};

// Найти площадь отсечения и длину сечения (если надо)
AnS clip_area_and_slice(const Polygon& poly, const Vector3d& p, const Vector3d& n, bool slice) {
    const Vector3d nan = {NAN, NAN, NAN};
    const Vector3d zero = Vector3d::Zero();

    if (poly.empty()) {
        return {0.0, nan, nan, false};
    }

    if (bad_normal(n)) {
        if (poly.inside(p)) {
            return {poly.area(), nan, nan, false};
        }
        else {
            return {0.0, nan, nan, false};
        }
    }

    int size = poly.size();

    // Индекс стороны многоугольника, на которой первая точка снаружи,
    // а вторая внутри, sec_out_in -- пересечение
    int idx_out_in = -1;
    Vector3d sec_out_in;

    // Индекс стороны многоугольника, на которой первая точка внутри,
    // а вторая снуружи, sec_in_out -- пересечение
    int idx_in_out = -1;
    Vector3d sec_in_out;

    double dot1 = (poly[size - 1] - p).dot(n);
    for (int i2 = 0; i2 < size; ++i2) {
        const auto &v2 = poly[i2];

        double dot2 = (v2 - p).dot(n);

        // Найдено пересечение
        if (dot1 * dot2 < 0.0) {
            int i1 = (i2 - 1 + size) % poly.size();

            const auto &v1 = poly[i1];
            double Det = (v2 - v1).dot(n);

            Vector3d sec = (dot2 * v1 - dot1 * v2) / Det;
            if (dot1 <= 0.0) {
                idx_in_out = i1;
                sec_in_out = sec;
            } else {
                idx_out_in = i1;
                sec_out_in = sec;
            }
        }

        dot1 = dot2;
    }

    // Нет пересечений, полностью внутри или снаружи
    if (idx_in_out < 0 || idx_out_in < 0) {
        if (dot1 < 0.0) {
            return {poly.area(), Vector3d::Zero(), Vector3d::Zero(), false};
        }
        else {
            return {0.0, zero, zero, false};
        }
    }

    AnS res = {0.0, zero, zero, false};
    if (slice) {
        res.p1 = sec_in_out;
        res.p2 = sec_out_in;
        res.tri_in = (idx_in_out - idx_out_in + size) % size == 1;
        res.tri_out = (idx_out_in - idx_in_out + size) % size == 1;
    }

    res.area += cross(poly.center(), sec_out_in, poly[(idx_out_in + 1) % size]);
    res.area += cross(poly.center(), poly[idx_in_out], sec_in_out);
    res.area += cross(poly.center(), sec_in_out, sec_out_in);

    if (idx_in_out < idx_out_in) {
        idx_in_out += size;
    }

    for (int i = idx_out_in + 1; i < idx_in_out; ++i) {
        const auto& v1 = poly[i % size];
        const auto& v2 = poly[(i + 1) % size];
        res.area += cross(poly.center(), v1, v2);
    }

    res.area *= 0.5;

    return res;
}

double Polygon::clip_area(const Vector3d& p, const Vector3d& n) const {
    return clip_area_and_slice(*this, p, n, false).area;
}

std::tuple<Vector3d, Vector3d> minmax(const Polygon& poly, const Vector3d& n) {
    if (poly.empty()) {
        return std::make_tuple(
                Vector3d{NAN, NAN, NAN},
                Vector3d{NAN, NAN, NAN});
    }

    Vector3d min = poly[0];
    Vector3d max = poly[0];
    double min_dot = min.dot(n);
    double max_dot = min_dot;
    for (int i = 1; i < poly.size(); ++i) {
        double dot = poly[i].dot(n);
        if (dot < min_dot) {
            min = poly[i];
            min_dot = dot;
        } else if (dot > max_dot) {
            max = poly[i];
            max_dot = dot;
        }
    }

    return std::make_tuple(min, max);
}

#define COUNT_ITERATIONS 0

// Найти сечение полигона прямой методом бисекции
Vector3d find_section_bisect(const Polygon& poly, const Vector3d& n, double alpha) {
#if COUNT_ITERATIONS
    static int n_starts = 0; ++n_starts;
    static int total_iterations = 0;
#endif

    auto [v_min, v_max] = minmax(poly, n);

    // Отсекаем сразу только очень близкие к 0.0 и 1.0
    if (alpha < 1.0e-14) {
        return v_min - 0.1 * n;
    } else if (alpha > 1.0 - 1.0e-14) {
        return v_max + 0.1 * n;
    }

    double S = poly.area();

    Vector3d v_avg = 0.5 * (v_min + v_max);

    double f_min = - alpha * S;
    double f_avg = poly.clip_area(v_avg, n) - alpha * S;
    double f_max = (1.0 - alpha) * S;

    int n_iterations = 0;
    double dS = 1.0e-6 * S;
    while (f_max - f_min > dS && n_iterations < 30) {
        if (f_avg < 0.0) {
            v_min = v_avg;
            f_min = f_avg;
        } else {
            v_max = v_avg;
            f_max = f_avg;
        }

        v_avg = 0.5 * (v_min + v_max);
        f_avg = poly.clip_area(v_avg, n) - alpha * S;
        ++n_iterations;
    }
    v_avg = v_min - f_min * (v_max - v_min) / (f_max - f_min);

#if COUNT_ITERATIONS
    total_iterations += n_iterations;
    if (n_starts % 100000 == 0) {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "      avg n_iters (bisect): " << double(total_iterations) / n_starts << "\n";
    }
#endif

    return v_avg;
}

// Сходится дольше, чем метод бисекции, почему?
// Надо ещё метод Брента попробовать
Polygon::section find_section_newton(const Polygon& poly, const Vector3d& n, double alpha) {
#if COUNT_ITERATIONS
    static int n_starts = 0; ++n_starts;
    static int total_iterations = 0;
#endif
    auto [v_min, v_max] = minmax(poly, n);

    // Отсекаем сразу только очень близкие к 0.0 и 1.0
    if (alpha < 1.0e-14) {
        Vector3d out = v_min - 0.1 * n;
        return {out, out};
    } else if (alpha > 1.0 - 1.0e-14) {
        Vector3d out = v_max + 0.1 * n;
        return {out, out};
    }

    double DV = (v_max - v_min).dot(n);
    Vector3d dv = DV * n;

    double S = poly.area();

    double x = alpha;

    AnS func = clip_area_and_slice(poly, v_min + x * dv, n, true);

    double err = std::abs(func.area / S - alpha);
    int n_iterations = 0;
    while (err > 1.0e-12 && n_iterations < 30){
        if (func.tri_in) {
            // Отсекается внутренний треугольник
            x = std::sqrt(2.0 * alpha * x * S / (DV * func.slice()));
        }
        else if (func.tri_out) {
            // Отсекается внешний треугольник
            x = 1.0 - std::sqrt(2.0 * (1.0 - alpha) * (1.0 - x) * S / (DV * func.slice()));
        }
        else {
            x -= (func.area - alpha * S) / (DV * func.slice());
        }
        x = std::max(0.0, std::min(x, 1.0));

        func = clip_area_and_slice(poly, v_min + x * dv, n, true);

        err = std::abs(func.area / S - alpha);

        ++n_iterations;
    }

#if COUNT_ITERATIONS
    total_iterations += n_iterations;
    if (n_starts % 100000 == 0) {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "      avg n_iters (newton): " << double(total_iterations) / n_starts << "\n";
    }
#endif

    if (alpha < 1.0e-14) {
        Vector3d out = v_min - 0.1 * n;
        return {out, out};
    } else if (alpha > 1.0 - 1.0e-14) {
        Vector3d out = v_max + 0.1 * n;
        return {out, out};
    }
    else {
        return {0.5 * (func.p1 + func.p2), func.p2};
    }
}

Polygon::section Polygon::find_section(const Vector3d& _n, double alpha) const {
    const Vector3d nan = {NAN, NAN, NAN};
    if (empty()) {
        return {nan, nan};
    }

    Vector3d n = _n.normalized();

    if (bad_normal(n)) {
        if (alpha < 0.5) {
            return {nan, nan};
        }
        else {
            return {center(), center()};
        }
    }

    //return {find_section_bisect(*this, n, alpha), Vector3d::Zero()};
    return find_section_newton(*this, n, alpha);
}

// Средний угол между векторами p1 - c, p2 - c.
inline double angle_avg(const Vector3d& c, const Vector3d& p1, const Vector3d& p2) {
    return std::atan2(0.5 * (p1.y() + p2.y()) - c.y(),
                      0.5 * (p1.x() + p2.x()) - c.x());
}

// Положительный угол между векторами p1 - c, p2 - c.
inline double angle(const Vector3d& c, const Vector3d& p1, const Vector3d& p2) {
    double d_phi = std::atan2(p2.y() - c.y(), p2.x() - c.x()) -
                   std::atan2(p1.y() - c.y(), p1.x() - c.x());

    if (std::abs(d_phi) >= M_PI) {
        // перескочили через ось
        if (d_phi < 0.0) {
            d_phi += 2.0 * M_PI;
        }
        else {
            d_phi -= 2.0 * M_PI;
        }
    }
    return d_phi;
}

double Polygon::disk_clip_area(const Vector3d& c, double R) const {
    obj::circle circle{c, R};

    double res = 0.0;
    int inside = 0;   // число граней полностью внутри
    int outside = 0;  // число граней полностью снаружи
    for (int i = 0; i < size(); ++i) {
        const auto &v1 = vs[i];
        const auto &v2 = vs[(i + 1) % size()];

        auto sec = intersection2D::find(circle, obj::segment{v1, v2});

        auto inner = [this, &c, &R](const Vector3d& p1, const Vector3d& p2) -> double {
            return cross(m_center, p1, p2);
        };

        auto outer = [this, &c, &R](const Vector3d& p1, const Vector3d& p2) -> double {
            return R * R * angle(c, p1, p2) + cross(p2 - p1, m_center - c);
        };

        double part = 0.0;
        if (sec.exist) {
            if (sec.t1 < 0.0) {
                if (sec.t2 < 0.0) {
                    // t1 < t2 < 0.0
                    part = outer(v1, v2);
                    ++outside;
                }
                else if (sec.t2 < 1.0) {
                    // t1 < 0.0 < t2 < 1.0
                    part = inner(v1, sec.p2) + outer(sec.p2, v2);
                }
                else {
                    // t1 < 0.0 < 1.0 < t2
                    part = inner(v1, v2);
                    ++inside;
                }
            }
            else if (sec.t1 < 1.0) {
                if (sec.t2 < 1.0) {
                    // 0.0 < t1 < t2 < 1.0
                    part = outer(v1, sec.p1) + inner(sec.p1, sec.p2) + outer(sec.p2, v2);
                }
                else {
                    // 0.0 < t1 < 1.0 < t2
                    part = outer(v1, sec.p1) + inner(sec.p1, v2);
                }
            }
            else {
                // 1.0 < t1 < t2
                part = outer(v1, v2);
                ++outside;
            }
        }
        else {
            // Нет пересечений
            part = outer(v1, v2);
            ++outside;
        }

        res += 0.5 * part;
    }

    double S = this->area();
    res = std::max(0.0, std::min(res, S));
    if (outside == this->size()) {
        return 0.0;
    }
    if (inside == this->size()) {
        return S;
    }
    return res;
}

Vector3d Polygon::disk_clip_normal(const Vector3d& c, double R) const {
    obj::circle circle{c, R};

    // Ищем точку снаружи
    int out_idx = 0;
    while (out_idx < size() && (vs[out_idx] - c).norm() < R) {
        ++out_idx;
    }

    // Все точки внутри
    if (out_idx == size()) {
        return Vector3d::Zero();
    }

    // Ищем пересечения
    std::vector<Vector3d> inout;
    for (int i = out_idx; i < out_idx + size(); ++i) {
        obj::segment segment{
            vs[i % size()],
                vs[(i + 1) % size()]
        };

        auto sec = intersection2D::find(circle, segment);

        if (sec.exist) {
            if (0.0 <= sec.t1 && sec.t1 < 1.0) {
                inout.push_back(sec.p1);
            }
            if (0.0 <= sec.t2 && sec.t2 < 1.0) {
                inout.push_back(sec.p2);
            }
        }
    }

    // Нет пересечений
    if (inout.empty()) {
        return Vector3d::Zero();
    }

    inout.insert(inout.begin(), inout.back());
    inout.pop_back();

    if (inout.size() % 2 != 0) {
        //throw std::runtime_error("Something is wrong");
    }

    double phi = 0.0;
    double max_delta = 0.0;
    for (int i = 0; i < size() / 2; i += 2) {
        Vector3d &v1 = inout[2 * i];
        Vector3d &v2 = inout[2 * i + 1];

        double dphi = angle(c, v1, v2);
        bool wide = cross(c, v1, v2) < 0.0;

        if (wide) {
            dphi += 2.0 * M_PI;
        }

        if (dphi > max_delta) {
            max_delta = dphi;
            phi = angle_avg(c, v1, v2);

            if (wide) {
                phi += M_PI;
            }
        }
    }

    return {std::cos(phi), std::sin(phi), 0.0};
}

double Polygon::volume_fraction(
        const std::function<bool(const Vector3d&)>& inside,
        int n_points) const {
    int n = std::max(n_points / size(), 1);
    double res = 0.0;
    for (int i = 0; i < size(); ++i) {
        int j = (i + 1) % size();
        Triangle tri(vs[i], vs[j], m_center);
        res += tri.volume_fraction(inside, n) * tri.area();
    }
    return res / area();
}

double Polygon::integrate_low(const std::function<double(const Vector3d&)>& func, int n) const {
    double res = 0.0;
    for (int i = 0; i < size(); ++i) {
        int j = (i + 1) % size();
        Triangle tri(vs[i], vs[j], m_center);
        res += tri.integrate_low(func, n);
    }
    return res;
}

double Polygon::integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const {
    double res = 0.0;
    for (int i = 0; i < size(); ++i) {
        int j = (i + 1) % size();
        Triangle tri(vs[i], vs[j], m_center);
        res += tri.integrate_mid(func, n);
    }
    return res;
}

double Polygon::integrate_high(const std::function<double(const Vector3d&)>& func, int n) const {
    double res = 0.0;
    for (int i = 0; i < size(); ++i) {
        int j = (i + 1) % size();
        Triangle tri(vs[i], vs[j], m_center);
        res += tri.integrate_high(func, n);
    }
    return res;
}

double Polygon::integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const {
    double res = 0.0;
    for (int i = 0; i < size(); ++i) {
        int j = (i + 1) % size();
        Triangle tri(vs[i], vs[j], m_center);
        res += tri.integrate_extra(func, n);
    }
    return res;
}

std::ostream& operator<<(std::ostream& os, const Polygon& poly) {
    if (poly.empty()) {
        os << "{ }";
        return os;
    }

    os << std::scientific << std::setprecision(6);

    os << "{\n";
    for (int i = 0; i < poly.size() - 1; ++i) {
        os << "    { " << poly[i].x() << ", " << poly[i].y() << ", " << poly[i].z() << " },\n";
    }
    Vector3d last = poly[poly.size() - 1];
    os << "    { " << last.x() << ", " << last.y() << ", " << last.z() << " }\n";
    os << "};";
    return os;
}

} // namespace zephyr::geom