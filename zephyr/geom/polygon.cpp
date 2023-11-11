#include <iomanip>
#include <zephyr/geom/polygon.h>

namespace zephyr::geom {

inline bool bad_normal(const Vector3d& n) {
    return n.hasNaN() || n.squaredNorm() == 0.0;
}

Polygon::Polygon(std::vector<Vector3d>& vertices)
    : vs(vertices.data()), m_size(vertices.size()) { }

Polygon::Polygon(Vector3d* buff, int size)
    : vs(buff), m_size(size) { }

Vector3d Polygon::center() const {
    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < m_size; ++i) {
        C += vs[i];
    }
    return C / m_size;
}

void Polygon::sort(const Vector3d& _c) {
    Vector3d c = _c.hasNaN() ? center() : _c;

    std::sort(vs, vs + m_size,
              [&c](const Vector3d &a, const Vector3d &b) -> bool {
                  double phi_a = std::atan2(a.y() - c.y(), a.x() - c.x());
                  double phi_b = std::atan2(b.y() - c.y(), b.x() - c.x());
                  return phi_a < phi_b;
              });
}

bool Polygon::inside(const Vector3d& p) const {
    const Vector3d n = {0.123, 0.346433, 0.0};

    int counter = 0;
    for (int i = 0; i < m_size; ++i) {
        const auto& v1 = vs[i];
        const auto& v2 = vs[(i + 1) % m_size];

        double dot1 = (v1 - p).dot(n);
        double dot2 = (v2 - p).dot(n);

        // Проверяем пересечение
        if (dot1 * dot2 < 0.0) {
            ++counter;
        }
    }

    return counter % 2 > 0;
}

double Polygon::area(const Vector3d& _c) const {
    Vector3d c = _c.hasNaN() ? center() : _c;

    double S = 0.0;
    for (int i = 0; i < m_size; ++i) {
        const Vector3d &v1 = vs[i];
        const Vector3d &v2 = vs[(i + 1) % m_size];

        Vector3d face_c = 0.5 * (v1 + v2);
        Vector3d normal = {v2.y() - v1.y(), v1.x() - v2.x(), 0.0};

        S += std::abs((face_c - c).dot(normal));
    }

    return 0.5 * S;
}

Vector3d Polygon::centroid(double S) const {
    if (S == 0.0) {
        S = area();
    }

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < m_size; ++i) {
        auto &v1 = vs[i];
        auto &v2 = vs[(i + 1) % m_size];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * S);

    return C;
}

double Polygon::volume_as() const {
    return NAN;
}

Polygon::MinMax Polygon::minmax(const Vector3d& n) const {
    auto comp = [&n](const Vector3d &a, const Vector3d &b) -> bool {
        return a.dot(n) < b.dot(n);
    };
    auto mm = std::minmax_element(vs, vs + m_size, comp);
    return {*mm.first, *mm.second};
}

void Polygon::clip(const Vector3d& p, const Vector3d& n, Polygon& part) const {
    // Подразумевается, что в part выделен массив вершин достаточного размера
    part.m_size = 0;
    if (empty()) {
        return;
    }

    for (int i = 0; i < m_size; ++i) {
        const auto& v1 = vs[i];
        const auto& v2 = vs[(i + 1) % m_size];

        double dot1 = (v1 - p).dot(n);
        double dot2 = (v2 - p).dot(n);
        double Det = (v2 - v1).dot(n);

        double eps = 0.0e-12;

        // Проверяем вершину (если внутри)
        if (dot1 <= eps) {
            part[part.m_size] = v1;
            ++part.m_size;
        }

        // Проверяем пересечение
        if ((dot1 - eps) * (dot2 - eps) < 0.0) {
            part[part.m_size] = (dot2 * v1 - dot1 * v2) / Det;
            ++part.m_size;
        }
    }
}

void Polygon::clip(const Vector3d& p, const Vector3d& n, Polygon& part, Line& line) const {
    // Подразумевается, что в part выделен массив вершин достаточного размера

    part.m_size = 0;
    line[0] = line[1] = Vector3d::Zero();

    if (empty()) {
        return;
    }

    // Число пересечений
    int n_inters = 0;
    for (int i = 0; i < m_size; ++i) {
        const auto& v1 = vs[i];
        const auto& v2 = vs[(i + 1) % m_size];

        double dot1 = (v1 - p).dot(n);
        double dot2 = (v2 - p).dot(n);
        double Det = (v2 - v1).dot(n);

        const double eps = 0.0;

        // Проверяем вершину (если внутри)
        if (dot1 <= eps) {
            part[part.m_size] = v1;
            ++part.m_size;
        }

        // Проверяем пересечение
        if ((dot1 - eps) * (dot2 - eps) < 0.0) {
            Vector3d inter = (dot2 * v1 - dot1 * v2) / Det;
            part[part.m_size] = inter;
            line[n_inters]  = inter;

            ++part.m_size;
            ++n_inters;
        }
    }
}

double Polygon::clip_area(const Vector3d& p, const Vector3d& n, Polygon& part) const {
    if (bad_normal(n)) {
        return inside(p) ? 1.0 : 0.0;
    }

    clip(p, n, part);
    return part.area();
}

Polygon::AnS Polygon::clip_area_and_slice(const Vector3d& p, const Vector3d& n, Polygon& part) const {
    if (bad_normal(n)) {
        return {.area=(inside(p) ? 1.0 : 0.0), .slice=NAN};
    }

    Line line(Vector3d::Zero(), Vector3d::Zero());
    clip(p, n, part, line);
    return {.area=part.area(), .slice=line.length()};
}

Vector3d Polygon::find_section(const Vector3d& n, double alpha, Polygon& part) const {
    if (bad_normal(n)) {
        return alpha < 0.5 ? vs[0] + 10.0 * (vs[1] - vs[0]) : center();
    }

    return find_section_bisect(n, alpha, part);
    //return find_section_newton(n, alpha, part);
}

#define COUNT_ITERATIONS 0

Vector3d Polygon::find_section_bisect(const Vector3d& n, double alpha, Polygon& part) const {
    // Допускается оптимизация, использовать метод Ньютона
    Vector3d C = center();
    double S = area();

    auto mm = minmax(n);
    Vector3d v_min = mm.min;
    Vector3d v_max = mm.max;

    // Отсекаем сразу только очень близкие к 0.0 и 1.0
    if (alpha < 1.0e-14) {
        return v_min - 0.1 * n;
    } else if (alpha > 1.0 - 1.0e-14) {
        return v_max + 0.1 * n;
    }

    Vector3d v_avg = 0.5 * (v_min + v_max);

    double f_min = - alpha * S;
    double f_avg = clip_area(v_avg, n, part) - alpha * S;
    double f_max = (1.0 - alpha) * S;

#if COUNT_ITERATIONS
    static int n_starts = 0; ++n_starts;
    static int total_iterations = 0;
#endif

    int n_iterations = 0;
    while (f_max - f_min > 1.0e-6 * S && n_iterations < 30) {
        if (f_avg < 0.0) {
            v_min = v_avg;
            f_min = f_avg;
        } else {
            v_max = v_avg;
            f_max = f_avg;
        }

        v_avg = 0.5 * (v_min + v_max);
        f_avg = clip_area(v_avg, n, part) - alpha * S;
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
Vector3d Polygon::find_section_newton(const Vector3d& n, double alpha, Polygon& part) const {
    // Допускается оптимизация, использовать метод Ньютона
    Vector3d C = center();
    double S = area();

    auto mm = minmax(n);
    Vector3d v_min = mm.min;
    Vector3d v_max = mm.max;

    double DV = (v_max - v_min).dot(n);
    Vector3d dv = DV * n;

    // Отсекаем сразу только очень близкие к 0.0 и 1.0
    if (alpha < 1.0e-14) {
        return v_min - 0.1 * n;
    } else if (alpha > 1.0 - 1.0e-14) {
        return v_max + 0.1 * n;
    }

#if COUNT_ITERATIONS
    static int n_starts = 0; ++n_starts;
    static int total_iterations = 0;
#endif

    double x = alpha;

    AnS    res = clip_area_and_slice(v_min + x * dv, n, part);
    double err = std::abs(res.area / S - 1.0);

    int n_iterations = 0;
    while (err > 1.0e-6 && n_iterations < 30) {
        x -= (res.area - alpha * S) / (DV * res.slice);
        res = clip_area_and_slice(v_min + x * dv, n, part);

        err = std::abs(clip_area(v_min + x * dv, n, part) / S - 1.0);
        ++n_iterations;
    }

#if COUNT_ITERATIONS
    total_iterations += n_iterations;
    if (n_starts % 100000 == 0) {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "      avg n_iters (newton): " << double(total_iterations) / n_starts << "\n";
    }
#endif

    if (x < 1.0e-14) {
        return v_min - 0.1 * n;
    } else if (x > 1.0 - 1.0e-14) {
        return v_max + 0.1 * n;
    } else {
        return v_min + x * dv;
    }
}

} // namespace zephyr::geom