#include <zephyr/geom/intersection.h>

namespace zephyr::geom::intersection2D {

inline double sqr(double x) {
    return x * x;
}

// vector product [v1, v2].z
inline double cross2D(const Vector3d& v1, const Vector3d& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

bool exist(const obj::ray& ray, const obj::segment& seg) {
    double dot1 = cross2D(seg.v1 - ray.p, ray.tau);
    double dot2 = cross2D(seg.v2 - ray.p, ray.tau);

    // Проверяем пересечение
    // Пересечение во второй точке не учитывается
    if (dot1 * dot2 <= 0.0 && dot2 != 0.0) {
        // det != 0.0, т.к. dot1 и dot2 разных знаков
        // и вместе не обращаются в ноль
        double det = dot2 - dot1;

        // Пересечение
        Vector3d vi = (dot2 * seg.v1 - dot1 * seg.v2) / det;
        if ((vi - ray.p).dot(ray.tau) > 0.0) {
            return true;
        }
    }

    return false;
}

bool exist(const obj::plane& plane, const obj::segment& seg) {
    double dot1 = plane.n.dot(seg.v1 - plane.p);
    double dot2 = plane.n.dot(seg.v2 - plane.p);

    return dot1 * dot2 <= 0.0;
}


// Пересечение двух прямых, параметризованных следующим образом:
// a(t) = a1 + tau1 * t
// b(s) = a2 + tau2 * s
// Возвращает два параметра (t, s) пересечения
// a(t) = b(s)
Vector3d find_fast(const obj::line& line1, const obj::line& line2) {

    double det = cross2D(line1.tau, line2.tau);
    Vector3d b = line2.p - line1.p;

    // Параметр на первой линии
    double t = cross2D(b, line2.tau) / det;
    return line1.get(t);
}

Vector3d find_fast(const obj::segment& seg1, const obj::segment& seg2) {
    return find_fast(obj::line{seg1.v1, seg1.v2 - seg1.v1},
                     obj::line{seg2.v2, seg2.v2 - seg2.v1});
}

Vector3d find_fast(const obj::line& line, const obj::segment& seg) {
    return find_fast(line, obj::line{seg.v2, seg.v2 - seg.v1});
}

circle_segment_intersection find(
        const obj::circle &circle, const obj::line &line) {

    circle_segment_intersection res;

    double A = line.tau.squaredNorm();
    double B = line.tau.dot(line.p - circle.c);
    double C = (line.p - circle.c).squaredNorm() - sqr(circle.r);

    double D = B * B - A * C;

    if (D <= 0.0) {
        res.exist = false;
    }
    else {
        res.exist = true;

        double SD = std::sqrt(D);
        res.t1 = (-B - SD) / A;
        res.t2 = (-B + SD) / A;

        res.p1 = line.get(res.t1);
        res.p2 = line.get(res.t2);
    }

    return res;
}

circle_segment_intersection find(
        const obj::circle &circle, const obj::segment &seg) {
    return intersection2D::find(circle, obj::line{seg.v1, seg.tau()});
}


}