/// @file Тестирование класса Polygon. Проверка функций интегрирования
/// по полигону, отсечения линией и кругом.

#include <iostream>

#include <zephyr/geom/polygon.h>
#include <zephyr/geom/triangle.h>

#include <zephyr/utils/matplotlib.h>

using namespace zephyr::geom;
namespace plt = zephyr::utils::matplotlib;

inline double cross(const Vector3d& v1, const Vector3d& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

// Пересечение двух прямых, параметризованных следующим образом:
// a(t) = a1 + tau1 * t
// b(s) = a2 + tau2 * s
// Возвращает два параметра (t, s) пересечения
// a(t) = b(s)
std::array<double, 2> intersection(
        const Vector3d& a1, const Vector3d& tau1,
        const Vector3d& a2, const Vector3d& tau2) {

    double det = cross(tau1, tau2);
    if (std::abs(det) < 1.0e-12) {
        return {NAN, NAN};
    }

    Vector3d b = a2 - a1;

    double t = cross(b, tau2) / det;
    double s = cross(b, tau1) / det;

    return {t, s};
}

// Пересечение двух прямых, параметризованных следующим образом:
// a(t) = a1 + tau1 * t
// b(s) = a2 + tau2 * s
// Возвращает два параметра (t, s) пересечения
// a(t) = b(s)
std::array<double, 2> intersection2(
        const Vector3d& a1, const Vector3d& tau1,
        const Vector3d& a2, const Vector3d& tau2, double w) {

    double axx = (+1.0 + 0.25 * w * w) * tau1.dot(tau1);
    double ayy = (+1.0 + 0.25 * w * w) * tau2.dot(tau2);
    double axy = (-1.0 + 0.25 * w * w) * tau1.dot(tau2);
    double bx = (a2 - a1).dot(tau1);
    double by = (a1 - a2).dot(tau2);

    double det = axx * ayy - axy * axy;

    double t = (bx * ayy - by * axy) / det;
    double s = (by * axx - bx * axy) / det;

    return {t, s};
}


struct FaceClip {
    double alpha;  // [0, 1]
    double sgn;    // {-1, +1},
                   // +1 от первой точки,
                   // -1 - от второй точки
};

// Отсечение полуплоскостью (p, n) от грани (f1, f2)
FaceClip clip_face(
        const Vector3d& p, const Vector3d& n,
        const Vector3d& f1, const Vector3d& f2) {
    Vector3d tau = f2 - f1;

    double det = tau.dot(n);
    if (std::abs(det) < 1.0e-12) {
        if (cross(tau, n) > 0.0) {
            return {1.0, 1.0};
        }
        else {
            return {0.0, 1.0};
        }
    }

    double t = n.dot(p - f1) / det;
    t = std::max(0.0, std::min(t, 1.0));

    if (tau.dot(n) > 0.0) {
        return {t, 1.0};
    }
    else {
        return {1.0 - t, -1.0};
    }
}

enum show {
    LESS,
    MORE,
    ALL
};

struct BezierSquare {
    Vector3d p0;
    Vector3d p1;
    Vector3d p2;

    Vector3d get(double t) const {
        return (1 - t) * (1 - t) * p0 + 2 * t * (1 - t) * p1 + t * t * p2;
    }

    std::vector<double> xs(int n, show s = ALL) const {
        std::vector<double> x;
        x.reserve(n);
        for (int i = 0; i < n; ++i) {
            auto p = get(i / (n - 1.0));
            bool add = s == show::ALL ||
                       (s == show::LESS && p.x() <= 0.0) ||
                       (s == show::MORE && p.x() >= 0.0);
            if (add) {
                x.push_back(p.x());
            }
        }
        return x;
    }

    std::vector<double> ys(int n, show s = ALL) const {
        std::vector<double> y;
        y.reserve(n);
        for (int i = 0; i < n; ++i) {
            auto p = get(i / (n - 1.0));

            bool add = s == show::ALL ||
                    (s == show::LESS && p.x() <= 0.0) ||
                    (s == show::MORE && p.x() >= 0.0);

            if (add) {
                y.push_back(p.y());
            }
        }
        return y;
    }
};

struct BezierCubic {
    Vector3d p0;
    Vector3d p1;
    Vector3d p2;
    Vector3d p3;

    Vector3d get(double t) const {
        return (1 - t) * (1 - t) * (1 - t) * p0 + 3 * t * (1 - t) * (1 - t) * p1 +
               3 * t * t * (1 - t) * p2 + t * t * t * p3;
    }

    std::vector<double> xs(int n) const {
        std::vector<double> x(n);
        for (int i = 0; i < n; ++i) {
            x[i] = get(i / (n - 1.0)).x();
        }
        return x;
    }

    std::vector<double> ys(int n) const {
        std::vector<double> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = get(i / (n - 1.0)).y();
        }
        return y;
    }
};

inline double sgn(double x) {
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

inline double some_weight_i(double x) {
    return x < 0.0 ? 0.0 : 3.0 * x * x / (1.0 + 2.0 * x * x * x);
}

inline double some_weight_f(double x) {
    return x < 0.0 ? 0.0 : 1.0 / (1.0 + x);
}

void test1(double phi1, double phi2, double alpha1, double alpha2) {
    Vector3d v1{-1.0, 0.0, 0.0};
    Vector3d v2{0.0, 0.0, 0.0};
    Vector3d v3{1.0, 0.0, 0.0};
    Vector3d v4{-1.0, 1.0, 0.0};
    Vector3d v5{0.0, 1.0, 0.0};
    Vector3d v6{1.0, 1.0, 0.0};

    Polygon cell_1 = {v1, v2, v5, v4};
    Polygon cell_2 = {v2, v3, v6, v5};

    Vector3d f1 = v2;
    Vector3d f2 = v5;
    Vector3d fc = 0.5 * (f1 + f2);

    Vector3d n1 = {std::cos(phi1), std::sin(phi1), 0.0};
    Vector3d n2 = {std::cos(phi2), std::sin(phi2), 0.0};

    auto[a1, b1] = cell_1.find_section(n1, alpha1);
    auto[a2, b2] = cell_2.find_section(n2, alpha2);

    Polygon clip1 = cell_1.clip(a1, n1);
    Polygon clip2 = cell_2.clip(a2, n2);


    plt::figure_size(14.0, 8.0);
    plt::title("Реконструкция по кривой Безье");
    plt::set_aspect_equal();

    plt::plot(cell_1.xs(), cell_2.ys(), "k-");
    plt::plot(cell_2.xs(), cell_2.ys(), "k-");

    plt::fill(clip1.xs(), clip1.ys(), {{"color", "#0000ff0f"}});
    plt::fill(clip2.xs(), clip2.ys(), {{"color", "#0000ff0f"}});

    plt::plot(
            std::vector<double>({2.0 * a1.x() - b1.x(), b1.x()}),
            std::vector<double>({2.0 * a1.y() - b1.y(), b1.y()}), "k--");
    plt::plot(
            std::vector<double>({2.0 * a2.x() - b2.x(), b2.x()}),
            std::vector<double>({2.0 * a2.y() - b2.y(), b2.y()}), "k--");

    auto[t1, t2] = intersection2(a1, b1 - a1, a2, b2 - a2, 0.001);

    Vector3d int1 = a1 + (b1 - a1) * t1;
    Vector3d int2 = a2 + (b2 - a2) * t2;

    plt::plot({int1.x()}, {int1.y()}, "bo");
    plt::plot({int2.x()}, {int2.y()}, "bo");

    auto[s1, q1] = intersection2(a1, b1 - a1, fc, 0.5 * (f2 - f1), 0.001);
    auto[s2, q2] = intersection2(a2, b2 - a2, fc, 0.5 * (f2 - f1), 0.001);

    Vector3d fint1 = fc + 0.5 * (f2 - f1) * q1;
    Vector3d fint2 = fc + 0.5 * (f2 - f1) * q2;

    {
        // Дотянуть пересечения в область около грани
        fint1 = fc + (fint1 - fc) * std::min(1.0, 3.0 * (f2 - f1).norm() / (fint1 - f1).norm());
        fint2 = fc + (fint2 - fc) * std::min(1.0, 3.0 * (f2 - f1).norm() / (fint2 - f1).norm());
    }

    plt::plot({fint1.x()}, {fint1.y()}, "go");
    plt::plot({fint2.x()}, {fint2.y()}, "go");

    if (t1 * t2 <= 0.0 && s1 * s2 <= 0.0) {
        plt::title("Хорошее пересечение интерфейсов и грани");
    } else if (t1 * t2 <= 0.0 && s1 * s2 > 0.0) {
        plt::title("Хорошее пересечение интерфейсов, плохое грани");
    } else if (t1 * t2 > 0.0 && s1 * s2 <= 0.0) {
        plt::title("Хорошее пересечение грани, плохое интерфейса");
    } else if (t1 * t2 > 0.0 && s1 * s2 > 0.0) {
        plt::title("Плохое пересечение всего");
    } else {
        plt::title("Непонятно что");
    }

    BezierSquare bezier_i1;
    BezierSquare bezier_i2;

    BezierSquare bezier_f1;
    BezierSquare bezier_f2;

    BezierSquare bezier_m1;
    BezierSquare bezier_m2;

    Vector3d p1, ii1, if1, if2, ii2, p2;
    p1 = a1;
    ii1 = int1;
    ii2 = int2;
    if1 = fint1;
    if2 = fint2;
    p2 = a2;

    bezier_i1.p0 = p1;
    bezier_i1.p1 = ii1;
    bezier_i1.p2 = p2;

    bezier_i2.p0 = p1;
    bezier_i2.p1 = ii2;
    bezier_i2.p2 = p2;

    bezier_f1.p0 = p1;
    bezier_f1.p1 = if1;
    bezier_f1.p2 = p2;

    bezier_f2.p0 = p1;
    bezier_f2.p1 = if2;
    bezier_f2.p2 = p2;

    // best is 1.0, more is OK (?)
    double Ci = some_weight_i(-t1 * t2);
    double Cf = some_weight_f(sgn(-s1 * s2) * std::abs(q1 * q2));
    double C = Ci / (Ci + Cf);
    std::cout << Ci << " " << Cf << " " << C << "\n";

    bezier_m1.p0 = C * bezier_i1.p0 + (1.0 - C) * bezier_f1.p0;
    bezier_m1.p1 = C * bezier_i1.p1 + (1.0 - C) * bezier_f1.p1;
    bezier_m1.p2 = C * bezier_i1.p2 + (1.0 - C) * bezier_f1.p2;

    bezier_m2.p0 = C * bezier_i2.p0 + (1.0 - C) * bezier_f2.p0;
    bezier_m2.p1 = C * bezier_i2.p1 + (1.0 - C) * bezier_f2.p1;
    bezier_m2.p2 = C * bezier_i2.p2 + (1.0 - C) * bezier_f2.p2;


    plt::plot(bezier_i1.xs(300, LESS), bezier_i1.ys(300, LESS), {{"color",     "blue"},
                                                                 {"linestyle", "dotted"},
                                                                 {"label",     "int square"}});
    plt::plot(bezier_i2.xs(300, MORE), bezier_i2.ys(300, MORE), {{"color",     "blue"},
                                                                 {"linestyle", "dotted"}});

    plt::plot(bezier_f1.xs(300, LESS), bezier_f1.ys(300, LESS), {{"color",     "green"},
                                                                 {"linestyle", "dotted"},
                                                                 {"label",     "face square"}});
    plt::plot(bezier_f2.xs(300, MORE), bezier_f2.ys(300, MORE), {{"color",     "green"},
                                                                 {"linestyle", "dotted"}});

    plt::plot(bezier_m1.xs(300, LESS), bezier_m1.ys(300, LESS), {{"color", "black"},
                                                                 {"label", "weighted"}});
    plt::plot(bezier_m2.xs(300, MORE), bezier_m2.ys(300, MORE), {{"color", "black"}});

    plt::text(-0.9, 0.9, "t1 * t2: " + std::to_string(t1 * t2));
    plt::text(-0.9, 0.7, "weight: " + std::to_string(C));


    /*
    BezierCubic bezier_i3;
    bezier_i3.p0 = p1;
    bezier_i3.p1 = 0.5*(ii1 + ii2);
    bezier_i3.p2 = 0.5*(ii1 + ii2);
    bezier_i3.p3 = p2;

    BezierCubic bezier_f3;
    bezier_f3.p0 = p1;
    bezier_f3.p1 = if1;
    bezier_f3.p2 = if2;
    bezier_f3.p3 = p2;

    double C1 = 0.5;
    BezierCubic bezier_if;
    bezier_if.p0 = C1 * bezier_i3.p0 + (1.0 - C1) * bezier_f3.p0;
    bezier_if.p1 = C1 * bezier_i3.p1 + (1.0 - C1) * bezier_f3.p1;
    bezier_if.p2 = C1 * bezier_i3.p2 + (1.0 - C1) * bezier_f3.p2;
    bezier_if.p3 = C1 * bezier_i3.p3 + (1.0 - C1) * bezier_f3.p3;

    //plt::plot(bezier_i3.xs(100), bezier_i3.ys(100), {{"color",     "blue"},
    //                                                 {"label",     "int cubic"},
    //                                                 {"linestyle", "solid"}});

    //plt::plot(bezier_f3.xs(100), bezier_f3.ys(100), {{"color",     "green"},
    //                                                 {"label",     "face cubic"},
    //                                                 {"linestyle", "solid"}});

    //plt::plot(bezier_if.xs(100), bezier_if.ys(100), {{"color",     "black"},
    //                                                 {"label",     "face cubic"},
    //                                                 {"linestyle", "dashed"}});
    */

    plt::legend();

    plt::tight_layout();
    plt::show();
}

void test2(double phi1, double phi2, double alpha1, double alpha2) {
    Vector3d v1{-1.0, 0.0, 0.0};
    Vector3d v2{0.0, 0.0, 0.0};
    Vector3d v3{1.0, 0.0, 0.0};
    Vector3d v4{-1.0, 1.0, 0.0};
    Vector3d v5{0.0, 1.0, 0.0};
    Vector3d v6{1.0, 1.0, 0.0};

    Polygon cell_1 = {v1, v2, v5, v4};
    Polygon cell_2 = {v2, v3, v6, v5};

    Vector3d f1 = v2;
    Vector3d f2 = v5;
    Vector3d fc = 0.5 * (f1 + f2);

    Vector3d n1 = {std::cos(phi1), std::sin(phi1), 0.0};
    Vector3d n2 = {std::cos(phi2), std::sin(phi2), 0.0};

    if (n1.hasNaN()) {
        n1 = Vector3d::Zero();
    }
    if (n2.hasNaN()) {
        n2 = Vector3d::Zero();
    }


    auto[a1, b1] = cell_1.find_section(n1, alpha1);
    auto[a2, b2] = cell_2.find_section(n2, alpha2);

    Vector3d tau1 = b1 - a1;
    Vector3d tau2 = b2 - a2;

    Polygon clip1 = cell_1.clip(a1, n1);
    Polygon clip2 = cell_2.clip(a2, n2);


    plt::figure_size(14.0, 8.0);
    plt::title("Реконструкция CRP");
    plt::set_aspect_equal();

    plt::plot(cell_1.xs(), cell_2.ys(), "k-");
    plt::plot(cell_2.xs(), cell_2.ys(), "k-");

    plt::fill(clip1.xs(), clip1.ys(), {{"color", "#0000ff0f"}});
    plt::fill(clip2.xs(), clip2.ys(), {{"color", "#0000ff0f"}});

    plt::plot(
            std::vector<double>({2.0 * a1.x() - b1.x(), b1.x()}),
            std::vector<double>({2.0 * a1.y() - b1.y(), b1.y()}), "k--");
    plt::plot(
            std::vector<double>({2.0 * a2.x() - b2.x(), b2.x()}),
            std::vector<double>({2.0 * a2.y() - b2.y(), b2.y()}), "k--");

    auto[as1, dir1] = clip_face(a1, n1, f1, f2);
    auto[as2, dir2] = clip_face(a2, n2, f1, f2);

    Vector3d fint1 = f1 + (f2 - f1) * (dir1 > 0.0 ? as1 : 1.0 - as1);
    Vector3d fint2 = f1 + (f2 - f1) * (dir2 > 0.0 ? as2 : 1.0 - as2);

    plt::plot({fint1.x()}, {fint1.y()}, {{"label", "left"}, {"color", "green"}, {"marker", "x"}});
    plt::plot({fint2.x()}, {fint2.y()}, {{"label", "right"}, {"color", "blue"}, {"marker", "x"}});

    std::cout << "as1, as2: " << as1 << " " << as2 << "\n";

    std::cout << "s1, s2: " << dir1 << " " << dir2 << "\n";

    double inter;
    double a_sig1, a_sig2;
    if (dir1 * dir2 > 0.0) {
        inter = std::min(as1, as2);
    }
    else {
        inter = std::max(0.0, as1 + as2 - 1.0);
    }
    a_sig1 = 0.5 * (as1 + inter);
    a_sig2 = 0.5 * (as2 + inter);
    std::cout << "inter: " << inter << "\n";

    Vector3d aint1 = f1 + (f2 - f1) * (dir1 > 0.0 ? a_sig1 : 1.0 - a_sig1);
    Vector3d aint2 = f1 + (f2 - f1) * (dir2 > 0.0 ? a_sig2 : 1.0 - a_sig2);

    plt::plot({aint1.x()}, {aint1.y()}, {{"color", "green"}, {"marker", "o"}});
    plt::plot({aint2.x()}, {aint2.y()}, { {"color", "blue"}, {"marker", "o"}});
    plt::legend();

    plt::tight_layout();
    plt::show();
}

int main() {
    test1(0.01 * M_PI, 0.51 * M_PI, 0.12, 0.1);
    //test2(0.0, NAN, 0.22, 0.0);

    return 0;
}