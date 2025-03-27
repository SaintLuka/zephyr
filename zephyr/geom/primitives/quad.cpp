#include <cstring>
#include <iostream>

#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/quad.h>

namespace zephyr::geom {

inline double sqr(double x) {
    return x * x;
}

// ============================================================================
//                                    QUAD
// ============================================================================

Quad::Quad(
        const Vector3d &v00,
        const Vector3d &v01,
        const Vector3d &v10,
        const Vector3d &v11)
        : verts({v00, v01, v10, v11}) {}

Vector3d Quad::operator()(double x, double y) const {
    return get(x, y);
}

Vector3d Quad::get(double x, double y) const {
    return 0.25 * (((1.0 - x) * (1.0 - y)) * vs<-1, -1>() +
                   ((1.0 + x) * (1.0 - y)) * vs<+1, -1>() +
                   ((1.0 - x) * (1.0 + y)) * vs<-1, +1>() +
                   ((1.0 + x) * (1.0 + y)) * vs<+1, +1>());
}

Vector3d Quad::tangent_x(double x, double y) const {
    return 0.25 * ((1.0 - y) * (vs<+1, -1>() - vs<-1, -1>()) +
                   (1.0 + y) * (vs<+1, +1>() - vs<-1, +1>()));
}

Vector3d Quad::tangent_y(double x, double y) const {
    return 0.25 * ((1.0 - x) * (vs<-1, +1>() - vs<-1, -1>()) +
                   (1.0 + x) * (vs<+1, +1>() - vs<+1, -1>()));
}

Vector3d Quad::normal(double x, double y, const Vector3d &c) const {
    Vector3d n = tangent_x(x, y).cross(tangent_y(x, y)).normalized();
    return (verts[0] - c).dot(n) > 0 ? n : -n;
}

double Quad::Jacobian(double x, double y) const {
    return tangent_x(x, y).cross(tangent_y(x, y)).norm();
}

Vector3d Quad::center() const {
    return 0.25 * (verts[0] + verts[1] + verts[2] + verts[3]);
}

Vector3d Quad::normal(const Vector3d &c) const {
    Vector3d n = (verts[3] - verts[0]).cross(verts[2] - verts[1]);
    n.normalize();
    return (verts[0] - c).dot(n) > 0.0 ? n : -n;
}

double Quad::area() const {
    return 0.5 * (verts[3] - verts[0]).cross(verts[2] - verts[1]).norm();
}

Vector3d Quad::area_n(const Vector3d& c) const {
    return Quad::area() * Quad::normal(c);
}

double Quad::volume_as() const {
    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    double V = 0.0;
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = verts[ord[i]];
        auto &v2 = verts[ord[(i + 1) % 4]];

        V -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }

    return V / 6.0;
}

Vector3d Quad::centroid(double area) const {
    if (area == 0.0) {
        area = Quad::area();
    }

    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = verts[ord[i]];
        auto &v2 = verts[ord[(i + 1) % 4]];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * area);

    return C;
}

Vector3d Quad::centroid_as(double vol_as) const {
    if (vol_as == 0.0) {
        vol_as = Quad::volume_as();
    }

    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = verts[ord[i]];
        auto &v2 = verts[ord[(i + 1) % 4]];

        C.x() += (+1.5 * std::pow(v1.y(), 2) + v1.y() * v2.y() + 0.5 * std::pow(v2.y(), 2)) * std::pow(v1.x(), 2) +
                 (-1.5 * std::pow(v2.y(), 2) - v1.y() * v2.y() - 0.5 * std::pow(v1.y(), 2)) * std::pow(v2.x(), 2) +
                 (v2.y() * v2.y() - v1.y() * v1.y()) * v1.x() * v2.x();

        C.y() -= (std::pow(v1.y(), 3) + std::pow(v1.y(), 2) * v2.y() +
                  std::pow(v2.y(), 3) + std::pow(v2.y(), 2) * v1.y()) * (v2.x() - v1.x());
    }
    C /= (12.0 * vol_as);

    return C;
}

// Реализация интегралов по Quad и SqQuad
namespace integral2D {

struct Node {
    double i, j, w;
};

// type Q is Quad or SqQuad
// N -- число ячеек, точность определения объемной доли ~ 1/N
template<typename Q>
double volume_fraction(const Q& quad, const std::function<bool(const Vector3d&)>& inside, int N) {
    int n = std::round(std::sqrt(N));
    double h = 1.0 / n;
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = (2 * i + 1) * h - 1.0;
        for (int j = 0; j < n; ++j) {
            double y = (2 * j + 1) * h - 1.0;
            if (inside(quad.get(x, y))) {
                res += quad.Jacobian(x, y);
            }
        }
    }
    return res * sqr(2.0 * h) / quad.area();
}

// type Q is Quad or SqQuad
template<typename Q>
double integrate_low(const Q &map2D, const std::function<double(const Vector3d &)> &func, int n) {
    double h = 1.0 / n;
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = (2 * i + 1) * h - 1.0;
        for (int j = 0; j < n; ++j) {
            double y = (2 * j + 1) * h - 1.0;
            res += func(map2D.get(x, y)) * map2D.Jacobian(x, y);
        }
    }
    return res * sqr(2.0 * h);
}

// type Q is Quad or SqQuad
template<typename Q>
double integrate_mid(const Q &map2D, const std::function<double(const Vector3d &)> &func, int n) {
    static const double cm = 1.0 - 1.0 / std::sqrt(3.0);
    static const double cp = 1.0 + 1.0 / std::sqrt(3.0);

    double h = 1.0 / n;
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x1 = (2 * i + cm) * h - 1.0;
        double x2 = (2 * i + cp) * h - 1.0;

        for (int j = 0; j < n; ++j) {
            double y1 = (2 * j + cm) * h - 1.0;
            double y2 = (2 * j + cp) * h - 1.0;

            res += func(map2D.get(x1, y1)) * map2D.Jacobian(x1, y1) +
                   func(map2D.get(x1, y2)) * map2D.Jacobian(x1, y2) +
                   func(map2D.get(x2, y1)) * map2D.Jacobian(x2, y1) +
                   func(map2D.get(x2, y2)) * map2D.Jacobian(x2, y2);
        }
    }
    return res * sqr(h);
}

// type Q is Quad or SqQuad
template<typename Q>
double integrate_high(const Q &map2D, const std::function<double(const Vector3d &)> &func, int n) {
    static const double a = std::sqrt((114.0 - 3.0 * std::sqrt(583.0)) / 287.0);
    static const double b = std::sqrt((114.0 + 3.0 * std::sqrt(583.0)) / 287.0);
    static const double c = std::sqrt(6.0 / 7.0);

    static const double wa = 307.0 / 810.0 + 923.0 / 270.0 / std::sqrt(583.0);
    static const double wb = 307.0 / 810.0 - 923.0 / 270.0 / std::sqrt(583.0);
    static const double wc = 98.0 / 405.0;

    const double h = 1.0 / n;
    const double ah = a * h;
    const double bh = b * h;
    const double ch = c * h;

    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x_c = (2 * i + 1) * h - 1.0;
        double x_ma = x_c - ah;
        double x_pa = x_c + ah;
        double x_mb = x_c - bh;
        double x_pb = x_c + bh;
        double x_mc = x_c - ch;
        double x_pc = x_c + ch;

        for (int j = 0; j < n; ++j) {
            double y_c = (2 * j + 1) * h - 1.0;
            double y_ma = y_c - ah;
            double y_pa = y_c + ah;
            double y_mb = y_c - bh;
            double y_pb = y_c + bh;
            double y_mc = y_c - ch;
            double y_pc = y_c + ch;

            res += wa * (func(map2D.get(x_ma, y_ma)) * map2D.Jacobian(x_ma, y_ma) +
                         func(map2D.get(x_ma, y_pa)) * map2D.Jacobian(x_ma, y_pa) +
                         func(map2D.get(x_pa, y_ma)) * map2D.Jacobian(x_pa, y_ma) +
                         func(map2D.get(x_pa, y_pa)) * map2D.Jacobian(x_pa, y_pa)) +
                   wb * (func(map2D.get(x_mb, y_mb)) * map2D.Jacobian(x_mb, y_mb) +
                         func(map2D.get(x_mb, y_pb)) * map2D.Jacobian(x_mb, y_pb) +
                         func(map2D.get(x_pb, y_mb)) * map2D.Jacobian(x_pb, y_mb) +
                         func(map2D.get(x_pb, y_pb)) * map2D.Jacobian(x_pb, y_pb)) +
                   wc * (func(map2D.get(x_mc, y_c)) * map2D.Jacobian(x_mc, y_c) +
                         func(map2D.get(x_pc, y_c)) * map2D.Jacobian(x_pc, y_c) +
                         func(map2D.get(x_c, y_mc)) * map2D.Jacobian(x_c, y_mc) +
                         func(map2D.get(x_c, y_pc)) * map2D.Jacobian(x_c, y_pc));
        }
    }
    return res * sqr(h);
}

// type Q is Quad or SqQuad
template<typename Q>
double integrate_extra(const Q &map2D, const std::function<double(const Vector3d &)> &func, int n) {
    static const std::array<Node, 37> weights = {
            Node{-9.1324172678365734e-01, -9.3149501298985338e-01, 9.7554807790077372e-03},
            {-9.7865852404709552e-01, -6.1165306951584819e-01, 1.0171445385633233e-02},
            {-9.3749714215853275e-01, -1.0520434537393597e-02, 2.1998021656947090e-02},
            {-9.6529260227686897e-01, 5.5065509547510505e-01, 1.3219343426915007e-02},
            {-9.2613868702909463e-01, 9.1983989754282880e-01, 9.7076558625866834e-03},
            {-5.1568654928683255e-01, -9.6524350839166417e-01, 1.3695752271564185e-02},
            {-7.1978981718941282e-01, -7.5183138235659153e-01, 2.8074929138724544e-02},
            {-8.0477746117089521e-01, -3.6511753726073226e-01, 3.0869749835561042e-02},
            {-7.4287042853522023e-01, 3.1819306172871120e-01, 3.9150952078167579e-02},
            {-7.3881591396021928e-01, 7.3892523474529592e-01, 2.6750632551567975e-02},
            {-5.8145986766020930e-01, 9.7609319698904362e-01, 1.0839251986837972e-02},
            {-2.8979445545058968e-01, -7.2372628537340544e-01, 3.9420078420011065e-02},
            {-4.8747307818939967e-01, -3.4316793927577516e-01, 3.7132288715576361e-02},
            {-4.8059467301081205e-01, 2.8020267677347777e-02, 4.4362000373970513e-02},
            {-3.1422994116334935e-01, 4.6772084934377323e-01, 4.2824895303643057e-02},
            {-3.6434592892724771e-01, 7.7527412891380232e-01, 2.9831598487120591e-02},
            {2.9287160126081346e-02, -9.2759405462351985e-01, 2.4102633262065763e-02},
            {-6.6026533930988895e-02, -4.4048747060251170e-01, 3.8331764055554537e-02},
            {2.2204460492503131e-16, 4.4408920985006262e-16, 5.9523052817090159e-02},
            {6.6026533930988229e-02, 4.4048747060251259e-01, 3.8331764055554475e-02},
            {-2.9287160126081568e-02, 9.2759405462351996e-01, 2.4102633262065739e-02},
            {3.6434592892724771e-01, -7.7527412891380221e-01, 2.9831598487120560e-02},
            {3.1422994116334890e-01, -4.6772084934377356e-01, 4.2824895303643078e-02},
            {4.8059467301081260e-01, -2.8020267677348443e-02, 4.4362000373970506e-02},
            {4.8747307818939878e-01, 3.4316793927577538e-01, 3.7132288715576459e-02},
            {2.8979445545058957e-01, 7.2372628537340589e-01, 3.9420078420011016e-02},
            {5.8145986766020896e-01, -9.7609319698904351e-01, 1.0839251986837982e-02},
            {7.3881591396021906e-01, -7.3892523474529614e-01, 2.6750632551567958e-02},
            {7.4287042853522034e-01, -3.1819306172871187e-01, 3.9150952078167538e-02},
            {8.0477746117089510e-01, 3.6511753726073182e-01, 3.0869749835561091e-02},
            {7.1978981718941282e-01, 7.5183138235659142e-01, 2.8074929138724562e-02},
            {5.1568654928683255e-01, 9.6524350839166417e-01, 1.3695752271564183e-02},
            {9.2613868702909441e-01, -9.1983989754282891e-01, 9.7076558625866817e-03},
            {9.6529260227686908e-01, -5.5065509547510527e-01, 1.3219343426914993e-02},
            {9.3749714215853275e-01, 1.0520434537393264e-02, 2.1998021656947069e-02},
            {9.7865852404709530e-01, 6.1165306951584797e-01, 1.0171445385633245e-02},
            {9.1324172678365745e-01, 9.3149501298985338e-01, 9.7554807790077407e-03}
    };

    const double h = 2.0 / n;

    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x1 = (i + 0) * h - 1.0;
        double x2 = (i + 1) * h - 1.0;

        for (int j = 0; j < n; ++j) {
            double y1 = (j + 0) * h - 1.0;
            double y2 = (j + 1) * h - 1.0;

            for (const Node &node: weights) {
                double x = 0.5 * ((1.0 - node.i) * x1 + (1.0 + node.i) * x2);
                double y = 0.5 * ((1.0 - node.j) * y1 + (1.0 + node.j) * y2);
                res += node.w * func(map2D.get(x, y)) * map2D.Jacobian(x, y);
            }
        }
    }
    return res * sqr(h);
}

}

double Quad::volume_fraction(const std::function<bool(const Vector3d&)>& inside, int N) const {
    return integral2D::volume_fraction(*this, inside, N);
}

double Quad::integrate_low(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_low(*this, func, n);
}

double Quad::integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_mid(*this, func, n);
}

double Quad::integrate_high(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_high(*this, func, n);
}

double Quad::integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_extra(*this, func, n);
}

// ============================================================================
//                                   SQ-QUAD
// ============================================================================

SqQuad::SqQuad(
        const Vector3d &v00, const Vector3d &v01,
        const Vector3d &v10, const Vector3d &v11)
        : verts({
                        v00, 0.5 * (v00 + v01), v01,
                        0.5 * (v00 + v10), 0.25 * (v00 + v01 + v10 + v11), 0.5 * (v01 + v11),
                        v10, 0.5 * (v10 + v11), v11
                }) {}

SqQuad::SqQuad(
        const Vector3d &v00, const Vector3d &v01, const Vector3d &v02,
        const Vector3d &v10, const Vector3d &v11, const Vector3d &v12,
        const Vector3d &v20, const Vector3d &v21, const Vector3d &v22)
        : verts({
                        v00, v01, v02,
                        v10, v11, v12,
                        v20, v21, v22
                }) {}

SqQuad::SqQuad(const Quad &quad)
        : verts({
                        quad.vs<-1, -1>(), quad(0.0, -1.0), quad.vs<+1, -1>(),
                        quad(-1.0, 0.0), quad(0.0, 0.0), quad(+1.0, 0.0),
                        quad.vs<-1, +1>(), quad(0.0, +1.0), quad.vs<+1, +1>()
                }) {}

Quad SqQuad::reduce() const {
    return Quad(vs<-1, -1>(), vs<+1, -1>(),
                vs<-1, +1>(), vs<+1, +1>());
}

Vector3d SqQuad::operator()(double x, double y) const {
    return SqQuad::get(x, y);
}

Vector3d SqQuad::get(double x, double y) const {
    // Операция упрощается до одномерных в каждом направлении
    // самый алгоритмически быстрый вариант
    SqLine sqline = {
            // Сплайны по строкам матрицы
            ((SqLine *) &verts[0])->get(x),
            ((SqLine *) &verts[3])->get(x),
            ((SqLine *) &verts[6])->get(x),
    };
    return sqline.get(y);

//    // Следующие эквивалентные формулы для наглядности,
//    // также эти формулы проще дифференцировать
//    SqLine sqline = {
//            // Сплайны по строкам матрицы
//            SqLine::get(verts[0], verts[1], verts[2], x),
//            SqLine::get(verts[3], verts[4], verts[5], x),
//            SqLine::get(verts[6], verts[7], verts[8], x)
//    };
//    return sqline.get(y);
//
//    SqLine sqline = {
//            // Сплайны по столбцам матрицы
//            SqLine::get(verts[0], verts[3], verts[6], y),
//            SqLine::get(verts[1], verts[4], verts[7], y),
//            SqLine::get(verts[2], verts[5], verts[8], y)
//    };
//    return sqline.get(x);
}

Vector3d SqQuad::tangent_x(double x, double y) const {
    SqLine sqline = {
            // Сплайны по столбцам матрицы
            SqLine::get(verts[0], verts[3], verts[6], y),
            SqLine::get(verts[1], verts[4], verts[7], y),
            SqLine::get(verts[2], verts[5], verts[8], y)
    };
    return sqline.tangent(x);
}

Vector3d SqQuad::tangent_y(double x, double y) const {
    SqLine sqline = {
            // Сплайны по строкам матрицы
            ((SqLine *) &verts[0])->get(x),
            ((SqLine *) &verts[3])->get(x),
            ((SqLine *) &verts[6])->get(x),
    };
    return sqline.tangent(y);
}

Vector3d SqQuad::normal(double x, double y, const Vector3d &c) const {
    Vector3d tau_x = tangent_x(x, y);
    Vector3d tau_y = tangent_y(x, y);
    Vector3d n = tau_x.cross(tau_y).normalized();
    return (verts[iss<0, 0>()] - c).dot(n) > 0.0 ? n : -n;
}

double SqQuad::Jacobian(double x, double y) const {
    Vector3d tau_x = tangent_x(x, y);
    Vector3d tau_y = tangent_y(x, y);
    return tau_x.cross(tau_y).norm();
}

const Vector3d &SqQuad::center() const {
    return verts[iss<0, 0>()];
}

Vector3d SqQuad::normal(const Vector3d &c) const {
    // Вроде так же, как для линейного четырехугольника
    // Перепроверить, как?
    Vector3d n = (verts[3] - verts[0]).cross(verts[2] - verts[1]);
    n.normalize();
    return (verts[0] - c).dot(n) > 0.0 ? n : -n;
}

static const
std::array<std::array<int, 3>, 4> faces2D = {
        std::array<int, 3>{6, 3, 0}, // left
        std::array<int, 3>{2, 5, 8}, // right
        std::array<int, 3>{0, 1, 2}, // bottom
        std::array<int, 3>{8, 7, 6}, // top
};

double SqQuad::area() const {
    double S = 0.0;
    for (auto &face: faces2D) {
        auto &v1 = verts[face[0]];
        auto &vc = verts[face[1]];
        auto &v2 = verts[face[2]];

        S += 4.0 * (v2.y() - v1.y()) * vc.x();
        S -= 4.0 * (v2.x() - v1.x()) * vc.y();
        S += v2.y() * (3.0 * v2.x() - v1.x()) - v1.y() * (3.0 * v1.x() - v2.x());
    }
    S /= 6.0;

    return S;
}

double SqQuad::volume_as() const {
    double V = 0.0;
    for (auto &face: faces2D) {
        auto &v1 = verts[face[0]];
        auto &vc = verts[face[1]];
        auto &v2 = verts[face[2]];

        V -= (8.0 * vc.y() * vc.y() - v1.y() * v2.y()) * (v2.x() - v1.x());
        V -= (5.0 * v1.y() * v1.y() + 6.0 * v1.y() * vc.y() - v2.y() * v2.y() - 2.0 * v2.y() * vc.y()) *
             (vc.x() - v1.x());
        V -= (5.0 * v2.y() * v2.y() + 6.0 * v2.y() * vc.y() - v1.y() * v1.y() - 2.0 * v1.y() * vc.y()) *
             (v2.x() - vc.x());
    }
    V /= 30.0;

    return V;
}

Vector3d SqQuad::centroid(double area) const {
    if (area == 0.0) {
        area = SqQuad::area();
    }

    Vector3d C = {0.0, 0.0, 0.0};
    for (auto &face: faces2D) {
        auto &v1 = verts[face[0]];
        auto &vc = verts[face[1]];
        auto &v2 = verts[face[2]];

        C.x() += (8.0 * vc.x() * vc.x() - v1.x() * v2.x()) * (v2.y() - v1.y());
        C.y() -= (8.0 * vc.y() * vc.y() - v1.y() * v2.y()) * (v2.x() - v1.x());

        C.x() += (5.0 * v1.x() * v1.x() + 6.0 * v1.x() * vc.x() - v2.x() * v2.x() - 2.0 * v2.x() * vc.x()) *
                 (vc.y() - v1.y());
        C.y() -= (5.0 * v1.y() * v1.y() + 6.0 * v1.y() * vc.y() - v2.y() * v2.y() - 2.0 * v2.y() * vc.y()) *
                 (vc.x() - v1.x());

        C.x() += (5.0 * v2.x() * v2.x() + 6.0 * v2.x() * vc.x() - v1.x() * v1.x() - 2.0 * v1.x() * vc.x()) *
                 (v2.y() - vc.y());
        C.y() -= (5.0 * v2.y() * v2.y() + 6.0 * v2.y() * vc.y() - v1.y() * v1.y() - 2.0 * v1.y() * vc.y()) *
                 (v2.x() - vc.x());
    }
    C /= (30.0 * area);

    return C;
}

std::array<SqQuad, 4> SqQuad::children() const {
    // @formatter:off
    return {
            SqQuad(vs<-1, -1>(),    get(-0.5, -1.0), vs<0, -1>(),
                   get(-1.0, -0.5), get(-0.5, -0.5), get(0.0, -0.5),
                   vs<-1, 0>(),     get(-0.5, 0.0), vs<0, 0>()),

            SqQuad(vs<0, -1>(),     get(0.5, -1.0), vs<1, -1>(),
                   get(0.0, -0.5),  get(0.5, -0.5), get(1.0, -0.5),
                   vs<0, 0>(),      get(0.5, 0.0), vs<1, 0>()),

            SqQuad(vs<-1, 0>(),     get(-0.5, 0.0), vs<0, 0>(),
                   get(-1.0, 0.5),  get(-0.5, 0.5), get(0.0, 0.5),
                   vs<-1, 1>(),     get(-0.5, 1.0), vs<0, 1>()),

            SqQuad(vs<0, 0>(),      get(0.5, 0.0), vs<1, 0>(),
                   get(0.0, 0.5),   get(0.5, 0.5), get(1.0, 0.5),
                   vs<0, 1>(),      get(0.5, 1.0), vs<1, 1>())
    };
    // @formatter:on
}

double SqQuad::volume_fraction(const std::function<bool(const Vector3d&)>& inside, int N) const {
    return integral2D::volume_fraction(*this, inside, N);
}

double SqQuad::integrate_low(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_low(*this, func, n);
}

double SqQuad::integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_mid(*this, func, n);
}

double SqQuad::integrate_high(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_high(*this, func, n);
}

double SqQuad::integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral2D::integrate_extra(*this, func, n);
}

} // namespace zephyr::geom