#include <zephyr/geom/maps.h>

namespace zephyr::geom {

// Перпендикуляр единичной длины 'n' к вектору 'tau', лежащий в плоскости
// векторов 'a' и 'tau', сонаправленный с вектором 'a', то есть
// (n, n) = 1       // единичный
// (tau, n) = 0     // перпендикуляр
// (tau, n, a) = 0  // в одной плоскости
// (tau, a) > 0     // сонаправлены
inline Vector3d perpendicular(const Vector3d &tau, const Vector3d &a) {
    return tau.cross(a).cross(tau).normalized();
}

// ============================================================================
//                                    LINE
// ============================================================================

Line::Line(const Vector3d &v1, const Vector3d &v2)
        : verts{v1, v2} {}

Vector3d Line::operator()(double x) const {
    return get(verts[0], verts[1], x);
}

Vector3d Line::get(const Vector3d &v1, const Vector3d &v2, double x) {
    return 0.5 * (1.0 - x) * v1 + 0.5 * (1.0 + x) * v2;
}

Vector3d Line::normal(const Vector3d &v1, const Vector3d &v2, const Vector3d &c) {
    return perpendicular(v2 - v1, v1 - c);
}

double Line::Jacobian() const {
    return 0.5 * (verts[1] - verts[0]).norm();
}

Vector3d Line::center() const {
    return 0.5 * (verts[0] + verts[1]);
}

Vector3d Line::normal(const Vector3d &c) const {
    return perpendicular(verts[1] - verts[0], verts[0] - c);
}

double Line::length() const {
    return (verts[1] - verts[0]).norm();
}

// ============================================================================
//                                  SQ-LINE
// ============================================================================

SqLine::SqLine(const Vector3d &v1, const Vector3d &v2)
        : verts({v1, 0.5 * (v1 + v2), v2}) {}

SqLine::SqLine(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3)
        : verts({v1, v2, v3}) {}

SqLine::SqLine(const Line &vs)
        : verts({vs[0], vs.center(), vs[1]}) {}

Vector3d SqLine::get(const Vector3d &v1, const Vector3d &vc, const Vector3d &v2, double x) {
    return vc + 0.5 * (v2 - v1) * x + 0.5 * (v2 - 2.0 * vc + v1) * x * x;
}

Vector3d SqLine::tangent(const Vector3d &v1, const Vector3d &vc, const Vector3d &v2, double x) {
    return 0.5 * (v2 - v1) + (v2 - 2.0 * vc + v1) * x;
}

Vector3d SqLine::projection(const Vector3d &a) const {
    const double eps = 1.0e-14;
    double L2 = (verts[2] - verts[0]).squaredNorm();

    // Перпендикуляр к плоскости, в которой лежит кривая
    // (ноль, если точки лежат на одной линии)
    Vector3d np = (verts[2] - verts[1]).cross(verts[0] - verts[1]);
    double np_norm = np.norm();

    if (np_norm < eps * L2) {
        // Считаем, что точки лежат на одной линии
        return a;
    } else {
        // Проекция точки 'a' на плоскость кривой
        np /= np_norm;
        return np.cross(a).cross(np);
    }
}

Vector3d SqLine::operator()(double x) const {
    return SqLine::get(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::get(double x) const {
    return SqLine::get(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::tangent(double x) const {
    return SqLine::tangent(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::normal(double x, const Vector3d &c) const {
    Vector3d a = projection(verts[1] - c);
    return perpendicular(tangent(x), a);
}

double SqLine::Jacobian(double x) const {
    return SqLine::tangent(verts[0], verts[1], verts[2], x).norm();
}

const Vector3d &SqLine::center() const {
    return verts[1];
}

Vector3d SqLine::normal(const Vector3d &c) const {
    // Это не заглушка, было установлено, что в расчете потока через
    // криволинейную грань следует использовать такую же нормаль,
    // как при расчете через обычную грань.
    // Единственное, здесь добавляется проекция на плоскость кривой.

    Vector3d a = projection(verts[1] - c);
    return perpendicular(verts[2] - verts[0], a);
}

double SqLine::length() const {
    // Это не заглушка, было установлено, что в расчете потока
    // по криволинейной грани следует использовать такую длину
    return (verts[2] - verts[0]).norm();
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

double Quad::volume_as() const {
    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    double V = 0.0;
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = verts[ord[i]];
        auto &v2 = verts[ord[(i + 1) % 4]];

        V -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    V /= 6.0;

    return V;
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

// ============================================================================
//                                  SQ-QUAD
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

// ============================================================================
//                                    CUBE
// ============================================================================

Cube::Cube(
        const Vector3d &v000, const Vector3d &v001,
        const Vector3d &v010, const Vector3d &v011,
        const Vector3d &v100, const Vector3d &v101,
        const Vector3d &v110, const Vector3d &v111)
        : verts({v000, v001, v010, v011, v100, v101, v110, v111}) {}

Vector3d Cube::operator()(double x, double y, double z) const {
    return Cube::get(x, y, z);
}

Vector3d Cube::get(double x, double y, double z) const {
    Vector3d res = {0.0, 0.0, 0.0};
    res += ((1 - x) * (1 - y) * (1 - z)) * vs<-1, -1, -1>();
    res += ((1 + x) * (1 - y) * (1 - z)) * vs<+1, -1, -1>();
    res += ((1 - x) * (1 + y) * (1 - z)) * vs<-1, +1, -1>();
    res += ((1 + x) * (1 + y) * (1 - z)) * vs<+1, +1, -1>();
    res += ((1 - x) * (1 - y) * (1 + z)) * vs<-1, -1, +1>();
    res += ((1 + x) * (1 - y) * (1 + z)) * vs<+1, -1, +1>();
    res += ((1 - x) * (1 + y) * (1 + z)) * vs<-1, +1, +1>();
    res += ((1 + x) * (1 + y) * (1 + z)) * vs<+1, +1, +1>();
    res /= 8.0;
    return res;
}

Vector3d Cube::tangent_x(double x, double y, double z) const {
    Vector3d res = {0.0, 0.0, 0.0};
    res += (1 - y) * (1 - z) * (vs<+1, -1, -1>() - vs<-1, -1, -1>());
    res += (1 + y) * (1 - z) * (vs<+1, +1, -1>() - vs<-1, +1, -1>());
    res += (1 - y) * (1 + z) * (vs<+1, -1, +1>() - vs<-1, -1, +1>());
    res += (1 + y) * (1 + z) * (vs<+1, +1, +1>() - vs<-1, +1, +1>());
    res /= 8.0;
    return res;
}

Vector3d Cube::tangent_y(double x, double y, double z) const {
    Vector3d res = {0.0, 0.0, 0.0};
    res += (1 - x) * (1 - z) * (vs<-1, +1, -1>() - vs<-1, -1, -1>());
    res += (1 + x) * (1 - z) * (vs<+1, +1, -1>() - vs<+1, -1, -1>());
    res += (1 - x) * (1 + z) * (vs<-1, +1, +1>() - vs<-1, -1, +1>());
    res += (1 + x) * (1 + z) * (vs<+1, +1, +1>() - vs<+1, -1, +1>());
    res /= 8.0;
    return res;
}

Vector3d Cube::tangent_z(double x, double y, double z) const {
    Vector3d res = {0.0, 0.0, 0.0};
    res += (1 - x) * (1 - y) * (vs<-1, -1, +1>() - vs<-1, -1, -1>());
    res += (1 + x) * (1 - y) * (vs<+1, -1, +1>() - vs<+1, -1, -1>());
    res += (1 - x) * (1 + y) * (vs<-1, +1, +1>() - vs<-1, +1, -1>());
    res += (1 + x) * (1 + y) * (vs<+1, +1, +1>() - vs<+1, +1, -1>());
    res /= 8.0;
    return res;
}

double Cube::Jacobian(double x, double y, double z) const {
    Vector3d tau_x = tangent_x(x, y, z);
    Vector3d tau_y = tangent_y(x, y, z);
    Vector3d tau_z = tangent_z(x, y, z);
    return tau_x.cross(tau_y).dot(tau_z);
}

Vector3d Cube::center() const {
    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < 8; ++i) {
        C += verts[i];
    }
    return C / 8.0;
}

double Cube::volume() const {
    // Действительно так сложно?
    Vector3d cell_c = Cube::center();

    std::array<Quad, 6> faces = {
            Quad({verts[4], verts[0], verts[6], verts[2]}), // left
            Quad({verts[5], verts[1], verts[7], verts[3]}), // right
            Quad({verts[4], verts[5], verts[0], verts[1]}), // bottom
            Quad({verts[6], verts[7], verts[2], verts[3]}), // top
            Quad({verts[0], verts[1], verts[2], verts[3]}), // back
            Quad({verts[4], verts[5], verts[6], verts[7]})  // front
    };

    double V = 0.0;
    for (auto &face: faces) {
        Vector3d face_c = face.center();
        Vector3d face_s = face.area() * face.normal(cell_c);
        V += (face_c - cell_c).dot(face_s);
    }
    V /= 3.0;

    return V;
}

Vector3d Cube::centroid(double volume) const {
    if (volume == 0.0) {
        volume = Cube::volume();
    }

    // Аппроксимация объемного интеграла методом Гаусса
    // демонстрирует огромную точность
    double a = 1.0 / std::sqrt(3.0);

    Vector3d C = {0.0, 0.0, 0.0};
    C += get(-a, -a, -a) * Jacobian(-a, -a, -a);
    C += get(+a, -a, -a) * Jacobian(+a, -a, -a);
    C += get(-a, +a, -a) * Jacobian(-a, +a, -a);
    C += get(+a, +a, -a) * Jacobian(+a, +a, -a);
    C += get(-a, -a, +a) * Jacobian(-a, -a, +a);
    C += get(+a, -a, +a) * Jacobian(+a, -a, +a);
    C += get(-a, +a, +a) * Jacobian(-a, +a, +a);
    C += get(+a, +a, +a) * Jacobian(+a, +a, +a);
    C /= volume;

    return C;
}

// ============================================================================
//                                  SQ-CUBE
// ============================================================================

SqCube::SqCube(const Vector3d &v000,
               const Vector3d &v002,
               const Vector3d &v020,
               const Vector3d &v022,
               const Vector3d &v200,
               const Vector3d &v202,
               const Vector3d &v220,
               const Vector3d &v222)
        : verts({
            // О ма гадабл, ты что крейзи?
                        v000, 0.5 * (v000 + v002), v002,
                        0.5 * (v000 + v020), 0.25 * (v000 + v002 + v020 + v022), 0.5 * (v002 + v022),
                        v020, 0.5 * (v020 + v022), v022,
                        0.5 * (v000 + v200), 0.25 * (v000 + v200 + v002 + v202), 0.5 * (v002 + v202),
                        0.25 * (v000 + v200 + v020 + v220),
                        0.125 * (v000 + v002 + v020 + v022 + v200 + v202 + v220 + v222),
                        0.25 * (v002 + v202 + v022 + v222),
                        0.5 * (v020 + v220), 0.25 * (v020 + v220 + v022 + v222), 0.5 * (v022 + v222),
                        v200, 0.5 * (v200 + v202), v202,
                        0.5 * (v200 + v220), 0.25 * (v200 + v202 + v220 + v222), 0.5 * (v202 + v222),
                        v220, 0.5 * (v220 + v222), v222
                }) {}

SqCube::SqCube(
        const Vector3d &v000, const Vector3d &v001, const Vector3d &v002,
        const Vector3d &v010, const Vector3d &v011, const Vector3d &v012,
        const Vector3d &v020, const Vector3d &v021, const Vector3d &v022,
        const Vector3d &v100, const Vector3d &v101, const Vector3d &v102,
        const Vector3d &v110, const Vector3d &v111, const Vector3d &v112,
        const Vector3d &v120, const Vector3d &v121, const Vector3d &v122,
        const Vector3d &v200, const Vector3d &v201, const Vector3d &v202,
        const Vector3d &v210, const Vector3d &v211, const Vector3d &v212,
        const Vector3d &v220, const Vector3d &v221, const Vector3d &v222)
        : verts({
                        v000, v001, v002, v010, v011, v012, v020, v021, v022,
                        v100, v101, v102, v110, v111, v112, v120, v121, v122,
                        v200, v201, v202, v210, v211, v212, v220, v221, v222
                }) {}

SqCube::SqCube(const Cube &cube)
        : SqCube(cube[0], cube[1], cube[2], cube[3],
                 cube[4], cube[5], cube[6], cube[7]) {}

SqCube::SqCube(const Quad& quad)
    : SqCube(SqQuad(quad)) { }

SqCube::SqCube(const SqQuad& quad) {
    std::memcpy((void *) verts.data(), (void *) &quad, 9 * sizeof(Vector3d));

    const Vector3d nanvec = {NAN, NAN, NAN};
    std::fill(verts.begin() + 9, verts.end(), nanvec);
}

Cube SqCube::reduce() const {
    return Cube(vs<-1, -1, -1>(), vs<+1, -1, -1>(),
                vs<-1, +1, -1>(), vs<+1, +1, -1>(),
                vs<-1, -1, +1>(), vs<+1, -1, +1>(),
                vs<-1, +1, +1>(), vs<+1, +1, +1>());
}

using LargeGrid3D = std::array<std::array<std::array<Vector3d, 5>, 5>, 5>;

template<int i, int j, int k>
inline SqCube sq_cube_from_table(const LargeGrid3D& grid) {
    static_assert(i == 0 || i == 2);
    static_assert(j == 0 || j == 2);
    static_assert(k == 0 || k == 2);

    return {
            grid[i + 0][j + 0][k + 0], grid[i + 1][j + 0][k + 0], grid[i + 2][j + 0][k + 0],
            grid[i + 0][j + 1][k + 0], grid[i + 1][j + 1][k + 0], grid[i + 2][j + 1][k + 0],
            grid[i + 0][j + 2][k + 0], grid[i + 1][j + 2][k + 0], grid[i + 2][j + 2][k + 0],

            grid[i + 0][j + 0][k + 1], grid[i + 1][j + 0][k + 1], grid[i + 2][j + 0][k + 1],
            grid[i + 0][j + 1][k + 1], grid[i + 1][j + 1][k + 1], grid[i + 2][j + 1][k + 1],
            grid[i + 0][j + 2][k + 1], grid[i + 1][j + 2][k + 1], grid[i + 2][j + 2][k + 1],

            grid[i + 0][j + 0][k + 2], grid[i + 1][j + 0][k + 2], grid[i + 2][j + 0][k + 2],
            grid[i + 0][j + 1][k + 2], grid[i + 1][j + 1][k + 2], grid[i + 2][j + 1][k + 2],
            grid[i + 0][j + 2][k + 2], grid[i + 1][j + 2][k + 2], grid[i + 2][j + 2][k + 2]
    };
}

std::array<SqCube, 8> SqCube::children() const {
    LargeGrid3D grid;
    for (int i = 0; i < 5; ++i) {
        double x = 0.5 * (i - 2);
        for (int j = 0; j < 5; ++j) {
            double y = 0.5 * (j - 2);
            for (int k = 0; k < 5; ++k) {
                double z = 0.5 * (k - 2);
                grid[i][j][k] = get(x, y, z);
            }
        }
    }

    return {
            sq_cube_from_table<0, 0, 0>(grid),
            sq_cube_from_table<2, 0, 0>(grid),
            sq_cube_from_table<0, 2, 0>(grid),
            sq_cube_from_table<2, 2, 0>(grid),
            sq_cube_from_table<0, 0, 2>(grid),
            sq_cube_from_table<2, 0, 2>(grid),
            sq_cube_from_table<0, 2, 2>(grid),
            sq_cube_from_table<2, 2, 2>(grid),
    };
}

Vector3d SqCube::operator()(double x, double y, double z) const {
    return SqCube::get(x, y, z);
}

Vector3d SqCube::operator()(double x, double y) const {
    return as2D().get(x, y);
}

Vector3d SqCube::get(double x, double y, double z) const {
    // Операция упрощается до одномерных в каждом направлении
    // самый алгоритмически быстрый вариант
    SqLine sqline = {
            // Сплайн по двумерным срезам
            ((SqQuad *) &verts[ 0])->get(x, y),
            ((SqQuad *) &verts[ 9])->get(x, y),
            ((SqQuad *) &verts[18])->get(x, y),
    };
    return sqline.get(z);
}

Vector3d SqCube::tangent_x(double x, double y, double z) const {
    // Для квадратичных зависимостей центральная разность
    // дает точный результат
    return get(x + 0.5, y, z) - get(x - 0.5, y, z);
}

Vector3d SqCube::tangent_y(double x, double y, double z) const {
    // Для квадратичных зависимостей центральная разность
    // дает точный результат
    return get(x, y + 0.5, z) - get(x, y - 0.5, z);
}

Vector3d SqCube::tangent_z(double x, double y, double z) const {
    // Для квадратичных зависимостей центральная разность
    // дает точный результат
    return get(x, y, z + 0.5) - get(x, y, z - 0.5);
}

double SqCube::Jacobian(double x, double y, double z) const {
    Vector3d tau_x = tangent_x(x, y, z);
    Vector3d tau_y = tangent_y(x, y, z);
    Vector3d tau_z = tangent_z(x, y, z);
    return tau_x.cross(tau_y).dot(tau_z);
}

Vector3d SqCube::center() const {
    return vs<0, 0, 0>();
}

double SqCube::volume() const {
    // Заплатка для линейных кубов
    return reduce().volume();
}

Vector3d SqCube::centroid(double volume) const {
    // Заплатка для линейных кубов
    return reduce().centroid(volume);
}

} // namespace zephyr::geom