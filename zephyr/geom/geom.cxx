#include <zephyr/geom/geom.h>
#include <zephyr/geom/maps.h>

namespace zephyr { namespace geom {

double length(const ShortList1D& vs) {
    return (vs[0] - vs[1]).norm();
}

double length(const LargeList1D& vs) {
    // Это не заглушка, было установлено, что в расчете потока
    // по криволинейной грани следует использовать такую длину
    return (vs[0] - vs[2]).norm();
}

Vector3d normal(const ShortList1D& vs, const Vector3d& center) {
    Vector3d n = {vs[1].y() - vs[0].y(), vs[0].x() - vs[1].x(), 0.0};
    n.normalize();
    return (center - vs[0]).dot(n) < 0.0 ? n : -n;
}

Vector3d normal(const LargeList1D& vs, const Vector3d& center) {
    // Это не заглушка, было установлено, что в расчете потока
    // по криволинейной грани следует использовать такую нормаль
    Vector3d n = {vs[2].y() - vs[0].y(), vs[0].x() - vs[2].x(), 0.0};
    n.normalize();
    return (center - vs[0]).dot(n) < 0.0 ? n : -n;
}

Vector3d normal(const ShortList2D& vs, const Vector3d& center) {
    Vector3d n = (vs[3] - vs[0]).cross(vs[2] - vs[1]);
    n.normalize();
    return (center - vs[0]).dot(n) < 0.0 ? n : -n;
}

double area(const ShortList2D &vs) {
    return 0.5 * (vs[3] - vs[0]).cross(vs[2] - vs[1]).norm();
}

double area(const LargeList2D &vs) {
    std::array<std::array<Vector3d, 3>, 4> faces = {
            std::array<Vector3d, 3>({vs[6], vs[3], vs[0]}), // left
            std::array<Vector3d, 3>({vs[2], vs[5], vs[8]}), // right
            std::array<Vector3d, 3>({vs[0], vs[1], vs[2]}), // bottom
            std::array<Vector3d, 3>({vs[8], vs[7], vs[6]})  // top
    };

    double S = 0.0;
    for (auto& face: faces) {
        auto &v1 = face[0];
        auto &vc = face[1];
        auto &v2 = face[2];

        S += 4.0 * (v2.y() - v1.y()) * vc.x();
        S -= 4.0 * (v2.x() - v1.x()) * vc.y();
        S += v2.y() * (3.0 * v2.x() - v1.x()) - v1.y() * (3.0 * v1.x() - v2.x());
    }
    S /= 6.0;

    return S;
}

double area(const VerticesList& vs) {
    int n_points = vs.size();
    Vector3d cell_c = center(vs);

    double S = 0.0;
    for (int i = 0; i < n_points; ++i) {
        const Vector3d &v1 = vs[i];
        const Vector3d &v2 = vs[(i + 1) % n_points];

        Vector3d face_c = 0.5 * (v1 + v2);
        Vector3d normal = {v2.y() - v1.y(), v1.x() - v2.x(), 0.0};

        S += std::abs((face_c - cell_c).dot(normal));
    }
    S /= 2.0;

    return S;
}

double volume_as(const ShortList2D &vs) {
    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    double V = 0.0;
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = vs[ord[i]];
        auto &v2 = vs[ord[(i + 1) % 4]];

        V -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    V /= 6.0;

    return V;
}

double volume_as(const LargeList2D &vs) {
    std::array<std::array<Vector3d, 3>, 4> faces = {
            std::array<Vector3d, 3>({vs[6], vs[3], vs[0]}), // left
            std::array<Vector3d, 3>({vs[2], vs[5], vs[8]}), // right
            std::array<Vector3d, 3>({vs[0], vs[1], vs[2]}), // bottom
            std::array<Vector3d, 3>({vs[8], vs[7], vs[6]})  // top
    };

    double V = 0.0;
    for (auto& face: faces) {
        auto &v1 = face[0];
        auto &vc = face[1];
        auto &v2 = face[2];

        V -= (8.0 * vc.y() * vc.y() - v1.y() * v2.y()) * (v2.x() - v1.x());
        V -= (5.0 * v1.y() * v1.y() + 6.0 * v1.y() * vc.y() - v2.y() * v2.y() - 2.0 * v2.y() * vc.y()) * (vc.x() - v1.x());
        V -= (5.0 * v2.y() * v2.y() + 6.0 * v2.y() * vc.y() - v1.y() * v1.y() - 2.0 * v1.y() * vc.y()) * (v2.x() - vc.x());
    }
    V /= 30.0;

    return V;
}

double volume(const ShortList3D &vs) {
    Vector3d cell_c = center(vs);

    std::array<ShortList2D, 6> faces = {
            ShortList2D({vs[4], vs[0], vs[6], vs[2]}), // left
            ShortList2D({vs[5], vs[1], vs[7], vs[3]}), // right
            ShortList2D({vs[4], vs[5], vs[0], vs[1]}), // bottom
            ShortList2D({vs[6], vs[7], vs[2], vs[3]}), // top
            ShortList2D({vs[0], vs[1], vs[2], vs[3]}), // back
            ShortList2D({vs[4], vs[5], vs[6], vs[7]})  // front
    };

    double V = 0.0;
    for (auto &face: faces) {
        Vector3d face_c = center(face);
        Vector3d face_s = area(face) * normal(face, cell_c);
        V += (face_c - cell_c).dot(face_s);
    }
    V /= 3.0;

    return V;
}

Vector3d center(const ShortList1D &vs) {
    return (vs[0] + vs[1]) / 2.0;
}

Vector3d center(const LargeList1D &vs) {
    return vs[1];
}
Vector3d center(const ShortList2D &vs) {
    return (vs[0] + vs[1] + vs[2] + vs[3]) / 4.0;
}

Vector3d center(const LargeList2D &vs) {
    return (vs[1] + vs[3] + vs[5] + vs[7]) / 4.0;
}

Vector3d center(const VerticesList& vs) {
    Vector3d C = {0.0, 0.0, 0.0};
    for (auto& v: vs) {
        C += v;
    }
    return C / vs.size();
}

Vector3d center(const ShortList3D &vs) {
    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < 8; ++i) {
        C += vs[i];
    }
    return C / 8.0;
}

Vector3d centroid(const ShortList2D &vs, double area) {
    if (area == 0.0) {
        area = geom::area(vs);
    }

    // Обход вершин против часовой стрелки
    int ord[4] = {0, 1, 3, 2};

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i: {0, 1, 2, 3}) {
        auto &v1 = vs[ord[i]];
        auto &v2 = vs[ord[(i + 1) % 4]];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * area);

    return C;
}

Vector3d centroid(const VerticesList& vs, double area) {
    if (area == 0.0) {
        area = geom::area(vs);
    }
    int n_points = vs.size();

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < n_points; ++i) {
        auto& v1 = vs[i];
        auto& v2 = vs[(i + 1) % n_points];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * area);

    return C;
}

Vector3d centroid(const LargeList2D &vs, double area) {
    if (area == 0.0) {
        area = geom::area(vs);
    }

    std::array<std::array<Vector3d, 3>, 4> faces = {
            std::array<Vector3d, 3>({vs[6], vs[3], vs[0]}), // left
            std::array<Vector3d, 3>({vs[2], vs[5], vs[8]}), // right
            std::array<Vector3d, 3>({vs[0], vs[1], vs[2]}), // bottom
            std::array<Vector3d, 3>({vs[8], vs[7], vs[6]})  // top
    };

    Vector3d C = {0.0, 0.0, 0.0};
    for (auto& face: faces) {
        auto &v1 = face[0];
        auto &vc = face[1];
        auto &v2 = face[2];

        C.x() += (8.0 * vc.x() * vc.x() - v1.x() * v2.x()) * (v2.y() - v1.y());
        C.y() -= (8.0 * vc.y() * vc.y() - v1.y() * v2.y()) * (v2.x() - v1.x());

        C.x() += (5.0 * v1.x() * v1.x() + 6.0 * v1.x() * vc.x() - v2.x() * v2.x() - 2.0 * v2.x() * vc.x()) * (vc.y() - v1.y());
        C.y() -= (5.0 * v1.y() * v1.y() + 6.0 * v1.y() * vc.y() - v2.y() * v2.y() - 2.0 * v2.y() * vc.y()) * (vc.x() - v1.x());

        C.x() += (5.0 * v2.x() * v2.x() + 6.0 * v2.x() * vc.x() - v1.x() * v1.x() - 2.0 * v1.x() * vc.x()) * (v2.y() - vc.y());
        C.y() -= (5.0 * v2.y() * v2.y() + 6.0 * v2.y() * vc.y() - v1.y() * v1.y() - 2.0 * v1.y() * vc.y()) * (v2.x() - vc.x());
    }
    C /= (30.0 * area);

    return C;
}

Vector3d centroid(const ShortList3D &vs, double volume) {
    if (volume == 0.0) {
        volume = geom::volume(vs);
    }

    // Аппроксимация объемного интеграла методом Гаусса
    // демонстрирует огромную точность
    Mapping3D map3D(vs);
    double a = 1.0 / std::sqrt(3.0);

    Vector3d C = {0.0, 0.0, 0.0};
    C += map3D(-a, -a, -a) * map3D.Jacobian(-a, -a, -a);
    C += map3D(+a, -a, -a) * map3D.Jacobian(+a, -a, -a);
    C += map3D(-a, +a, -a) * map3D.Jacobian(-a, +a, -a);
    C += map3D(+a, +a, -a) * map3D.Jacobian(+a, +a, -a);
    C += map3D(-a, -a, +a) * map3D.Jacobian(-a, -a, +a);
    C += map3D(+a, -a, +a) * map3D.Jacobian(+a, -a, +a);
    C += map3D(-a, +a, +a) * map3D.Jacobian(-a, +a, +a);
    C += map3D(+a, +a, +a) * map3D.Jacobian(+a, +a, +a);
    C /= volume;

    return C;
}

} // geom
} // zephyr