#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/surface/solid_body_3d.h>

namespace zephyr::geom {
// двумерная матрица поворота на угол phi
inline Matrix3d rotation_matrix(const Vector3d& axis, double phi) {
    Vector3d n = axis.normalized();
    double x = n.x();
    double y = n.y();
    double z = n.z();

    double cos = std::cos(phi);
    double sin = std::sin(phi);
    Matrix3d R;
    R << cos + (1 - cos)*x*x, (1 - cos)*x*y - sin*z, (1-cos)*x*z + sin*y,
         (1 - cos)*y*x + sin*z, cos + (1 - cos)*y*y, (1-cos)*y*z - sin*x,
         (1 - cos)*z*x - sin*y, (1 - cos)*z*y + sin*x, cos + (1 - cos)*z*z;
    return R;
}

SolidBody3D::SolidBody3D() = default;

SolidBody3D::SolidBody3D(const Vector3d& center)
    : SolidBody(center) { }

void SolidBody3D::rotation_local(const Vector3d& axis, double angle) {
    SolidBody::rotation_local(rotation_matrix(axis, angle));
}

void SolidBody3D::rotation_global(const Vector3d& axis, double angle) {
    SolidBody::rotation_global(rotation_matrix(axis, angle));
}

void SolidBody3D::rotation_relative(const Vector3d& axis, double angle, const Vector3d& c) {
    SolidBody::rotation_relative(rotation_matrix(axis, angle), c);
}

std::vector<std::array<Vector3d, 3>> SolidBody3D::triangulation(int n_elements) const {
    auto tris = local_triangulation(n_elements);
    for (auto& tri: tris) {
        for (auto& v: tri) {
            v.applyOnTheLeft(m_rotation);
            v += m_center;
        }
    }
    return tris;
}

// ============================================================================
//                                 BALL
// ============================================================================

void BodyBall::init_inside() {
    m_inside = [this](const Vector3d &v) -> bool {
        return (v - m_center).norm() < m_radius;
    };
}

BodyBall::BodyBall(double radius) : m_radius(radius) {
    init_inside();
}

BodyBall::BodyBall(double radius, const Vector3d& center)
    : SolidBody3D(center), m_radius(radius) {
    init_inside();
}

inline std::vector<std::array<Vector3d, 3>> icosahedron() {
    const double phi = 0.5 * (1.0 + std::sqrt(5.0));

    Vector3d vs[] = {
        {-1,  phi, 0}, // 0
        { 1,  phi, 0}, // 1
        {-1, -phi, 0}, // 2
        { 1, -phi, 0}, // 3
        {0, -1,  phi}, // 4
        {0,  1,  phi}, // 5
        {0, -1, -phi}, // 6
        {0,  1, -phi}, // 7
        { phi, 0, -1}, // 8
        { phi, 0,  1}, // 9
        {-phi, 0, -1}, //10
        {-phi, 0,  1}  //11
    };

    for (auto& vv: vs) {
        vv.normalize();
    }

    return {
        {vs[0],  vs[11], vs[5]},
        {vs[0],  vs[5],  vs[1]},
        {vs[0],  vs[1],  vs[7]},
        {vs[0],  vs[7],  vs[10]},
        {vs[0],  vs[10], vs[11]},
        {vs[1],  vs[5],  vs[9]},
        {vs[5],  vs[11], vs[4]},
        {vs[11], vs[10], vs[2]},
        {vs[10], vs[7],  vs[6]},
        {vs[7],  vs[1],  vs[8]},
        {vs[3],  vs[9],  vs[4]},
        {vs[3],  vs[4],  vs[2]},
        {vs[3],  vs[2],  vs[6]},
        {vs[3],  vs[6],  vs[8]},
        {vs[3],  vs[8],  vs[9]},
        {vs[4],  vs[9],  vs[5]},
        {vs[2],  vs[4],  vs[11]},
        {vs[6],  vs[2],  vs[10]},
        {vs[8],  vs[6],  vs[7]},
        {vs[9],  vs[8],  vs[1]}
    };
}

void subdivide_triangle(const std::array<Vector3d, 3>& tri, int n,
    std::vector<std::array<Vector3d, 3>>& triangles) {

    triangles.reserve(triangles.size() + n * n);

    // Вершины исходного треугольника
    const Vector3d& A = tri[0];
    const Vector3d& B = tri[1];
    const Vector3d& C = tri[2];

    // Генерация сетки точек на треугольнике
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n - i; ++j) {
            Vector3d p1 = A + ((i + 0.0) / n) * (B - A)
                            + ((j + 0.0) / n) * (C - A);

            Vector3d p2 = A + ((i + 1.0) / n) * (B - A)
                            + ((j + 0.0) / n) * (C - A);

            Vector3d p3 = A + ((i + 0.0) / n) * (B - A)
                            + ((j + 1.0) / n) * (C - A);

            triangles.push_back({p1, p2, p3});
            if (i + j < n - 1) {
                Vector3d p4 = A + ((i + 1.0) / n) * (B - A)
                                + ((j + 1.0) / n) * (C - A);
                triangles.push_back({p2, p4, p3});
            }
        }
    }
}

std::vector<std::array<Vector3d, 3>> BodyBall::local_triangulation(int n_elements) const {
    const auto ico = icosahedron();

    std::vector<std::array<Vector3d, 3>> triangles;

    int n = std::floor(std::sqrt(n_elements / 20.0));
    if (n <= 1) {
        triangles = ico;
    }
    else {
        triangles.reserve(20 * n * n);
        for (const auto& tri: ico) {
            subdivide_triangle(tri, n, triangles);
        }
    }

    for (auto& tri: triangles) {
        for (auto& v: tri) {
            v.normalize();
            v *= m_radius;
        }
    }
    return triangles;
}

// ============================================================================
//                                  CUBE
// ============================================================================

void BodyCube::init_inside() {
    m_inside = [this](const Vector3d &v) -> bool {
        Vector3d p = this->in_local(v);
        return p.cwiseAbs().maxCoeff() < 0.5 * m_length;
    };
}

BodyCube::BodyCube(double length) : m_length(length) {
    init_inside();
}

BodyCube::BodyCube(double length, const Vector3d& center)
    : SolidBody3D(center), m_length(length) {
    init_inside();
}

std::vector<std::array<Vector3d, 3>> BodyCube::local_triangulation(int n_elements) const {
    double a = 0.5 * m_length;
    const Vector3d v[8] = {
        Vector3d{-a, -a, -a},
        Vector3d{ a, -a, -a},
        Vector3d{-a,  a, -a},
        Vector3d{ a,  a, -a},
        Vector3d{-a, -a,  a}, 
        Vector3d{ a, -a,  a},
        Vector3d{-a,  a,  a},
        Vector3d{ a,  a,  a},
    };

    constexpr int indices[12][3] = {
        {4, 5, 7}, {4, 7, 6},
        {0, 2, 3}, {0, 3, 1},
        {1, 3, 7}, {1, 7, 5},
        {0, 4, 6}, {0, 6, 2},
        {2, 6, 7}, {2, 7, 3},
        {0, 1, 5}, {0, 5, 4}
    };

    std::vector<std::array<Vector3d, 3>> triangles;
    triangles.reserve(12);
    
    for (const auto& idx: indices) {
        triangles.push_back({v[idx[0]], v[idx[1]], v[idx[2]]});
    }
    return triangles;
}

} // namespace zephyr::geom