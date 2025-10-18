#include <cstring>
#include <iostream>

#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>

namespace zephyr::geom {

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

// Реализация интегралов по Quad и SqQuad
namespace integral3D {

// type Map3D is Cube or SqCube
// N -- число ячеек, точность определения объемной доли ~ 1/N
template <typename Map3D>
double volume_fraction(const Map3D& cube, const std::function<bool(const Vector3d&)>& inside, int N) {
    int n = std::max(int(std::round(std::cbrt(N))), 1);
    double h = 1.0 / n;
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = (2 * i + 1) * h - 1.0;
        for (int j = 0; j < n; ++j) {
            double y = (2 * j + 1) * h - 1.0;
            for (int k = 0; k < n; ++k) {
                double z = (2 * k + 1) * h - 1.0;
                if (inside(cube.get(x, y, z))) {
                    res += cube.Jacobian(x, y, z);
                }
            }
        }
    }
    return res * std::pow(2.0 * h, 3) / cube.volume();
}

// type Map3D is Cube or SqCube
template <typename Map3D>
double integrate_low(const Map3D& cube, const std::function<double(const Vector3d &)> &func, int n) {
    double h = 1.0 / n;
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = (2 * i + 1) * h - 1.0;
        for (int j = 0; j < n; ++j) {
            double y = (2 * j + 1) * h - 1.0;
            for (int k = 0; k < n; ++k) {
                double z = (2 * k + 1) * h - 1.0;
                res += func(cube.get(x, y, z)) * cube.Jacobian(x, y, z);
            }
        }
    }
    return res * std::pow(2.0 * h, 3);
}

// type Map3D is Cube or SqCube
template <typename Map3D>
double integrate_mid(const Map3D& cube, const std::function<double(const Vector3d &)> &func, int n) {
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
            
            for (int k = 0; k < n; ++k) {
                double z1 = (2 * k + cm) * h - 1.0;
                double z2 = (2 * k + cp) * h - 1.0;

                res += func(cube.get(x1, y1, z1)) * cube.Jacobian(x1, y1, z1) +
                       func(cube.get(x1, y2, z1)) * cube.Jacobian(x1, y2, z1) +
                       func(cube.get(x2, y1, z1)) * cube.Jacobian(x2, y1, z1) +
                       func(cube.get(x2, y2, z1)) * cube.Jacobian(x2, y2, z1) +
                       func(cube.get(x1, y1, z2)) * cube.Jacobian(x1, y1, z2) +
                       func(cube.get(x1, y2, z2)) * cube.Jacobian(x1, y2, z2) +
                       func(cube.get(x2, y1, z2)) * cube.Jacobian(x2, y1, z2) +
                       func(cube.get(x2, y2, z2)) * cube.Jacobian(x2, y2, z2);
            }
        }
    }
    return res * std::pow(h, 3);
}

} // integral3D

double Cube::volume_fraction(const std::function<bool(const Vector3d&)>& inside, int n_points) const {
    return integral3D::volume_fraction(*this, inside, n_points);
}

double Cube::integrate_low(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral3D::integrate_low(*this, func, n);
}

double Cube::integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const {
    return integral3D::integrate_mid(*this, func, n);
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