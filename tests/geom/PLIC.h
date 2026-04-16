// Header для тестов PLIC_2D и PLIC_3D
#pragma once

#include <iomanip>
#include <atomic>

#include <zephyr/math/funcs.h>

#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/io/pvd_file.h>

using namespace zephyr;
using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Cuboid;
using generator::Rectangle;

// Радиус и угол в декартовы координаты
inline Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

// Характеристическая функция (функция-индикатор)
using InFunction = std::function<bool(const Vector3d &)>;

// Пространственная функция
using SpFunction = std::function<double(const Vector3d &)>;

using triplet = std::array<double, 3>;
using table_3x3 = std::array<triplet, 3>;
using table_3x3x3 = std::array<table_3x3, 3>;

// Конечно-разностные производные на сетке из трёх узлов
inline triplet derivatives(const triplet& f, double h) {
    return {
        (f[1] - f[0]) / h,
        (f[2] - f[0]) / (2 * h),
        (f[2] - f[1]) / h
    };
}

// Конечно-разностные производные на сетке 3x3
inline std::tuple<triplet, triplet> derivatives(const table_3x3& f, double h1, double h2) {
    triplet df_d1 = {
        (f[1][1] - f[0][1]) / h1,
        (f[2][1] - f[0][1]) / (2 * h1),
        (f[2][1] - f[1][1]) / h1
    };
    triplet df_d2 = {
        (f[1][1] - f[1][0]) / h2,
        (f[1][2] - f[1][0]) / (2 * h2),
        (f[1][2] - f[1][1]) / h2
    };
    return {df_d1, df_d2};
}

// Объемные доли на шаблоне 3 x 3
class Stencil2D {
    table_3x3 arr;

public:
    Stencil2D(EuCell& cell, Storable<double> a) {
        for (int i: {-1, 0, 1}) {
            for (int j: {-1, 0, 1}) {
                C(i, j) = cell.neib(i, j)[a];
            }
        }
    }

    double& C(int i, int j)       { return arr[i + 1][j + 1]; }
    double  C(int i, int j) const { return arr[i + 1][j + 1]; }

    double& operator()(int i, int j)       { return C(i, j); }
    double  operator()(int i, int j) const { return C(i, j); }

    Vector3d Youngs(double hx, double hy) const {
        Vector3d norm = Vector3d::Zero();
        norm.x() = ((C(-1, +1) + 2 * C(-1, 0) + C(-1, -1)) - (C(+1, +1) + 2 * C(+1, 0) + C(+1, -1))) / hx;
        norm.y() = ((C(+1, -1) + 2 * C(0, -1) + C(-1, -1)) - (C(+1, +1) + 2 * C(0, +1) + C(-1, +1))) / hy;
        return norm.normalized();
    }

    Vector3d ELVIRA(double hx, double hy) const {
        // Функции строк x(y)
        triplet func_x = {
            (C(-1, -1) + C(0, -1) + C(+1, -1)) * hx,
            (C(-1,  0) + C(0,  0) + C(+1,  0)) * hx,
            (C(-1, +1) + C(0, +1) + C(+1, +1)) * hx
        };
        triplet dx_dy = derivatives(func_x, hy);

        // Функции столбцов y(x)
        triplet func_y = {
            (C(-1, -1) + C(-1, 0) + C(-1, 1)) * hy,
            (C( 0, -1) + C( 0, 0) + C( 0, 1)) * hy,
            (C(+1, -1) + C(+1, 0) + C(+1, 1)) * hy
        };
        triplet dy_dx = derivatives(func_y, hx);

        // Нормали для тестов
        Vector3d YN = Youngs(hx, hy);
        double sgn_nx = math::sign(YN.x());
        double sgn_ny = math::sign(YN.y());
        std::array ns = {
            Vector3d{-dy_dx[0], sgn_ny, 0.0}.normalized(),
            Vector3d{-dy_dx[1], sgn_ny, 0.0}.normalized(),
            Vector3d{-dy_dx[2], sgn_ny, 0.0}.normalized(),
            Vector3d{sgn_nx, -dx_dy[0], 0.0}.normalized(),
            Vector3d{sgn_nx, -dx_dy[1], 0.0}.normalized(),
            Vector3d{sgn_nx, -dx_dy[2], 0.0}.normalized(),
        };

        std::array<double, 6> errors;
        for (int i = 0; i < 6; ++i) {
            Vector3d n = ns[i];
            double p = quad_find_section(C(0, 0), n, hx, hy);

            errors[i] = 0.0;
            errors[i] += std::pow(C(-1,  0) - quad_volume_fraction(p + hx * n.x(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1,  0) - quad_volume_fraction(p - hx * n.x(), n, hx, hy), 2);
            errors[i] += std::pow(C( 0, -1) - quad_volume_fraction(p + hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C( 0, +1) - quad_volume_fraction(p - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(-1, -1) - quad_volume_fraction(p + hx * n.x() + hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1, +1) - quad_volume_fraction(p - hx * n.x() - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(-1, +1) - quad_volume_fraction(p + hx * n.x() - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1, -1) - quad_volume_fraction(p - hx * n.x() + hy * n.y(), n, hx, hy), 2);
        }

        int err_min = std::min_element(errors.begin(), errors.end()) - errors.begin();
        return ns[err_min];
    }
};

// Объемные доли на шаблоне 3 x 3 x 3
class Stencil3D {
    table_3x3x3 arr;

public:
    Stencil3D(EuCell& cell, Storable<double> a) {
        for (int i: {-1, 0, 1}) {
            for (int j: {-1, 0, 1}) {
                for (int k: {-1, 0, 1}) {
                    C(i, j, k) = cell.neib(i, j, k)[a];
                }
            }
        }
    }

    double& C(int i, int j, int k)       { return arr[i + 1][j + 1][k + 1]; }
    double  C(int i, int j, int k) const { return arr[i + 1][j + 1][k + 1]; }

    double& operator()(int i, int j, int k)       { return C(i, j, k); }
    double  operator()(int i, int j, int k) const { return C(i, j, k); }

    Vector3d Youngs(double hx, double hy, double hz) const {
        Vector3d norm = Vector3d::Zero();
        norm.x()  = C(+1,-1,-1) + C(+1,-1,1) + C(+1,1,-1) + C(+1,1,1) + 2*(C(+1,0,1) + C(+1,0,-1) + C(+1,1,0) + C(+1,-1,0)) + 4*C(+1,0,0);
        norm.x() -= C(-1,-1,-1) + C(-1,-1,1) + C(-1,1,-1) + C(-1,1,1) + 2*(C(-1,0,1) + C(-1,0,-1) + C(-1,1,0) + C(-1,-1,0)) + 4*C(-1,0,0);
        norm.x() /= hx;
        norm.y()  = C(-1,+1,-1) + C(-1,+1,1) + C(1,+1,-1) + C(1,+1,1) + 2*(C(0,+1,1) + C(0,+1,-1) + C(1,+1,0) + C(-1,+1,0)) + 4*C(0,+1,0);
        norm.y() -= C(-1,-1,-1) + C(-1,-1,1) + C(1,-1,-1) + C(1,-1,1) + 2*(C(0,-1,1) + C(0,-1,-1) + C(1,-1,0) + C(-1,-1,0)) + 4*C(0,-1,0);
        norm.y() /= hy;
        norm.z()  = C(-1,-1,+1) + C(-1,1,+1) + C(1,-1,+1) + C(1,1,+1) + 2*(C(0,1,+1) + C(0,-1,+1) + C(1,0,+1) + C(-1,0,+1)) + 4*C(0,0,+1);
        norm.z() -= C(-1,-1,-1) + C(-1,1,-1) + C(1,-1,-1) + C(1,1,-1) + 2*(C(0,1,-1) + C(0,-1,-1) + C(1,0,-1) + C(-1,0,-1)) + 4*C(0,0,-1);
        norm.z() /= hz;
        return -norm.normalized();
    }

    Vector3d ELVIRA(double hx, double hy, double hz) const {
        Vector3d YN = Youngs(hx, hy, hz);

        // Height-functions: x(y, z), y(z, x), z(x, y)
        table_3x3 func_x, func_y, func_z;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                func_x[i][j] = func_y[i][j] = func_z[i][j] = 0.0;
                for (int k = 0; k < 3; ++k) {
                    func_x[i][j] += arr[k][i][j] * hx;
                    func_y[i][j] += arr[j][k][i] * hy;
                    func_z[i][j] += arr[i][j][k] * hz;
                }
            }
        }

        // Производные функций
        auto[dx_dy, dx_dz] = derivatives(func_x, hy, hz);
        auto[dy_dz, dy_dx] = derivatives(func_y, hz, hx);
        auto[dz_dx, dz_dy] = derivatives(func_z, hx, hy);

        double min_error = 1e300;
        Vector3d norm = Vector3d::Zero();

        Vector3d n_sgn = {math::sign(YN.x()), math::sign(YN.y()), math::sign(YN.z())};

        for (int k1 = 0; k1 < 3; ++k1) {
            for (int k2 = 0; k2 < 3; ++k2) {
                Vector3d test_n_x = Vector3d{n_sgn.x(), -dx_dy[k1], -dx_dz[k2]}.normalized();
                Vector3d test_n_y = Vector3d{-dy_dx[k1], n_sgn.y(), -dy_dz[k2]}.normalized();
                Vector3d test_n_z = Vector3d{-dz_dx[k1], -dz_dy[k2], n_sgn.z()}.normalized();

                double test_p_x = cube_find_section(C(0,0,0), test_n_x, hx, hy, hz);
                double test_p_y = cube_find_section(C(0,0,0), test_n_y, hx, hy, hz);
                double test_p_z = cube_find_section(C(0,0,0), test_n_z, hx, hy, hz);

                // Проверить сечения 9 столбцов
                double err_x{0.0}, err_y{0.0}, err_z{0.0};
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        Vector3d r_x = {0.0, (i - 1) * hy, (j - 1) * hz};
                        Vector3d r_y = {(j - 1) * hx, 0.0, (i - 1) * hz};
                        Vector3d r_z = {(i - 1) * hx, (j - 1) * hy, 0.0};

                        double vf_x = cube_volume_fraction(test_p_x - r_x.dot(test_n_x), test_n_x, 3*hx, hy, hz);
                        double vf_y = cube_volume_fraction(test_p_y - r_y.dot(test_n_y), test_n_y, hx, 3*hy, hz);
                        double vf_z = cube_volume_fraction(test_p_z - r_z.dot(test_n_z), test_n_z, hx, hy, 3*hz);

                        err_x += std::pow(func_x[i][j]/(3*hx) - vf_x, 2);
                        err_y += std::pow(func_y[i][j]/(3*hy) - vf_y, 2);
                        err_z += std::pow(func_z[i][j]/(3*hz) - vf_z, 2);
                    }
                }

                if (err_x < min_error) {
                    min_error = err_x;
                    norm = test_n_x;
                }
                if (err_y < min_error) {
                    min_error = err_y;
                    norm = test_n_y;
                }
                if (err_z < min_error) {
                    min_error = err_z;
                    norm = test_n_z;
                }
            }
        }
        return norm;
    }
};
