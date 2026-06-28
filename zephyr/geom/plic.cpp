#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/polyhedron.h>
#include <zephyr/geom/sections.h>
#include <zephyr/math/funcs.h>
#include <zephyr/geom/plic.h>

namespace zephyr::geom {

using mesh::EuCell;

Plic::plane_t Plic::plane(EuCell& cell, int idx) const {
    return m_find_plane(cell, idx);
}


Plic::Plic() {
    m_find_plane = [](const EuCell& cell, int idx) ->plane_t {
        return {0.0, Vector3d::Zero()};
    };
}

namespace {

// Простая производная (формула Гаусса)
Vector3d grad_normal(const EuCell& cell, int idx, const Plic::get_fraction_t& get_vf) {
    double a0 = get_vf(cell, idx);
    Vector3d grad = Vector3d::Zero();
    for (auto face: cell.faces()) {
        double a_f = 0.5 * (a0 + get_vf(face.neib(), idx));
        grad += a_f * face.area_n();
    }
    return -grad.normalized();
}

// Двумерная формула CSIR
Vector3d csir2_normal(const EuCell& cell, int idx, const Plic::get_fraction_t& get_vf) {
    double a0 = get_vf(cell, idx);
    Vector3d grad = Vector3d::Zero();
    for (auto face: cell.faces()) {
        double a_f = face_fraction(a0, get_vf(face.neib(), idx));
        grad += a_f * face.area_n();
    }
    return -grad.normalized();
}

// Трёхмерная формула CSIR
Vector3d csir3_normal(const EuCell& cell, int idx, const Plic::get_fraction_t& get_vf) {
    // Собрать доли в соседних ячейках
    std::array<double, Side3D::count()> a_neib;
    for (auto side: Side3D::items()) {
        a_neib[side] = get_vf(cell.face(side).neib(), idx);
    }

    // Объемные доли на гранях
    double a0 = get_vf(cell, idx);
    auto a_f = face_fractions(a0, a_neib);

    // И обычный Гаусс
    Vector3d grad = Vector3d::Zero();
    for (auto face: cell.faces()) {
        grad += a_f[face.side()] * face.area_n();
    }
    return -grad.normalized();
}

double polygon_section(const EuCell& cell, int idx, const Vector3d& n, const Plic::get_fraction_t& get_vf) {
    Vector3d P = cell.polygon().find_section(n, get_vf(cell, idx)).p1;
    return (P - cell.center()).dot(n);
}

double polyhedron_section(const EuCell& cell, int idx, const Vector3d& n, const Plic::get_fraction_t& get_vf) {
    Vector3d P = cell.polyhedron().find_section(n, get_vf(cell, idx));
    return (P - cell.center()).dot(n);
}

double quad_section(const EuCell& cell, int idx, const Vector3d& n, const Plic::get_fraction_t& get_vf) {
    return quad_find_section(get_vf(cell, idx), n, cell.hx(), cell.hy());
}

double cube_section(const EuCell& cell, int idx, const Vector3d& n, const Plic::get_fraction_t& get_vf) {
    return cube_find_section(get_vf(cell, idx), n, cell.hx(), cell.hy(), cell.hz());
}

using triplet = std::array<double, 3>;
using table_3x3 = std::array<triplet, 3>;
using table_3x3x3 = std::array<table_3x3, 3>;

// Конечно-разностные производные на сетке из трёх узлов
triplet derivatives(const triplet& f, double h) {
    return {
        (f[1] - f[0]) / h,
        (f[2] - f[0]) / (2 * h),
        (f[2] - f[1]) / h
    };
}

// Конечно-разностные производные на сетке 3x3
std::tuple<triplet, triplet> derivatives(const table_3x3& f, double h1, double h2) {
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
    Stencil2D(const EuCell& cell, int idx, const Plic::get_fraction_t& get_vf) {
        for (int i: {-1, 0, 1}) {
            for (int j: {-1, 0, 1}) {
                C(i, j) = get_vf(cell.neib(i, j), idx);
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
    Stencil3D(const EuCell& cell, int idx, const Plic::get_fraction_t& get_vf) {
        for (int i: {-1, 0, 1}) {
            for (int j: {-1, 0, 1}) {
                for (int k: {-1, 0, 1}) {
                    C(i, j, k) = get_vf(cell.neib(i, j, k), idx);
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

} // anonymous namespace

Plic::Plic(int dim, bool cartesian, Type type, const get_fraction_t& get_vf) {
    if (type == GRAD) {
        if (!cartesian) {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = grad_normal(cell, idx, get_vf);
                    double p = polygon_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = grad_normal(cell, idx, get_vf);
                    double p = polyhedron_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #1");
            }
        }
        else {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = grad_normal(cell, idx, get_vf);
                    double p = quad_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = grad_normal(cell, idx, get_vf);
                    double p = cube_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #2");
            }
        }
    }
    else if (type == PnY) {
        if (cartesian) {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Stencil2D arr(cell, idx, get_vf);
                    Vector3d n = arr.Youngs(cell.hx(), cell.hy());
                    double p = quad_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Stencil3D arr(cell, idx, get_vf);
                    Vector3d n = arr.Youngs(cell.hx(), cell.hy(), cell.hz());
                    double p = cube_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #2");
            }
        }
        else {
            throw std::runtime_error("Plic::PnY is not implemented for general meshes");
        }
    }
    else if (type == ELVIRA) {
        if (cartesian) {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Stencil2D arr(cell, idx, get_vf);
                    Vector3d n = arr.ELVIRA(cell.hx(), cell.hy());
                    double p = quad_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Stencil3D arr(cell, idx, get_vf);
                    Vector3d n = arr.ELVIRA(cell.hx(), cell.hy(), cell.hz());
                    double p = cube_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #3");
            }
        }
        else {
            throw std::runtime_error("Plic::ELVIRA is not implemented for general meshes");
        }
    }
    else if (type == CSIR) {
        if (cartesian) {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = csir2_normal(cell, idx, get_vf);
                    double p = quad_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = csir3_normal(cell, idx, get_vf);
                    double p = cube_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #4");
            }
        }
        else {
            throw std::runtime_error("Plic::CSIR is not implemented for general meshes");
        }
    }
    else if (type == CSIR_2D) {
        if (cartesian) {
            if (dim == 2) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = csir2_normal(cell, idx, get_vf);
                    double p = quad_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else if (dim == 3) {
                m_find_plane = [get_vf](const EuCell& cell, int idx) -> plane_t {
                    Vector3d n = csir2_normal(cell, idx, get_vf);
                    double p = cube_section(cell, idx, n, get_vf);
                    return {p, n};
                };
            }
            else {
                throw std::runtime_error("Plic setup: bad dimension #5");
            }
        }
        else {
            throw std::runtime_error("Plic::CSIR is not implemented for general meshes");
        }
    }
}

} // namespace zephyr::geom