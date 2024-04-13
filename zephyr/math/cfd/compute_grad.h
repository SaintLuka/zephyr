#pragma once

#include <functional>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

template<int row, int col>
void check_matrix(const Matrix<double, row, col> &m, const std::string &name) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (std::isnan(m(i, j)) || std::isinf(m(i, j))) {
                std::cerr << "Bad matrix " << name << " in compute grad:\n" << m << '\n';
                throw std::runtime_error("Bad Matrix");
            }
        }
    }
}

template<class T>
std::array<T, 3> compute_grad_gauss(Cell &cell, const std::function<T(Cell &)> &to_state) {
    std::array<T, 3> z_xyz;
    T zc = to_state(cell);

    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn = to_state(neib);

        Eigen::Vector3d S = face.normal() * face.area();
        double d1 = abs(1 / (face.center() - cell.center()).dot(face.normal()));
        double d2 = abs(1 / (face.center() - neib.center()).dot(face.normal()));
        double a1 = d1 / (d1 + d2), a2 = d2 / (d1 + d2);
        T state = a1 * zc.vec() + a2 * zn.vec();

        z_xyz[0].vec() += state.vec() * S.x();
        z_xyz[1].vec() += state.vec() * S.y();
        z_xyz[2].vec() += state.vec() * S.z();
    }

    for (auto &z: z_xyz)
        z.vec() /= cell.volume();

    return z_xyz;
}

template<class T>
std::array<T, 3> compute_gradient_LSM(Cell &cell,
                                      const std::function<T(Cell &)> &to_state,
                                      const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    T zc = to_state(cell);
    constexpr int rows = 1;
    constexpr int cols = T::size();
    using StateVec = Matrix<double, rows, cols>;

    StateVec Fx = StateVec::Zero(), Fy = StateVec::Zero(), Fz = StateVec::Zero();
    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn(zc);

        if (!face.is_boundary()) {
            zn = to_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr * dr.transpose();

        StateVec dF = zn.vec() - zc.vec();

        Fx += weight * dF * dr.x();
        Fy += weight * dF * dr.y();
        Fz += weight * dF * dr.z();
    }

    if (Fz.squaredNorm() < 1e-10 && abs(A(2, 2)) < 1e-10) {
        A(2, 2) = 1;
    }
    check_matrix(A, "A");
    check_matrix(Fx, "Fx");
    check_matrix(Fy, "Fy");
    check_matrix(Fz, "Fz");

    Matrix3d B = A.inverse();
    check_matrix(B, "B");

    std::array<T, 3> z_xyz;
    z_xyz[0] = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    z_xyz[1] = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    z_xyz[2] = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;

    return z_xyz;
}

template<class T>
std::array<T, 3> gradient_limiting(Cell &cell, const std::array<T, 3> &grad,
                                   const std::function<T(Cell &)> &to_state,
                                   const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    T zc = to_state(cell);

    const double epsilon = std::numeric_limits<double>::min();
    constexpr int rows = 1;
    constexpr int cols = T::size();
    using StateVec = Matrix<double, rows, cols>;

    StateVec Fx = StateVec::Zero(), Fy = StateVec::Zero(), Fz = StateVec::Zero();
    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn(zc);

        if (!face.is_boundary()) {
            zn = to_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr * dr.transpose();

        StateVec dF = zn.vec() - zc.vec();
        StateVec theta = 2 * (grad[0].vec() * dr.x() + grad[1].vec() * dr.y() + grad[2].vec() * dr.z());
        StateVec dF_lim = StateVec::Zero();
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) {
                theta(i, j) = limiters::MC(theta(i, j) / (dF(i, j) + epsilon) - 1);
                dF_lim(i, j) = dF(i, j) * theta(i, j);
            }

        Fx += weight * dF_lim * dr.x();
        Fy += weight * dF_lim * dr.y();
        Fz += weight * dF_lim * dr.z();
    }

    if (Fz.squaredNorm() < 1e-10 && abs(A(2, 2)) < 1e-10) {
        A(2, 2) = 1;
    }
    check_matrix(A, "A");
    check_matrix(Fx, "Fx");
    check_matrix(Fy, "Fy");
    check_matrix(Fz, "Fz");

    Matrix3d B = A.inverse();
    check_matrix(B, "B");

    std::array<T, 3> z_xyz;
    z_xyz[0] = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    z_xyz[1] = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    z_xyz[2] = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;

    return z_xyz;
}

template<class T>
T compute_gradient_LSM_1D(Cell &cell,
                          const std::function<T(Cell &)> &to_state,
                          const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    T zc = to_state(cell);
    constexpr int rows = 1;
    constexpr int cols = T::size();
    using StateVec = Matrix<double, rows, cols>;

    StateVec Fx = StateVec::Zero();
    double A = 0;
    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn(zc);

        if (!face.is_boundary()) {
            zn = to_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr.x() * dr.x();

        StateVec dF = zn.vec() - zc.vec();

        Fx += weight * dF * dr.x();
    }

    check_matrix(Fx, "Fx");

    return Fx / A;
}

template<class T>
T gradient_limiting_1D(Cell &cell, const T &grad,
                       const std::function<T(Cell &)> &to_state,
                       const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    T zc = to_state(cell);

    const double epsilon = std::numeric_limits<double>::min();
    constexpr int rows = 1;
    constexpr int cols = T::size();
    using StateVec = Matrix<double, rows, cols>;

    StateVec Fx = StateVec::Zero();
    double A = 0;
    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn(zc);

        if (!face.is_boundary()) {
            zn = to_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr.x() * dr.x();

        StateVec dF = zn.vec() - zc.vec();
        StateVec theta = 2 * grad.vec() * dr.x();
        StateVec dF_lim = StateVec::Zero();
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) {
                theta(i, j) = limiters::MC(theta(i, j) / (dF(i, j) + epsilon) - 1);
                dF_lim(i, j) = dF(i, j) * theta(i, j);
            }

        Fx += weight * dF_lim * dr.x();
    }

    check_matrix(Fx, "Fx");

    return Fx / A;
}


/*

void compute_gradient_LSM(ICell::Ref cell) {
    auto zc = get_component(cell, Part::INIT);
    svector dzx;
    svector dzy;
    svector dzz;

    Point cell_c = cell->center();

    svector Fx = svector::zero();
    svector Fy = svector::zero();
    svector Fz = svector::zero();

    matrix3D A(0.0);

    //swap over neibhbors
    for (auto &face: cell->faces()) {

        auto flag = face->flag();

        auto neib = face->neibhbor(cell);
        auto normal = face->normal(neib);
        Point face_c = face->center();

        Point neib_c = neib->center();

        auto zn = get_component(neib, Part::INIT);

        if (face->is_boundary()) {
            Point face_c = face->center();
            neib_c = get_symm_point(cell_c, face_c, normal);
            zn = boundary_value(zc, normal, neib_c, flag);
        }

        //ORIGINAL; by default
        Point dr = neib_c - cell_c;
        double weight = 1.0 / length_squared(dr);

        //MRV; by default
        if (is_AMR_slope_corr() && is_AMR_slope_corr_MRVtype()) {
            weight = face->area() / length_squared(dr);
        }

        A += weight * tensor_dot(dr, dr);

        Fx += weight * (zn - zc) * dr.x;
        Fy += weight * (zn - zc) * dr.y;
        Fz += weight * (zn - zc) * dr.z;
    }

    // Вычисляем обычные производные
    if (is_2D()) {
        A.zz = 1.0;
    }
    matrix3D B = A.inverse();
    dzx = B.xx * Fx + B.xy * Fy + B.xz * Fz;
    dzy = B.yx * Fx + B.yy * Fy + B.yz * Fz;
    dzz = B.zx * Fx + B.zy * Fy + B.zz * Fz;

    if (dzx.hasNaN() || dzy.hasNaN() || dzz.hasNaN()) {
        //throw runtime_error("compute_gradient(cell) error: d/dx, d/dy or d/dz is NaN");
        dzx = svector::zero();
        dzy = dzx;
        dzz = dzx;
    }

    set_component(cell, dzx, Part::D_DX);
    set_component(cell, dzy, Part::D_DY);
    set_component(cell, dzz, Part::D_DZ);
}

void gradient_limiting(ICell::Ref cell) {

    const double epsilon = std::numeric_limits<double>::min();

    auto zc = get_component(cell, Part::INIT);
    auto dzx = get_component(cell, Part::D_DX);
    auto dzy = get_component(cell, Part::D_DY);
    auto dzz = get_component(cell, Part::D_DZ);

    Point cell_c = cell->center();

    auto Fx = svector::zero();
    auto Fy = svector::zero();
    auto Fz = svector::zero();
    matrix3D A(0.0);

    //swap over neibs
    for (auto &face: cell->faces()) {
        auto flag = face->flag();
        auto neib = face->neibhbor(cell);
        auto normal = face->normal(neib);
        Point face_c = face->center();

        Point neib_c = neib->center();
        auto zn = get_component(neib, Part::INIT);

        if (face->is_boundary()) {
            Point face_c = face->center();
            neib_c = get_symm_point(cell_c, face_c, normal);
            zn = boundary_value(zc, normal, neib_c, flag);
        }

        //ORIGINAL; by default
        Point dr = neib_c - cell_c;
        double weight = 1.0 / length_squared(dr);

        //MRV; by default
        if (is_AMR_slope_corr() && is_AMR_slope_corr_MRVtype()) {
            weight = face->area() / length_squared(dr);
        }

        A += weight * tensor_dot(dr, dr);

        //zn - уже может быnm
        svector dF = zn - zc;
        for (size_t i = 0; i < dF.size(); ++i) {

            if (std::abs(dF[i]) <= 1.0e-12 * std::max(std::abs(zn[i]), std::abs(zc[i]))) {
                dF[i] = 0.0;
            }
        }

        svector theta = 2.0 * (dzx * dr.x + dzy * dr.y + dzz * dr.z) / (dF + epsilon) - 1.0;
        svector lim = limiter(theta);
        svector dF_lim = lim * dF;

        Fx += weight * dF_lim * dr.x;
        Fy += weight * dF_lim * dr.y;
        Fz += weight * dF_lim * dr.z;
    }

    //Лимитированные производные
    if (is_2D()) {
        A.zz = 1.0;
    }
    auto B = A.inverse();
    dzx = B.xx * Fx + B.xy * Fy + B.xz * Fz;
    dzy = B.yx * Fx + B.yy * Fy + B.yz * Fz;
    dzz = B.zx * Fx + B.zy * Fy + B.zz * Fz;

    if (dzx.hasNaN() || dzy.hasNaN() || dzz.hasNaN()) {
        //throw runtime_error("compute_gradient(cell) error: d/dx, d/dy or d/dz is NaN");
        dzx = svector::zero();
        dzy = dzx;
        dzz = dzx;
    }

    set_component(cell, dzx, Part::D_DX);
    set_component(cell, dzy, Part::D_DY);
    set_component(cell, dzz, Part::D_DZ);
}

*/

}