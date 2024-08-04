#pragma once

#include <functional>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math::gradient {

template<class EigenType>
void check_matrix(const EigenType &m, const std::string &name) {
    if (m.hasNaN()) {
        std::cerr << "Bad matrix " << name << " in compute grad:\n" << m << '\n';
        throw std::runtime_error("Bad Matrix");
    }
}

template<class T>
std::array<T, 3> gauss(Cell &cell, const std::function<T(Cell &)> &get_state) {
    std::array<T, 3> grad;
    T zc = get_state(cell);

    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn = get_state(neib);

        Vector3d S = face.normal() * face.area();

        double d1 = std::abs(1.0 / (face.center() - cell.center()).dot(face.normal()));
        double d2 = std::abs(1.0 / (face.center() - neib.center()).dot(face.normal()));
        double a1 = d1 / (d1 + d2);
        double a2 = d2 / (d1 + d2);

        // Значение на грани
        T zf = a1 * zc.vec() + a2 * zn.vec();

        grad[0].arr() += zf.arr() * S.x();
        grad[1].arr() += zf.arr() * S.y();
        grad[2].arr() += zf.arr() * S.z();
    }

    for (auto &z: grad) {
        z.arr() /= cell.volume();
    }

    return grad;
}

template<class T>
std::array<T, 3> LSM(Cell &cell,
        const std::function<T(Cell &)> &get_state,
        const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    using Array = ei_arr<T>;

    T zc = get_state(cell);

    Array Fx = Array::Zero();
    Array Fy = Array::Zero();
    Array Fz = Array::Zero();

    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        T zn;
        if (!face.is_boundary()) {
            auto neib = face.neib();
            zn = get_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2.0 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr * dr.transpose();

        Array dF = zn.arr() - zc.arr();

        Fx += weight * dF * dr.x();
        Fy += weight * dF * dr.y();
        Fz += weight * dF * dr.z();
    }

    if (std::abs(A(1, 1)) < 1e-14) {
        A(1, 1) = 1.0;
    }
    if (std::abs(A(2, 2)) < 1e-14) {
        A(2, 2) = 1.0;
    }

    check_matrix(A, "A");
    check_matrix(Fx, "Fx");
    check_matrix(Fy, "Fy");
    check_matrix(Fz, "Fz");

    Matrix3d B = A.inverse();
    check_matrix(B, "B");

    std::array<T, 3> grad;
    grad[0] = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    grad[1] = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    grad[2] = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;

    return grad;
}

template<class T>
std::array<T, 3> limiting(Cell &cell,
        const std::array<T, 3> &grad,
        const std::function<T(Cell &)> &get_state,
        const std::function<T(const T &, const Vector3d &, Boundary)> &boundary_value) {
    using Array = ei_arr<T>;

    T zc = get_state(cell);

    const double epsilon = 1.0e-14;
    constexpr int rows = 1;
    constexpr int cols = T::size();
    using Array = ei_arr<T>;

    Array Fx = Array::Zero();
    Array Fy = Array::Zero();
    Array Fz = Array::Zero();

    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        T zn;
        if (!face.is_boundary()) {
            auto neib = face.neib();
            zn = get_state(neib);
        } else {
            zn = boundary_value(zc, face.normal(), face.flag());
        }

        Vector3d dr = 2.0 * (face.center() - cell.center());
        double weight = face.area() / dr.squaredNorm();

        A += weight * dr * dr.transpose();

        Array dF = zn.vec() - zc.vec();
        Array theta = 2 * (grad[0].vec() * dr.x() + grad[1].vec() * dr.y() + grad[2].vec() * dr.z());
        Array dF_lim = Array::Zero();
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) {
                theta(i, j) = limiters::MC(theta(i, j) / (dF(i, j) + epsilon) - 1);
                dF_lim(i, j) = dF(i, j) * theta(i, j);
            }

        Fx += weight * dF_lim * dr.x();
        Fy += weight * dF_lim * dr.y();
        Fz += weight * dF_lim * dr.z();
    }

    if (std::abs(A(1, 1)) < 1e-14) {
        A(1, 1) = 1.0;
    }
    if (std::abs(A(2, 2)) < 1e-14) {
        A(2, 2) = 1.0;
    }

    check_matrix(A, "A");
    check_matrix(Fx, "Fx");
    check_matrix(Fy, "Fy");
    check_matrix(Fz, "Fz");

    Matrix3d B = A.inverse();
    check_matrix(B, "B");

    std::array<T, 3> lim_grad;
    lim_grad[0] = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    lim_grad[1] = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    lim_grad[2] = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;

    return lim_grad;
}

} // namespace zephyr::math::gradient