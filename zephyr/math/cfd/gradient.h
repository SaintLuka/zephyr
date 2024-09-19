#pragma once

#include <functional>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math::gradient {

/// @brief Получить вектор состояния типа T из ячейки
template <class T>
using GetState = std::function<T(Cell &)>;

/// @brief Получить вектор состояния типа T через грань
/// с граничным условием
template <class T>
using GetBoundary = std::function<T(
        const T &,         // Вектор состояния в ячейке
        const Vector3d &,  // Внешняя нормаль к грани
        Boundary           // Флаг граничного условия
)>;

/// @brief Простенькая структура для хранения градиента
template <class State>
struct Grad {
    /// @brief Векторная версия типа
    using Array = ei_arr<State>;

    /// @brief Компоненты градиента
    State x, y, z;

    Grad() {
        (Array &) x = Array::Zero();
        (Array &) y = Array::Zero();
        (Array &) z = Array::Zero();
    }

    // Приведение компонент к Eigen массивам
    inline Array &x_arr() { return (Array &) x; }

    inline Array &y_arr() { return (Array &) y; }

    inline Array &z_arr() { return (Array &) z; }

    inline const Array &x_arr() const { return (Array &) x; }

    inline const Array &y_arr() const { return (Array &) y; }

    inline const Array &z_arr() const { return (Array &) z; }
};

/// @brief Расчет градиента методом Гаусса
template<class State>
Grad<State> gauss(Cell &cell,
                  const GetState<State> &get_state,
                  const GetBoundary<State> &boundary_value) {
    using Array = ei_arr<State>;

    Grad<State> grad;

    State zc = get_state(cell);
    Vector3d cell_c = cell.center();
    for (auto &face: cell.faces()) {
        const Vector3d &normal = face.normal();
        const Vector3d &face_c = face.center();

        State zn;
        Vector3d neib_c;
        if (!face.is_boundary()) {
            auto neib = face.neib();
            zn = get_state(neib);
            neib_c = face.neib().center();
        } else {
            zn = boundary_value(zc, normal, face.flag());
            neib_c = face.symm_point(cell_c);
        }

        Vector3d S = face.area() * normal;

        double d1 = std::abs(1.0 / (face_c - cell_c).dot(normal));
        double d2 = std::abs(1.0 / (face_c - neib_c).dot(normal));
        double a1 = d1 / (d1 + d2);
        double a2 = d2 / (d1 + d2);

        // Значение на грани
        Array zf = a1 * (Array &) zc + a2 * (Array &) zn;

        grad.x_arr() += zf * S.x();
        grad.y_arr() += zf * S.y();
        grad.z_arr() += zf * S.z();
    }

    double V = cell.volume();
    grad.x_arr() /= V;
    grad.y_arr() /= V;
    grad.z_arr() /= V;

    return grad;
}

/// @brief Метод наименьших квадратов (классическая версия)
template<class State>
Grad<State> LSM_orig(Cell &cell,
                const GetState<State> &get_state,
                const GetBoundary<State> &boundary_value) {
    using Array = ei_vec<State>;

    State zc = get_state(cell);
    Vector3d cell_c = cell.center();

    Array Fx = Array::Zero();
    Array Fy = Array::Zero();
    Array Fz = Array::Zero();
    Matrix3d A = Matrix3d::Zero();

    for (auto &face: cell.faces()) {
        const Vector3d &normal = face.normal();

        State zn;
        Vector3d dr;
        if (!face.is_boundary()) {
            auto neib = face.neib();
            zn = get_state(neib);
            dr = neib.center() - cell_c;
        } else {
            zn = boundary_value(zc, normal, face.flag());
            dr = face.symm_point(cell_c) - cell_c;
        }

        double w = 1.0 / dr.squaredNorm();

        A += w * dr * dr.transpose();

        Array dF = (Array &) zn - (Array &) zc;

        Fx += w * dF * dr.x();
        Fy += w * dF * dr.y();
        Fz += w * dF * dr.z();
    }

    // Матрица A безразмерная, коэффициенты порядка единиц
    if (std::abs(A(2, 2)) < 1e-12) {
        A(2, 2) = 1.0; // Случай 2D
        if (std::abs(A(1, 1)) < 1e-12) {
            A(1, 1) = 1.0; // Случай 1D
        }
    }

    Matrix3d B = A.inverse();

    Grad<State> grad;
    grad.x_arr() = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    grad.y_arr() = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    grad.z_arr() = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;
    return grad;
}

/// @brief Ограничитель градиента (классическая версия)
template<class State>
Grad<State> limiting_orig(Cell &cell, const Limiter& limiter,
                     const Grad<State> &grad,
                     const GetState<State> &get_state,
                     const GetBoundary<State> &boundary_value) {
    using Array = ei_arr<State>;

    State zc = get_state(cell);
    Vector3d cell_c = cell.center();

    Array Fx = Array::Zero();
    Array Fy = Array::Zero();
    Array Fz = Array::Zero();

    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        const Vector3d &normal = face.normal();

        State zn;
        Vector3d dr;
        if (!face.is_boundary()) {
            auto neib = face.neib();
            zn = get_state(neib);
            dr = neib.center() - cell_c;
        } else {
            zn = boundary_value(zc, normal, face.flag());
            dr = face.symm_point(cell_c) - cell_c;
        }

        double w = 1.0 / dr.squaredNorm();

        A += w * dr * dr.transpose();

        Array dF_p = (Array &) zn - (Array &) zc;
        Array du_n = grad.x_arr() * dr.x() + grad.y_arr() * dr.y() + grad.z_arr() * dr.z();
        Array dF_m = 2.0 * du_n - dF_p;

        Array dF_lim = limiter(dF_m, dF_p) * dF_p;

        Fx += w * dF_lim * dr.x();
        Fy += w * dF_lim * dr.y();
        Fz += w * dF_lim * dr.z();
    }

    // Матрица A безразмерная, коэффициенты порядка единиц
    if (std::abs(A(2, 2)) < 1e-12) {
        A(2, 2) = 1.0; // Случай 2D
        if (std::abs(A(1, 1)) < 1e-12) {
            A(1, 1) = 1.0; // Случай 1D
        }
    }

    Matrix3d B = A.inverse();

    Grad<State> lim_grad;
    lim_grad.x_arr() = B(0, 0) * Fx + B(0, 1) * Fy + B(0, 2) * Fz;
    lim_grad.y_arr() = B(1, 0) * Fx + B(1, 1) * Fy + B(1, 2) * Fz;
    lim_grad.z_arr() = B(2, 0) * Fx + B(2, 1) * Fy + B(2, 2) * Fz;
    return lim_grad;
}

} // namespace zephyr::math::gradient