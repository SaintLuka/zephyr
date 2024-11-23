#include <iostream>

#include <zephyr/math/sle/tridiagonal.h>

namespace zephyr::math::tridiagonal {

inline bool diff_size(const array& x, const array& y) {
    return x.size() != y.size();
}

inline bool diff_size(const array& x, const array& y, const array& z, const array& w) {
    return diff_size(x, y) || diff_size(x, z) || diff_size(x, w);
}

array solve(const array& A, const array& B, const array& C, const array& F) {
    if (diff_size(A, B, C, F)) {
        throw std::runtime_error("tridiagonal::solve error: different size");
    }

    int n = A.size();

    // Создает два дополнительных массива
    array c(n), x(n);

    c[0] = C[0] / B[0];
    x[0] = F[0] / B[0];

    for (int i = 1; i < n; ++i) {
        double div = 1.0 / (B[i] - A[i] * c[i - 1]);
        c[i] = C[i] * div;
        x[i] = (F[i] - A[i] * x[i - 1]) * div;
    }

    for (int i = n - 2; i >= 0; --i) {
        x[i] -= c[i] * x[i + 1];
    }

    return x;
}

void solve_l(const array& A, const array& B, const array& C, array& F) {
    array C2(C);
    solve_m(A, B, C2, F);
}

void solve_m(const array& A, const array& B, array& C, array& F) {
    if (diff_size(A, B, C, F)) {
        throw std::runtime_error("tridiagonal::solve error: different size");
    }

    int n = A.size();

    C[0] /= B[0];
    F[0] /= B[0];
    for (int i = 1; i < n; ++i) {
        double div = 1.0 / (B[i] - A[i] * C[i - 1]);
        C[i] *= div;
        F[i] = (F[i] - A[i] * F[i - 1]) * div;
    }

    for (int i = n - 2; i >= 0; --i) {
        F[i] -= C[i] * F[i + 1];
    }
}

array solve_cyclic(const array& A, const array& B, const array& C, const array& F) {
    if (diff_size(A, B, C, F)) {
        throw std::runtime_error("tridiagonal::solve_cyclic error: different size");
    }

    int n = A.size();

    double vn = -A[0] / B[0];

    array c(n), q(n), x(n);

    c[0] = +0.5 * C[0] / B[0];
    q[0] = -0.5;
    x[0] = +0.5 * F[0] / B[0];

    for (int i = 1; i < n - 1; ++i) {
        double div = 1.0 / (B[i] - A[i] * c[i - 1]);
        c[i] = C[i] * div;
        q[i] = -A[i] * q[i - 1] * div;
        x[i] = (F[i] - A[i] * x[i - 1]) * div;
    }

    double xi = 1.0 / (B[n - 1] - C[n - 1] * vn - A[n - 1] * c[n - 2]);

    q[n - 1] = (C[n - 1] - A[n - 1] * q[n - 2]) * xi;
    x[n - 1] = (F[n - 1] - A[n - 1] * x[n - 2]) * xi;

    for (int i = n - 2; i >= 0; --i) {
        q[i] -= c[i] * q[i + 1];
        x[i] -= c[i] * x[i + 1];
    }

    double fact = (x[0] + vn * x[n - 1]) / (1.0 + q[0] + vn * q[n - 1]);

    for (int i = 0; i < n; ++i) {
        x[i] -= q[i] * fact;
    }

    return x;
}

void solve_cyclic_l(const array& A, const array& B, const array& C, array& F) {
    array B2(B);
    array C2(C);
    solve_cyclic_m(A, B2, C2, F);
}

/// Работет магическим образом, не пытайтесь понять
void solve_cyclic_m(const array& A, array& B, array& C, array& F) {
    if (diff_size(A, B, C, F)) {
        throw std::runtime_error("tridiagonal::solve_cyclic error: different size");
    }

    int n = A.size();

    double vn = -A[0] / B[0];

    C[0] /= 2.0 * B[0];
    F[0] /= 2.0 * B[0];
    B[0] = -0.5;

    for (int i = 1; i < n - 1; ++i) {
        double div = 1.0 / (B[i] - A[i] * C[i - 1]);
        C[i] *= div;
        B[i] = -A[i] * B[i - 1] * div;
        F[i] = (F[i] - A[i] * F[i - 1]) * div;
    }

    double xi = 1.0 / (B[n - 1] - C[n - 1] * vn - A[n - 1] * C[n - 2]);

    B[n - 1] = (C[n - 1] - A[n - 1] * B[n - 2]) * xi;
    F[n - 1] = (F[n - 1] - A[n - 1] * F[n - 2]) * xi;

    for (int i = n - 2; i >= 0; --i) {
        B[i] -= C[i] * B[i + 1];
        F[i] -= C[i] * F[i + 1];
    }

    double fact = (F[0] + vn * F[n - 1]) / (1.0 + B[0] + vn * B[n - 1]);

    for (int i = 0; i < n; ++i) {
        F[i] -= B[i] * fact;
    }
}

} // namespace zephyr::math::tridiagonal