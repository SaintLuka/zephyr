#include <zephyr/utils/numpy.h>

namespace zephyr::np {

template <>
Array<np::real> zeros<np::real>(size_t N) {
    return std::vector<np::real>(N, 0.0);
}

template <>
Array<np::vec3> zeros<np::vec3>(size_t N) {
    return std::vector<np::vec3>(N, vec3::Zero());
}

Array2d zeros(size_t N, size_t M) {
    return std::vector<std::vector<real>>(N, std::vector<real>(M, 0.0));
}

Array1d zeros_like(const Array1d& arr) {
    return zeros(arr.size());
}

Array2d zeros_like(const Array2d& arr) {
    if (arr.empty()) { return {}; }
    return zeros(arr.size(), arr[0].size());
}

Array<real> linspace(real t1, real t2, size_t N) {
    Array1d res(N);
    real h = (t2 - t1) / (N - 1.0);
    for (size_t i = 0; i < N; ++i) {
        res[i] = t1 + i * h;
    }
    return res;
}

std::tuple<Array2d,
        Array2d> meshgrid(
        const Array1d& x,
        const Array1d& y) {

    std::vector X(x.size(), Array1d(y.size()));
    std::vector Y(x.size(), Array1d(y.size()));

    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < y.size(); ++j) {
            X[i][j] = x[i];
            Y[i][j] = y[j];
        }
    }
    return {X, Y};
}

Array1v linspace(
        const vec3& v1, const vec3& v2, size_t N) {
    Array1v res(N);
    vec3 tau = (v2 - v1) / (N - 1.0);
    for (size_t i = 0; i < N; ++i) {
        res[i] = v1 + tau * i;
    }
    return res;
}

Array1d get_x(const Array1v& arr) {
    Array1d res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].x();
    }
    return res;
}

Array1d get_y(const Array1v& arr) {
    Array1d res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].y();
    }
    return res;
}

Array1d get_z(const Array1v& arr) {
    Array1d res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].z();
    }
    return res;
}

Array1v zip(const Array1d& xs,
                                 const Array1d& ys) {
    if (xs.size() != ys.size()) {
        throw std::runtime_error("zip error: different array size");
    }

    Array1v res(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        res[i] = {xs[i], ys[i], 0.0};
    }
    return res;
}

Array1v zip(const Array1d& xs,
                                 const Array1d& ys,
                                 const Array1d& zs) {
    if (xs.size() != ys.size() || xs.size() != zs.size()) {
        throw std::runtime_error("zip error: different array size");
    }

    Array1v res(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        res[i] = {xs[i], ys[i], zs[i]};
    }
    return res;
}

} // namespace zephyr::np


