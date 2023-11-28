#pragma once

#include <zephyr/configuration.h>

#ifdef ZEPHYR_ENABLE_EIGEN

#define VECTORIZE(CustomVector) \
    static constexpr int size() { \
        return sizeof(CustomVector) / sizeof(double); } \
    \
    template<typename OtherDerived> \
    CustomVector(const Eigen::MatrixBase<OtherDerived>& other) { \
        for(int i = 0; i < CustomVector::size(); ++i){ \
            ((double*)this)[i] = other(i); \
        } \
    } \
    \
    template<typename OtherDerived> \
    CustomVector& operator=(const Eigen::MatrixBase<OtherDerived>& other) { \
        *reinterpret_cast<Eigen::Matrix<double, CustomVector::size(), 1>*>(this) = other; \
        return *this; \
    } \
    \
    const auto &vec() const { \
        return *reinterpret_cast<const Eigen::Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    auto &vec() { \
        return *reinterpret_cast<Eigen::Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    inline static auto Zero() { \
        return Eigen::Matrix<double, CustomVector::size(), 1>::Zero(); \
    }

#else

#define VECTORIZE(CustomVector) \
    static constexpr int size() { \
        return sizeof(CustomVector) / sizeof(double); } \
    \
    template<typename _Scalar, int _Rows, int _Cols> \
    CustomVector(const Matrix<_Scalar, _Rows, _Cols>& other) { \
        for(int i = 0; i < CustomVector::size(); ++i){ \
            ((double*)this)[i] = other.get(i); \
        } \
    } \
    \
    const auto &vec() const { \
        return *reinterpret_cast<const Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    auto &vec() { \
        return *reinterpret_cast<Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    inline static auto Zero() { \
        return Matrix<double, CustomVector::size(), 1>::Zero(); \
    }

#endif