#pragma once

#include <zephyr/configuration.h>

#define VECTORIZE(CustomVector) \
    static constexpr int size() { \
        return sizeof(CustomVector) / sizeof(double); } \
    \
    template<typename OtherDerived> \
    CustomVector(const Eigen::ArrayBase<OtherDerived>& other) { \
        for(int i = 0; i < CustomVector::size(); ++i){ \
            ((double*)this)[i] = other(i); \
        } \
    } \
    \
    template<typename OtherDerived> \
    CustomVector& operator=(const Eigen::ArrayBase<OtherDerived>& other) { \
        *reinterpret_cast<Eigen::Array<double, CustomVector::size(), 1>*>(this) = other; \
        return *this; \
    } \
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
    const auto &arr() const { \
        return *reinterpret_cast<const Eigen::Array<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    auto &arr() { \
        return *reinterpret_cast<Eigen::Array<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    inline static auto Zero() { \
        return Eigen::Array<double, CustomVector::size(), 1>::Zero(); \
    } \
    \
    inline static auto NaN() { \
        Eigen::Array<double, CustomVector::size(), 1> res(NAN); \
        return res; \
    }
