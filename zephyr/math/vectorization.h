#pragma once

#define VECTORIZE(CustomVector) \
    static constexpr int size() { \
        return sizeof(CustomVector) / sizeof(double); } \
    \
    template<typename OtherDerived> \
    CustomVector(const Eigen::MatrixBase<OtherDerived>& other) { \
        *reinterpret_cast<Eigen::Matrix<double, CustomVector::size(), 1>*>(this) = other; \
    } \
    \
    template<typename OtherDerived> \
    CustomVector& operator=(const Eigen::MatrixBase<OtherDerived>& other) { \
        *reinterpret_cast<Eigen::Matrix<double, CustomVector::size(), 1>*>(this) = other; \
        return *this; \
    } \
    \
    const auto vec() const { \
        return *reinterpret_cast<const Eigen::Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    auto vec() { \
        return *reinterpret_cast<Eigen::Matrix<double, CustomVector::size(), 1> *>(this); \
    } \
    \
    inline static auto Zero() { \
        return Eigen::Matrix<double, CustomVector::size(), 1>::Zero(); \
    }
