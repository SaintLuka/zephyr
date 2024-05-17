#pragma once

#include <array>
#include <cmath>

namespace zephyr::math {

/// @brief Интерполяционный шаблон на ячейках с индексами от b до e 
/// включительно. К примеру, WENO<-1, 1> -- центральный шаблон на 
/// трех ячейках с индексами {-1, 0, 1}.
/// Справедливо для интерполяции на равномерной сетке.
template<int b, int e>
class WENO {
public:
    /// @brief Длина шаблона (число ячеек)
    static const int len = e + 1 - b;

    /// @brief Конструктор
    /// @param u Массив усредненных значений в ячейках
    WENO(const std::array<double, len> &u)
            : m_data(u) {};

    /// @brief Интерполяция на правой грани ячейки,
    /// то есть значение функции в точке u(+1/2).
    double p() const;

    /// @brief Интерполяция на левой грани ячейки, 
    /// то есть значение функции в точке u(-1/2).
    double m() const;

    /// @brief Индикатор гладкости на шаблоне,
    /// для гладких функций O(h^2).
    double beta() const;

private:
    std::array<double, len> m_data;
};


template<>
double WENO<-2, 0>::p() const {
    const double w1 = +1.0 / 3.0;
    const double w2 = -7.0 / 6.0;
    const double w3 = 11.0 / 6.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<-1, 1>::p() const {
    const double w1 = -1.0 / 6.0;
    const double w2 = +5.0 / 6.0;
    const double w3 = +1.0 / 3.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<0, 2>::p() const {
    const double w1 = +1.0 / 3.0;
    const double w2 = +5.0 / 6.0;
    const double w3 = -1.0 / 6.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<-2, 0>::m() const {
    const double w1 = -1.0 / 6.0;
    const double w2 = +5.0 / 6.0;
    const double w3 = +1.0 / 3.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<-1, 1>::m() const {
    const double w1 = +1.0 / 3.0;
    const double w2 = +5.0 / 6.0;
    const double w3 = -1.0 / 6.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<0, 2>::m() const {
    const double w1 = 11.0 / 6.0;
    const double w2 = -7.0 / 6.0;
    const double w3 = +1.0 / 3.0;
    return w1 * m_data[0] + w2 * m_data[1] + w3 * m_data[2];
}

template<>
double WENO<-2, 0>::beta() const {
    double d1 = m_data[0] - 4.0 * m_data[1] + 3.0 * m_data[2];
    double d2 = m_data[0] - 2.0 * m_data[1] + m_data[2];
    return 13.0 / 12.0 * d2 * d2 + 1.0 / 4.0 * d1 * d1;
}

template<>
double WENO<-1, 1>::beta() const {
    double d1 = m_data[0] - m_data[2];
    double d2 = m_data[0] - 2.0 * m_data[1] + m_data[2];
    return 13.0 / 12.0 * d2 * d2 + 1.0 / 4.0 * d1 * d1;
}

template<>
double WENO<0, 2>::beta() const {
    double d1 = 3.0 * m_data[0] - 4.0 * m_data[1] + m_data[2];
    double d2 = m_data[0] - 2.0 * m_data[1] + m_data[2];
    return 13.0 / 12.0 * d2 * d2 + 1.0 / 4.0 * d1 * d1;
}

inline double sqr(double x) {
    return x * x;
}

/// @brief WENO-реконструкция на пятиточечном шаблоне
class WENO5 {
public:

    WENO5(double u1, double u2, double u3, double u4, double u5)
            : we2({u1, u2, u3}), we1({u2, u3, u4}), we0({u3, u4, u5}) {

    }

    /// @brief Интерполяция на правой грани ячейки,
    /// то есть значение функции в точке u(+1/2).
    double p() const {
        const double d0_c = 3.0 / 10.0;
        const double d1_c = 6.0 / 10.0;
        const double d2_c = 1.0 / 10.0;

        double b0 = we0.beta();
        double b1 = we1.beta();
        double b2 = we2.beta();

        double d0 = d0_c / sqr(1.0e-8 + b0);
        double d1 = d1_c / sqr(1.0e-8 + b1);
        double d2 = d2_c / sqr(1.0e-8 + b2);

        double w0 = d0 / (d0 + d1 + d2);
        double w1 = d1 / (d0 + d1 + d2);
        double w2 = d2 / (d0 + d1 + d2);

        return w0 * we0.p() + w1 * we1.p() + w2 * we2.p();
    }

    /// @brief Интерполяция на левой грани ячейки,
    /// то есть значение функции в точке u(-1/2).
    double m() const {
        const double d0_c = 1.0 / 10.0;
        const double d2_c = 3.0 / 10.0;
        const double d1_c = 6.0 / 10.0;

        double b0 = we0.beta();
        double b1 = we1.beta();
        double b2 = we2.beta();

        double d0 = d0_c / sqr(1.0e-8 + b0);
        double d1 = d1_c / sqr(1.0e-8 + b1);
        double d2 = d2_c / sqr(1.0e-8 + b2);

        double w0 = d0 / (d0 + d1 + d2);
        double w1 = d1 / (d0 + d1 + d2);
        double w2 = d2 / (d0 + d1 + d2);

        return w0 * we0.m() + w1 * we1.m() + w2 * we2.m();
    }

public:
    WENO<-2, 0> we2;
    WENO<-1, 1> we1;
    WENO<0, 2> we0;
};

}