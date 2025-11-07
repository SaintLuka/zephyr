/// @brief Численное дифференцирование.

#pragma once

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <functional>

namespace zephyr::math {

/// @brief Центральная разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
/// @tparam F Целевая функция F(double) -> double
/// @param x Аргумент
/// @param h Величина шага разностного шаблона
template <int ord, int acc=6, typename F>
double derivative(F&& f, double x, double h);

/// @brief Правая разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
template <int ord, int acc=4, typename F>
double derivative_r(F&& f, double x, double h);

/// @brief Левая разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
template <int ord, int acc=4, typename F>
double derivative_l(F&& f, double x, double h);




// Центральные разностные производные

template <typename F>
double derivative_c_o1a2(F&& f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

template <typename F>
double derivative_c_o1a4(F&& f, double x, double h) {
    return (8 * (f(x + h) - f(x - h)) - (f(x + 2 * h) - f(x - 2 * h))) / (12 * h);
}

template <typename F>
double derivative_c_o1a6(F&& f, double x, double h) {
    return (45 * (f(x + h) - f(x - h)) - 9 * (f(x + 2 * h) - f(x - 2 * h)) +
            (f(x + 3 * h) - f(x - 3 * h))) / (60 * h);
}

template <typename F>
double derivative_c_o2a2(F&& f, double x, double h) {
    return (-2 * f(x) + (f(x + h) + f(x - h))) / (h * h);
}

template <typename F>
double derivative_c_o2a4(F&& f, double x, double h) {
    return (-30 * f(x) + 16 * (f(x + h) + f(x - h)) - (f(x + 2 * h) + f(x - 2 * h))) / (12 * h * h);
}

template <typename F>
double derivative_c_o2a6(F&& f, double x, double h) {
    return (-245 * f(x) + 135 * (f(x + h) + f(x - h)) - 13.5 * (f(x + 2 * h) + f(x - 2 * h)) +
            (f(x + 3 * h) + f(x - 3 * h))) / (90 * h * h);
}

// Правые разностные производные

template <typename F>
double derivative_r_o1a1(F&& f, double x, double h) {
    const double c0 = -1.0;
    const double c1 = +1.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h)) * v;
}

template <typename F>
double derivative_r_o1a2(F&& f, double x, double h) {
    const double c0 = -3.0 / 2.0;
    const double c1 = +2.0;
    const double c2 = -1.0 / 2.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h)) * v;
}

template <typename F>
double derivative_r_o1a3(F&& f, double x, double h) {
    const double c0 = -11.0 / 6.0;
    const double c1 = +3.0;
    const double c2 = -3.0 / 2.0;
    const double c3 = +1.0 / 3.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h)) * v;
}

template <typename F>
double derivative_r_o1a4(F&& f, double x, double h) {
    const double c0 = -25.0 / 12.0;
    const double c1 = +4.0;
    const double c2 = -3.0;
    const double c3 = +4.0 / 3.0;
    const double c4 = -1.0 / 4.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) +
            c3 * f(x + 3 * h) + c4 * f(x + 4 * h)) * v;
}

template <typename F>
double derivative_r_o1a5(F&& f, double x, double h) {
    const double c0 = -137.0 / 60.0;
    const double c1 = +5.0;
    const double c2 = -5.0;
    const double c3 = +10.0 / 3.0;
    const double c4 = -5.0 / 4.0;
    const double c5 = +1.0 / 5.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h) +
            c4 * f(x + 4 * h) + c5 * f(x + 5 * h)) * v;
}

template <typename F>
double derivative_r_o1a6(F&& f, double x, double h) {
    const double c0 = -49.0 / 20.0;
    const double c1 = +6.0;
    const double c2 = -15.0 / 2.0;
    const double c3 = +20.0 / 3.0;
    const double c4 = -15.0 / 4.0;
    const double c5 = +6.0 / 5.0;
    const double c6 = -1.0 / 6.0;
    const double v = 1.0 / h;
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h) +
            c4 * f(x + 4 * h) + c5 * f(x + 5 * h) + c6 * f(x + 6 * h)) * v;
}

template <typename F>
double derivative_r_o2a1(F&& f, double x, double h) {
    const double c0 = +1.0;
    const double c1 = -2.0;
    const double c2 = +1.0;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h)) * v;
}

template <typename F>
double derivative_r_o2a2(F&& f, double x, double h) {
    const double c0 = +2.0;
    const double c1 = -5.0;
    const double c2 = +4.0;
    const double c3 = -1.0;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h)) * v;
}

template <typename F>
double derivative_r_o2a3(F&& f, double x, double h) {
    const double c0 = +35.0 / 12.0;
    const double c1 = -26.0 / 3.0;
    const double c2 = +19.0 / 2.0;
    const double c3 = -14.0 / 3.0;
    const double c4 = +11.0 / 12;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) +
            c3 * f(x + 3 * h) + c4 * f(x + 4 * h)) * v;
}

template <typename F>
double derivative_r_o2a4(F&& f, double x, double h) {
    const double c0 = +15.0 / 4.0;
    const double c1 = -77.0 / 6.0;
    const double c2 = +107.0 / 6.0;
    const double c3 = -13.0;
    const double c4 = +61.0 / 12.0;
    const double c5 = -5.0 / 6.0;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h) +
            c4 * f(x + 4 * h) + c5 * f(x + 5 * h)) * v;
}

template <typename F>
double derivative_r_o2a5(F&& f, double x, double h) {
    const double c0 = +203.0 / 45.0;
    const double c1 = -87.0 / 5.0;
    const double c2 = +117.0 / 4.0;
    const double c3 = -254.0 / 9.0;
    const double c4 = +33.0 / 2.0;
    const double c5 = -27.0 / 5.0;
    const double c6 = +137.0 / 180.0;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h) +
            c4 * f(x + 4 * h) + c5 * f(x + 5 * h) + c6 * (f(x + 6 * h))) * v;
}

template <typename F>
double derivative_r_o2a6(F&& f, double x, double h) {
    const double c0 = +469.0 / 90.0;
    const double c1 = -223 / 10.0;
    const double c2 = +879.0 / 20.0;
    const double c3 = -949.0 / 18.0;
    const double c4 = +41.0;
    const double c5 = -201.0 / 10.0;
    const double c6 = +1019.0 / 180.0;
    const double c7 = -7.0 / 10.0;
    const double v = 1.0 / (h * h);
    return (c0 * f(x) + c1 * f(x + h) + c2 * f(x + 2 * h) + c3 * f(x + 3 * h) +
            c4 * f(x + 4 * h) + c5 * f(x + 5 * h) + c6 * f(x + 6 * h) + c7 * f(x + 7 * h)) * v;
}




/// @brief Центральная разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
template <int ord, int acc, typename F>
double derivative(F&& f, double x, double h) {
    if (ord < 2) {
        switch (acc) {
            case 1:
            case 2:
                return derivative_c_o1a2(f, x, h);
            case 3:
            case 4:
                return derivative_c_o1a4(f, x, h);
            default:
                return derivative_c_o1a6(f, x, h);
        }
    }
    else if (ord < 3) {
        switch (acc) {
            case 1:
            case 2:
                return derivative_c_o2a2(f, x, h);
            case 3:
            case 4:
                return derivative_c_o2a4(f, x, h);
            default:
                return derivative_c_o2a6(f, x, h);
        }
    }
    else {
        throw std::runtime_error("derivative third order not implemented");
    }
}

/// @brief Левая разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
template <int ord, int acc, typename F>
double derivative_l(F&& f, double x, double h) {
    if (ord < 2) {
        switch (acc) {
            case 1:
                return derivative_r_o1a1(f, x, -h);
            case 2:
                return derivative_r_o1a2(f, x, -h);
            case 3:
                return derivative_r_o1a3(f, x, -h);
            case 4:
                return derivative_r_o1a4(f, x, -h);
            case 5:
                return derivative_r_o1a5(f, x, -h);
            default:
                return derivative_r_o1a6(f, x, -h);
        }
    }
    else if (ord < 3) {
        switch (acc) {
            case 1:
                return derivative_r_o2a1(f, x, -h);
            case 2:
                return derivative_r_o2a2(f, x, -h);
            case 3:
                return derivative_r_o2a3(f, x, -h);
            case 4:
                return derivative_r_o2a4(f, x, -h);
            case 5:
                return derivative_r_o2a5(f, x, -h);
            default:
                return derivative_r_o2a6(f, x, -h);
        }
    }
    else {
        throw std::runtime_error("derivative third order not implemented");
    }
}

/// @brief Правая разностная производная
/// @tparam ord Порядок производной
/// @tparam acc Порядок аппроксимации не меньше
template <int ord, int acc, typename F>
double derivative_r(F&& f, double x, double h) {
    if (ord < 2) {
        switch (acc) {
            case 1:
                return derivative_r_o1a1(f, x, h);
            case 2:
                return derivative_r_o1a2(f, x, h);
            case 3:
                return derivative_r_o1a3(f, x, h);
            case 4:
                return derivative_r_o1a4(f, x, h);
            case 5:
                return derivative_r_o1a5(f, x, h);
            default:
                return derivative_r_o1a6(f, x, h);
        }
    }
    else if (ord < 3) {
        switch (acc) {
            case 1:
                return derivative_r_o2a1(f, x, h);
            case 2:
                return derivative_r_o2a2(f, x, h);
            case 3:
                return derivative_r_o2a3(f, x, h);
            case 4:
                return derivative_r_o2a4(f, x, h);
            case 5:
                return derivative_r_o2a5(f, x, h);
            default:
                return derivative_r_o2a6(f, x, h);
        }
    }
    else {
        throw std::runtime_error("derivative third order not implemented");
    }
}


enum class DiffStyle : int {
    Central  = 0,
    Forward  = 1,
    Backward = 2,
};

// DERIVEST calculates numeric derivative of a function.
// Arguments (input)
//   fun              - function to differentiate
//   x0               - point at which to differentiate fun
//   derivative_order - specifies the derivative order estimated.
//                      Must be a positive integer from the set {1, 2, 3, 4}.
//   method_order     - specifies the order of the basic method used for
//                      the estimation.
//                      For central methods, methods_order must be 2 or 4;
//                      otherwise can be 1, 2, 3 or 4
//   style            - specifies the style of the basic method used for the
//                      estimation. 'central', 'forward', or 'backwards'
//                      difference methods are used.
//   romberg_terms    - Allows the user to specify the generalized Romberg
//                      extrapolation method used, or turn it off completely.
//                      Must be a positive integer from the set {0, 1, 2, 3}.
//
// Arguments (output)
//   &der             - store the result (derivative)
//   &err             - store error estimate
//   &final_delta     - store the final overall stepsize chosen by DERIVEST
//
// Return value: true if succeeded, false otherwise.
//
// Remark. Recommended params are as follows:
//   diff_style    = DiffStyle::Central
//   method_order  = 4
//   romberg_terms = 2
bool derivest_c(const std::function<double(double)>& fun, double x0,
                int derivative_order, int method_order,
                DiffStyle diff_style, int romberg_terms,
                double &der, double& err, double& final_delta);

/// @brief Вызов функции derivest_c на манер C++
/// По хорошему derivest_c тоже можно переписать на шаблон. Но зачем?
/// Так симпатичнее же? df = derivest<1>(func, x);
template<int der_ord, int acc_ord = 4, DiffStyle style = DiffStyle::Central>
double derivest(const std::function<double(double)>& func, double x) {
    double res{NAN}, err{NAN}, delta{NAN};
    bool conv = derivest_c(func, x, der_ord, acc_ord, style, 2, res, err, delta);
    return conv ? res : NAN;
}

} // namespace zephyr::math