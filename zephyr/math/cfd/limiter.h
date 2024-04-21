#pragma once

#include <string>

namespace zephyr { namespace math {

/// @brief Класс-функтор, предоставляет функционал ограничителей
/// наклонов/потоков для скалярных и векторных функций
class Limiter {
public:

    /// @brief Ограничитель по умолчанию minmod
    Limiter();

    /// @brief Создать ограничитель по названию
    Limiter(const char* name);

    /// @brief Создать ограничитель по названию
    Limiter(const std::string& name);

    /// @brief Название ограничителя
    const std::string& name() const;

    /// @brief Ограничитель включен?
    bool is_on() const;

    /// @brief Унарная функция ограничителя
    double operator()(const double& val) const;

    /// @brief Бинарная функция ограничителя = limiter(val1 / val2)
    double operator()(const double& val1, const double& val2) const;

    /// @brief Размер произвольного типа в double
    template <class T>
    static int SizeOf() {
        return int(sizeof(T) / sizeof(double));
    }

    /// @brief Унарная функция ограничителя для векторных функций
    template <class T>
    T operator()(const T& val) const {
        T result;
        auto out = (double *) &result;
        auto in1 = (const double *) &val;

        for (int i = 0; i < SizeOf<T>(); ++i) {
            out[i] = func1(in1[i]);
        }
        return result;
    }

    /// @brief Бинарная функция ограничителя для векторных функций
    template <class T>
    T operator()(const T& val1, const T& val2) const {
        T result;
        auto out = (double *) &result;
        auto in1 = (const double *) &val1;
        auto in2 = (const double *) &val2;

        for (int i = 0; i < SizeOf<T>(); ++i) {
            out[i] = func2(in1[i], in2[i]);
        }
        return result;
    }

    /// @brief Тестирует все ограничители
    static void check();

private:
    bool m_active;                    ///< Ограничитель выключен
    std::string m_name;               ///< Название ограничителя
    double (*func1)(double);          ///< Унарная функция ограничителя
    double (*func2)(double, double);  ///< Бинарная функция ограничителя
};

/// @brief Коллекция ограничителей
namespace limiters {

double minmod(double r);
double minmod(double a, double b);

double MC(double r);
double MC(double a, double b);

double vanLeer(double r);
double vanLeer(double a, double b);

double vanAlbada1(double r);
double vanAlbada1(double a, double b);

double vanAlbada2(double r);
double vanAlbada2(double a, double b);

double CHARM(double r);
double CHARM(double a, double b);

double HCUS(double r);
double HCUS(double a, double b);

double HQUICK(double r);
double HQUICK(double a, double b);

double Koren(double r);
double Koren(double a, double b);

double smart(double r);
double smart(double a, double b);

double superbee(double r);
double superbee(double a, double b);

double atvl(double r);
double atvl(double a, double b);

}

} // math
} // zephyr