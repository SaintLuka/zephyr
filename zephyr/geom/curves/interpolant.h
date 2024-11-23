#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom::curves {

/// @brief Типы границ сплайна и экстраполяции за границу.
/// Для разных сплайнов могут быть недоступны все варианты границ.
enum class SplineBound {
    Crop,      ///< Завершить на последнем значении
    Free,      ///< Линейное продолжение
    Periodic,  ///< Периодическое замыкание
    Warning,   ///< Завершить на последнем, выводить предупреждение
    ///  при выходе за границы диапаза
};

/// @brief Параметризация параметрических кривых.
/// При любой параметризации первый узел = 0, последний = 1.
enum class Parametrization {
    Uniform,    ///< Равномерно между узлами
    Chord,      ///< Пропорционально длине хорд между узлами
    Chebyshev,  ///< Корни полинома Чебышёва
};


/// @brief Поместить параметр x в период между [x_min, x_max]
void fit_into_period(double x_min, double x_max, double& x);

/// @brief Найти номер сегмента n: xs[n] < x < xs[n + 1]
int find_segment(const std::vector<double>& xs, double x);

/// @brief Проверить массив на строгую монотоность
/// @return +1 -- строго возрастает,
///         -1 -- строго убывает,
///          0 -- иначе.
int monotonic(const std::vector<double>& arr);

/// @brief Параметризация по длине отрезков (хорд)
std::vector<double> chord_parametrization(
        const std::vector<double>& xs,
        const std::vector<double>& ys,
        const std::vector<double>& zs);

/// @brief Параметризация по корням полиномов Чебышёва.
/// Корни масштабируются на отрезок [0, 1], крайние корни на границах.
std::vector<double> chebyshev_parametrization(int n);


/// @brief Абстрактный класс интерполяции.
/// Изначально хотел назвать Spline, но от него наследуется класс Lagrange,
/// с полиномиальной интерполяцией Лагранжа, который не является сплайном.
class Interpolant {
public:
    /// @brief Конструктор по умолчанию
    Interpolant() = default;


    /// @brief Минимальное значение аргумента
    double x_min() const;

    /// @brief Максимальное значение аргумента
    double x_max() const;


    /// @brief Основная функция. Значение функции от аргумента
    virtual double get(double x) const = 0;

    /// @brief Значение функции от аргумента (синоним get)
    inline double y(double x) const { return get(x); }

    /// @brief Значение функции от аргумента (синоним get)
    inline double operator()(double x) const { return get(x); };


    /// @brief Значения аргумента в узлах
    const std::vector<double>& xs() const;

    /// @brief Значения функции в узлах
    const std::vector<double>& ys() const;


    /// @brief Массив N значений аргумента на отрезке [x_min, x_max]
    std::vector<double> xs(int N) const;

    /// @brief Массив N значений фукнции для аргументов на отрезке [x_min, x_max]
    std::vector<double> ys(int N) const;


    /// @brief Массив N значений аргумента на отрезке [x1, x2]
    std::vector<double> xs(int N, double x1, double x2) const;

    /// @brief Массив N значений фукнции для аргументов на отрезке [x1, x2]
    std::vector<double> ys(int N, double x1, double x2) const;

protected:
    /// @brief Левая и правая граница сплайна
    SplineBound m_left, m_right;
    std::vector<double> m_xs;  ///< x-координаты узлов
    std::vector<double> m_ys;  ///< y-координаты узлов
};


/// @brief Абстрактный класс параметрической интерполяции.
/// Набор параметрических кривых: x(t), y(t), z(t), значения параметра t ∈ [0, 1].
class PInterpolant {
public:
    /// @brief Конструктор по умолчанию
    PInterpolant() = default;


    /// @brief Минимальное значение параметра
    double t_min() const;

    /// @brief Максимальное значение параметра
    double t_max() const;


    //// @brief Получить точку на кривой
    Vector3d get(double t) const;

    /// @brief Получить точку на ломаной (синоним get)
    Vector3d v(double t) const { return get(t); };

    /// @brief Получить точку на ломаной (синоним get)
    Vector3d operator()(double t) const { return get(t); };


    /// @brief Получить точку на кривой
    virtual double x(double t) const = 0;

    /// @brief Получить точку на кривой
    virtual double y(double t) const = 0;

    /// @brief Получить точку на кривой
    virtual double z(double t) const = 0;


    /// @brief Узлы ломаной
    std::vector<Vector3d> vs() const;

    /// @brief Значения параметра в узлах
    const std::vector<double>& ts() const;

    /// @brief Значения x-координат в узлах
    const std::vector<double>& xs() const;

    /// @brief Значения y-координат в узлах
    const std::vector<double>& ys() const;

    /// @brief Значения z-координат в узлах
    const std::vector<double>& zs() const;


    /// @brief Массив N точек на кривой, параметр t ∈ [t_min, t_max]
    std::vector<Vector3d> vs(int N) const;

    /// @brief Массив N значений параметров, параметр t ∈ [t_min, t_max]
    std::vector<double> ts(int N) const;

    /// @brief Массив N значений x-координат, параметр t ∈ [t_min, t_max]
    std::vector<double> xs(int N) const;

    /// @brief Массив N значений y-координат, параметр t ∈ [t_min, t_max]
    std::vector<double> ys(int N) const;

    /// @brief Массив N значений z-координат, параметр t ∈ [t_min, t_max]
    std::vector<double> zs(int N) const;


    /// @brief Массив N точек на кривой, параметр t ∈ [t1, t2]
    std::vector<Vector3d> vs(int N, double t1, double t2) const;

    /// @brief Массив N значений параметров, параметр t ∈ [t1, t2]
    std::vector<double> ts(int N, double t1, double t2) const;

    /// @brief Массив N значений x-координат, параметр t ∈ [t1, t2]
    std::vector<double> xs(int N, double t1, double t2) const;

    /// @brief Массив N значений y-координат, параметр t ∈ [t1, t2]
    std::vector<double> ys(int N, double t1, double t2) const;

    /// @brief Массив N значений z-координат, параметр t ∈ [t1, t2]
    std::vector<double> zs(int N, double t1, double t2) const;

protected:
    /// @brief Левая и правая границы сплайна
    SplineBound m_left, m_right;

    std::vector<double> m_ts;  ///< Значения параметра
    std::vector<double> m_xs;  ///< Координаты x узлов
    std::vector<double> m_ys;  ///< Координаты y узлов
    std::vector<double> m_zs;  ///< Координаты z узлов
};

} // namespace zephyr::geom::curves