#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom::curves {

/// @brief Кривая Безье, параметризация на отрезке [0, 1].
class Bezier {
public:

    /// @brief Конструктор кривой Безье по узлам
    Bezier(const std::vector<Vector3d>& vs);

    /// @brief Конструктор кривой Безье по координатам узлов в нескольких массивах
    Bezier(std::vector<double> xs, std::vector<double> ys);

    /// @brief Конструктор кривой Безье по координатам узлов в нескольких массивах
    Bezier(std::vector<double> xs, std::vector<double> ys, std::vector<double> zs);


    /// @brief Получить точку на кривой
    Vector3d get(double t) const;

    /// @brief Касательный вектор к кривой
    Vector3d tangent(double t) const;

    /// @brief Нормаль к кривой
    /// @param c Точка, от которой направлена нормаль, если не задана,
    /// то используется правая нормаль (в двумерном случае).
    Vector3d normal(double t, Vector3d p = {NAN, NAN, NAN}) const;


    /// @brief x-координаты узлов кривой Безье
    const std::vector<double>& xs() const;

    /// @brief y-координаты узлов кривой Безье
    const std::vector<double>& ys() const;

    /// @brief z-координаты узлов кривой Безье
    const std::vector<double>& zs() const;

    /// @brief Узлы кривой Безье
    std::vector<Vector3d> vs() const;


    /// @brief Массив N значений x-координат, параметр t на отрезке [0, 1]
    std::vector<double> xs(int N) const;

    /// @brief Массив N значений y-координат, параметр t на отрезке [0, 1]
    std::vector<double> ys(int N) const;

    /// @brief Массив N значений z-координат, параметр t на отрезке [0, 1]
    std::vector<double> zs(int N) const;

    /// @brief Массив N точек на кривой, параметр t на отрезке [0, 1]
    std::vector<Vector3d> vs(int N) const;


    /// @brief Массив N значений x-координат, параметр t на отрезке [t1, t2]
    std::vector<double> xs(int N, double t1, double t2) const;

    /// @brief Массив N значений y-координат, параметр t на отрезке [t1, t2]
    std::vector<double> ys(int N, double t1, double t2) const;

    /// @brief Массив N значений z-координат, параметр t на отрезке [t1, t2]
    std::vector<double> zs(int N, double t1, double t2) const;

    /// @brief Массив N точек на кривой, параметр t на отрезке [t1, t2]
    std::vector<Vector3d> vs(int N, double t1, double t2) const;

private:
    std::vector<double> m_xs;  ///< x-координаты узлов кривой Безье
    std::vector<double> m_ys;  ///< y-координаты узлов кривой Безье
    std::vector<double> m_zs;  ///< z-координаты узлов кривой Безье
};

} // namespace zephyr::geom::curves