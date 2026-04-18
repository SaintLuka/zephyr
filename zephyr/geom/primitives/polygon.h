#pragma once

#include <span>
#include <vector>
#include <functional>

#include <zephyr/geom/box.h>
#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @addtogroup geom-primitives
/// @{

/// @brief Плоский многоугольник в плоскости Oxy или полигон в пространстве.
/// Для треугольника лучше использовать класс Triangle.
/// @details По умолчанию вершины полигона не сортируются, если порядок
/// вершин в конструкторе неизвестен, то лучше вызвать функцию Polygon::sort().
/// Полигон можно построить перемещением вершин из std::vector.
class Polygon {
public:
    /// @brief Пустой полигон
    Polygon();

    /// @brief Вершины полигона можно напрямую задавать в фигурных скобках
    /// Vector3d v1, v2, v3, v4;
    /// ...
    /// Polygon poly = {v1, v2, v3, v4};
    Polygon(std::initializer_list<Vector3d> list);

    /// @brief Создание полигона на основе имеющегося массива вершин
    /// заданного размера. Вершины копируются в буфер полигона
    explicit Polygon(std::span<const Vector3d> vertices);

    /// @brief Создание полигона на основе имеющегося массива вершин
    explicit Polygon(std::vector<Vector3d>&& vertices);

    /// @brief Полигон заглушка?
    bool empty() const { return vs.empty(); }

    /// @brief Число вершин в буфере
    int size() const { return static_cast<int>(vs.size()); }

    /// @brief Оператор доступа по индексу
    const Vector3d& operator[](int idx) const { return vs[idx]; }

    /// @brief Координаты x вершин многоугольника (с замыканием)
    std::vector<double> xs() const;

    /// @brief Координаты y вершин многоугольника (с замыканием)
    std::vector<double> ys() const;

    /// @brief Массив вершин многоугольника
    const std::vector<Vector3d>& vertices() const { return vs; }

    /// @brief Ограничивающий прямоугольник
    Box bbox() const;

    /// @brief Центр полигона (среднее вершин)
    Vector3d center() const;

    /// @brief Сортировать вершины против часовой стрелки (полигон в плоскости Oxy)
    Polygon& sort();

    /// @brief Сортировать вершины по правилу правой руки: нормалью наружу из точки обзора
    Polygon& sort(const Vector3d& view);

    /// @brief Точка внутри полигона?
    bool inside(const Vector3d& p) const;

    /// @brief Внешняя нормаль к полигону из точки view
    Vector3d normal(const Vector3d& view) const;

    /// @brief Площадь многоугольника.
    double area() const;

    /// @brief Барицентр многоугольника.
    /// @param area Площадь многоугольника (если известна)
    Vector3d centroid(double area = 0.0) const;

    /// @brief Объем фигуры вращения (вокруг оси x)
    double volume_as() const;

    /// @brief Площадь внутри многоугольника на пересечении с характеристической
    /// функцией inside. Вычисляется приближенно с точностью ~ 1 / N
    /// @param inside Характеристическая функция, возвращает true, если точка p
    /// находится внутри области, иначе -- false.
    /// @param N Число пробных точек
    double clip_area(const std::function<bool(const Vector3d& p)>& inside, int N = 10000) const;

    /// @brief Отсечь от полигона часть с помощью прямой с внешней нормалью n,
    /// проходящей через точку p
    /// @return Отсеченный многоугольник
    Polygon clip(const Vector3d& p, const Vector3d& n) const;

    /// @brief Площадь многоугольника, отсекаемого прямой с внешней нормалью n,
    /// проходящей через точку p
    double clip_area(const Vector3d& p, const Vector3d& n) const;

    /// @brief Две точки прямой (p1, p2). Точка p1 находится в центре
    /// сечения, точка p2 на многоугольнике, таким образом, что при обходе
    /// в направлении p1 -> p2 область остается слева.
    struct section {
        Vector3d p1, p2;

        section(const Vector3d& v1, const Vector3d& v2) : p1(v1), p2(v2) { }

        operator const Vector3d&() const { return p1; }
    };

    /// @brief Находит отсечение от полигона с заданной объемной долей
    /// @param n Внешняя нормаль прямой
    /// @param alpha Объемная доля
    section find_section(const Vector3d& n, double alpha) const;

    /// @brief Площадь пересечения многоугольника с кругом
    /// @param c, R Центр круга, радиус круга
    double disk_clip_area(const Vector3d& c, double R) const;

    /// @brief Нормаль в точке окружности при пересечении многоугольника с кругом
    /// @param c, R Центр круга, радиус круга
    /// @return Средняя нормаль к окружности внутри многоугольника, если
    /// окружность пересекает многоугольник, нулевая нормаль в остальных случаях.
    Vector3d disk_clip_normal(const Vector3d& c, double R) const;

    /// @brief Доля полигона, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция: true, если точка находится
    /// внутри тела, иначе false
    /// @param n_points Число тестовых точек, погрешность ~ 1/N.
    double volume_fraction(const std::function<bool(const Vector3d&)>& inside,
                           int n_points = 10000) const;

    /// @brief Интеграл скалярной функции по полигону элементу
    /// @param n Разбиение по сторонам
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 3-го порядка по 6 узлам (middle accuracy order)
    double integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 6-го порядка по 12 узлам (high accuracy order)
    double integrate_high(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 13-го порядка по 37 узлам (extra-high accuracy order)
    double integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const;

protected:
    /// @brief Пересчитывает центр полигона
    void calc_params();

    bool plane_oxy{true};                 ///< Плоский полигон в плоскости Oxy
    std::vector<Vector3d> vs{};           ///< Массив вершин
    Vector3d m_center{Vector3d::Zero()};  ///< Центр полигона
};

/// @brief Вывод многоугольника в консоль
static std::ostream& operator<<(std::ostream& os, const Polygon& poly);

/// @}

} // namespace zephyr::geom