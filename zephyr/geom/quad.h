#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Представление четырехугольника. Содержит несколько полезных
/// функций, а также отвечает за линейное отображение квадрата
/// (x, y) in [-1, 1]^2 на заданый четырехугольник.
struct Quad {
protected:
    /// @brief Вершины в Z-порядке (таблица 2 x 2)
    std::array<Vector3d, 4> verts;

public:
    /// @brief Конструктор по угловым точкам (Z-порядок)
    Quad(const Vector3d& v00,
         const Vector3d& v01,
         const Vector3d& v10,
         const Vector3d& v11);

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..3]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..3]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индексы {-1, +1}^2 в [0..3]
    template <int i, int j>
    inline static constexpr int iss() {
        static_assert(i * i == 1 && j * j == 1,
                      "Available indices: {-1, +1}");
        return j + 1 + (i + 1) / 2;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, +1}
    /// @tparam j in {-1, +1}
    template <int i, int j>
    inline Vector3d &vs() {
        return verts[iss<i, j>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, +1}
    /// @tparam j in {-1, +1}
    template <int i, int j>
    inline const Vector3d &vs() const {
        return verts[iss<i, j>()];
    }

    /// @brief Отображение из квадрата [-1, 1]^2
    /// @details В угловых точках совпадает с операторами выше
    Vector3d operator()(double x, double y) const;

    /// @brief Отображение из квадрата [-1, 1]^2
    /// @details Дублирует оператор доступа operator()(x, y)
    Vector3d get(double x, double y) const;

    /// @brief Касательный вектор в точке отображения
    /// (производная отображения по параметру x)
    Vector3d tangent_x(double x, double y) const;

    /// @brief Касательный вектор в точке отображения
    /// (производная отображения по параметру y)
    Vector3d tangent_y(double x, double y) const;

    /// @brief Нормаль к поверхности в точке, направлена от точки c
    Vector3d normal(double x, double y, const Vector3d& c) const;

    /// @brief Якобиан отображения в точке
    double Jacobian(double x, double y) const;

    /// @brief Центр простой двумерной ячейки
    Vector3d center() const;

    /// @brief Внешняя нормаль к простой двумерной грани
    /// @param c Центр ячейки
    Vector3d normal(const Vector3d& c) const;

    /// @brief "Площадь" простой двумерной грани. Равна | int vec(n) dS |,
    /// то есть модулю поверхностного интеграла 2-го рода.
    double area() const;

    /// @brief Объем обычной двумерной ячейки в осесимметричной постновке
    /// Вершины должны располагаться в плоскости (x, y), где y - радиус
    double volume_as() const;

    /// @brief Барицентр простой двумерной ячейки.
    /// Вершины должны располагаться в плоскости (x, y).
    /// @param area Площадь ячейки (если известна)
    Vector3d centroid(double area = 0.0) const;

    /// @brief Посчитать объемную долю, которая отсекается от ячейки некоторым
    /// телом, точки которого определяются характеристической функцией inside
    /// @param inside Характеристическая функция: true, если точка находится
    /// внутри тела, иначе false
    /// @param n_points Число тестовых точек, погрешность ~ 1/N.
    double volume_fraction(const std::function<bool(const Vector3d&)>& inside,
            int n_points = 10000) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод Гаусса 4-го порядка (middle accuracy order)
    double integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод 8-го порядка по 12 узлам (high accuracy order)
    double integrate_high(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод 13-го порядка по 37 узлам (extra-high  accuracy order)
    double integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const;
};

/// @brief Квадратичное отображение на четырехугольник
struct SqQuad {
protected:
    /// @brief Вершины (таблица 3 x 3)
    std::array<Vector3d, 9> verts;

public:
    /// @brief Конструктор по умолчанию
    SqQuad() = default;

    /// @brief Конструктор по угловым точкам
    SqQuad(const Vector3d& v00,
           const Vector3d& v01,
           const Vector3d& v10,
           const Vector3d& v11);

    /// @brief Конструктор по полному набору вершин
    SqQuad(const Vector3d& v00, const Vector3d& v01, const Vector3d& v02,
           const Vector3d& v10, const Vector3d& v11, const Vector3d& v12,
           const Vector3d& v20, const Vector3d& v21, const Vector3d& v22);

    /// @brief Конструктор из простого четырехугольника
    SqQuad(const Quad& quad);

    /// @brief Удалить промежуточные вершины и вернуть
    /// простой четырехугольник
    Quad reduce() const;

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..8]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..8]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индексы {-1, 0, +1}^2 в [0..8]
    template<int i, int j>
    inline static constexpr int iss() {
        static_assert(i * i <= 1 && j * j <= 1,
                      "Available indices: {-1, 0, 1}");
        return 3 * j + i + 4;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @param i in {-1, 0, +1}
    /// @param j in {-1, 0, +1}
    template <int i, int j>
    inline Vector3d &vs() {
        return verts[iss<i, j>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @param i in {-1, 0, +1}
    /// @param j in {-1, 0, +1}
    template <int i, int j>
    inline const Vector3d &vs() const {
        return verts[iss<i, j>()];
    }

    /// @brief Отображение из квадрата [-1, 1]^2
    /// @details В угловых точках совпадает с операторами выше
    Vector3d operator()(double x, double y) const;

    /// @brief Отображение из квадрата [-1, 1]^2
    /// @details Дублирует оператор доступа operator()(x, y)
    Vector3d get(double x, double y) const;

    /// @brief Касательный вектор в точке отображения
    /// (производная отображения по параметру x)
    Vector3d tangent_x(double x, double y) const;

    /// @brief Касательный вектор в точке отображения
    /// (производная отображения по параметру y)
    Vector3d tangent_y(double x, double y) const;

    /// @brief Нормаль к поверхности в точке, направлена от точки c
    Vector3d normal(double x, double y, const Vector3d& c) const;

    /// @brief Якобиан отображения в точке
    double Jacobian(double x, double y) const;

    /// @brief Центр криволинейной двумерной ячейки
    const Vector3d& center() const;

    /// @brief Внешняя нормаль к простой двумерной грани
    /// @param c Центр ячейки
    Vector3d normal(const Vector3d& c) const;

    /// @brief Площадь криволинейной двумерной ячейки.
    /// Вершины должны располагаться в плоскости (x, y).
    double area() const;

    /// @brief Объем криволинейной двумерной ячейки в осесимметричной постновке
    /// Вершины должны располагаться в плоскости (x, y), где y - радиус
    double volume_as() const;

    /// @brief Барицентр криволинейной двумерной ячейки.
    /// Вершины должны располагаться в плоскости (x, y).
    /// @param area Площадь криволинейной ячейки
    Vector3d centroid(double area = 0.0) const;

    /// @brief Дочерние ячейки после операции split на адаптивных сетках
    std::array<SqQuad, 4> children() const;

    /// @brief Посчитать объемную долю, которая отсекается от ячейки некоторым
    /// телом, точки которого определяются характеристической функцией inside
    /// @param inside Характеристическая функция: true, если точка находится
    /// внутри тела, иначе false
    /// @param n_points Число тестовых точек, погрешность ~ 1/N.
    double volume_fraction(const std::function<bool(const Vector3d&)>& inside,
                           int n_points = 10000) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод Гаусса 4-го порядка (middle accuracy order)
    double integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод 8-го порядка по 12 узлам (high accuracy order)
    double integrate_high(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по квадратному элементу
    /// @param n Число подъячеек по осям
    /// @details Метод 13-го порядка по 37 узлам (extra-high  accuracy order)
    double integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const;
};

} // namespace zephyr::geom