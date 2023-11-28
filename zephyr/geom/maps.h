#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Представление отрезка. Содержит несколько полезных функций,
/// а также отвечает за отображение параметра x in [-1, 1] на отрезок
/// и сопутствующие отображению функции.
struct Line {
protected:
    /// @brief Концы отрезка
    std::array<Vector3d, 2> verts;

public:
    /// @brief Конструктор по двум крайним точкам
    Line(const Vector3d &v1, const Vector3d &v2);

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..1]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..1]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индекс {-1, +1} в {0, 1}
    /// @details Данная функция может показаться избыточной,
    /// она существует для унификации интерфейса в линейке классов:
    /// Line, Quad, Cube, SqLine, SqQuad, SqCube.
    template <int i>
    inline static constexpr int iss() {
        static_assert(i * i == 1, "Available indices: {-1, +1}");
        return (i + 1) / 2;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 1}
    template <int i>
    inline Vector3d &vs() {
        return verts[iss<i>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 1}
    template <int i>
    inline const Vector3d &vs() const {
        return verts[iss<i>()];
    }

    /// @brief Получить точку на отрезке
    /// @param x in [-1, 1]
    Vector3d operator()(double x) const;

    /// @brief Получить точку на отрезке без создания экземпляра класса
    static Vector3d get(const Vector3d &v1, const Vector3d &v2, double x);

    /// @brief Нормаль к отрезку, располагается в плоскости с отрезком
    /// и точкой 'c', направлена от точки 'c'.
    static Vector3d normal(const Vector3d& v1, const Vector3d& v2,
                           const Vector3d& c);

    /// @brief Якобиан отображения в точке x.
    /// Фактически половина длины отрезка.
    double Jacobian() const;

    /// @brief Центр отрезка
    Vector3d center() const;

    /// @brief Нормаль к отрезку, располагается в плоскости с отрезком
    /// и точкой 'c', направлена от точки 'c'.
    Vector3d normal(const Vector3d &c) const;

    /// @brief Длина отрезка
    double length() const;
};

/// @brief Криволинейный отрезок: строится по трем точкам в виде
/// квадратичного отображения. Класс содержит несколько полезных функций,
/// а также отвечает за отображение параметра x in [-1, 1] на кривую
/// и сопутствующие отображению функции.
struct SqLine {
protected:
    /// @brief Данные
    std::array<Vector3d, 3> verts;

public:
    /// @brief Конструктор по крайним точкам
    SqLine(const Vector3d& v1, const Vector3d& v2);

    /// @brief Конструктор по трем точкам
    SqLine(const Vector3d& v1, const Vector3d& vc, const Vector3d& v2);

    /// @brief Конструктор из простого отрезка
    SqLine(const Line &vs);

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..2]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..2]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индекс {-1, 0, +1} в {0, 1, 2}
    /// @details Данная функция может показаться избыточной,
    /// она существует для унификации интерфейса в линейке классов:
    /// Line, Quad, Cube, SqLine, SqQuad, SqCube.
    template <int i>
    inline static constexpr int iss() {
        static_assert(i * i <= 1, "Available indices: {-1, 0, +1}");
        return i + 1;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 0, 1}
    template <int i>
    inline Vector3d &vs() {
        return verts[iss<i>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 0, 1}
    template <int i>
    inline const Vector3d &vs() const {
        return verts[iss<i>()];
    }

    /// @brief Получить точку на кривой
    /// @param x in [-1, 1]
    Vector3d operator()(double x) const;

    /// @brief Получить точку на кривой
    /// @details Дубликат оператора доступа operator()(x)
    Vector3d get(double x) const;

    /// @brief Получить точку на кривой без создания экземпляра класса
    static Vector3d get(const Vector3d &v1, const Vector3d &vc,
                        const Vector3d &v2, double x);

    /// @brief Касательный вектор в точке x (производная отображения)
    Vector3d tangent(double x) const;

    /// @brief Получить касательный вектор без создания экземпляра класса
    static Vector3d tangent(const Vector3d &v1, const Vector3d &vc,
                            const Vector3d &v2, double x);

    /// @brief Нормаль к кривой в точке x, лежит в одной плоскости с кривой
    /// и вместе с точкой 'c' (если возможно), направлена от точки 'c'.
    Vector3d normal(double x, const Vector3d& c) const;

    /// @brief Якобиан отображения в точке x. Фактически модуль производной.
    double Jacobian(double x) const;

    /// @brief Центральный узел криволинейного отрезка,
    /// соответствует параметру x = 0.
    const Vector3d& center() const;

    /// @brief Нормализованный криволинейный интеграл второго рода
    /// int vec(n) dl, где выбирается нормаль, которая лежит в одной
    /// плоскости  вместе с кривой и вместе с точкой 'c' (если возможно),
    /// направлена от точки 'c'.
    Vector3d normal(const Vector3d& c) const;

    /// @brief "Длина" криволинейной одномерной грани. Равна |int vec(n) dl|,
    /// то есть модулю криволинейного интеграла 2-го рода.
    double length() const;

protected:
    /// @brief Возвращает проецию вектора 'a' на плоскость, в которой лежит
    /// кривая. Если все точки кривой лежат на одной прямой, то возвращает
    /// сам вектор 'a' без изменений.
    Vector3d projection(const Vector3d& a) const;
};

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

    std::array<SqQuad, 4> children() const;
};

/// @brief Представление шестигранника
struct Cube {
protected:
    /// @brief Вешины (таблица 2 x 2 x 2)
    std::array<Vector3d, 8> verts;

public:
    /// @brief Конструктор по умолчанию
    Cube() = default;

    /// @brief Конструктор по угловым точкам
    Cube(const Vector3d &v000,
         const Vector3d &v001,
         const Vector3d &v010,
         const Vector3d &v011,
         const Vector3d &v100,
         const Vector3d &v101,
         const Vector3d &v110,
         const Vector3d &v111);

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..7]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..7]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индексы {-1, +1}^3 в [0..7]
    template <int i, int j, int k>
    inline static constexpr int iss() {
        static_assert(i * i == 1 && j * j == 1 && k * k == 1,
                      "Available indices: {-1, +1}");
        return 2*(k + 1) + j + 1 + (i + 1) / 2;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, +1}
    /// @tparam j in {-1, +1}
    /// @tparam k in {-1, +1}
    template <int i, int j, int k>
    inline Vector3d &vs() {
        return verts[iss<i, j, k>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, +1}
    /// @tparam j in {-1, +1}
    /// @tparam k in {-1, +1}
    template <int i, int j, int k>
    inline const Vector3d &vs() const {
        return verts[iss<i, j, k>()];
    }

    /// @brief Непосредственно отображение
    Vector3d operator()(double x, double y, double z) const;

    /// @brief Непосредственно отображение
    Vector3d get(double x, double y, double z) const;

    Vector3d tangent_x(double x, double y, double z) const;

    Vector3d tangent_y(double x, double y, double z) const;

    Vector3d tangent_z(double x, double y, double z) const;

    /// @brief Якобиан отображения
    double Jacobian(double x, double y, double z) const;

    /// @brief Центр трехмерной ячейки
    Vector3d center() const;

    /// @brief Объем простой трехмерной ячейки
    double volume() const;

    /// @brief Барицентр трехмерной ячейки
    /// @param volume Объем ячейки
    Vector3d centroid(double volume = 0.0) const;
};

/// @brief Представление квадратичного шестигранника
struct SqCube {
protected:
    /// @brief Вешины (таблица 3 x 3 x 3)
    std::array<Vector3d, 27> verts;

public:
    /// @brief Конструктор по умолчанию
    SqCube() = default;

    /// @brief Конструктор по угловым точкам
    SqCube(const Vector3d &v000, const Vector3d &v002,
           const Vector3d &v020, const Vector3d &v022,
           const Vector3d &v200, const Vector3d &v202,
           const Vector3d &v220, const Vector3d &v222);

    /// @brief Конструктор по полному набору узлов
    SqCube(const Vector3d &v000, const Vector3d &v001, const Vector3d &v002,
           const Vector3d &v010, const Vector3d &v011, const Vector3d &v012,
           const Vector3d &v020, const Vector3d &v021, const Vector3d &v022,
           const Vector3d &v100, const Vector3d &v101, const Vector3d &v102,
           const Vector3d &v110, const Vector3d &v111, const Vector3d &v112,
           const Vector3d &v120, const Vector3d &v121, const Vector3d &v122,
           const Vector3d &v200, const Vector3d &v201, const Vector3d &v202,
           const Vector3d &v210, const Vector3d &v211, const Vector3d &v212,
           const Vector3d &v220, const Vector3d &v221, const Vector3d &v222);

    /// @brief Конструктор по угловым точкам
    SqCube(const Cube &cube);

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    SqCube(const Quad& quad);

    /// @brief Неполный конструктор. Инициализация только первого слоя
    /// вершин для хранения представления двумерной ячейки.
    SqCube(const SqQuad& quad);

    /// @brief Интерпретировать первый слой как двумерное отображение
    inline SqQuad& as2D() {
        return *reinterpret_cast<SqQuad *>(verts.data());
    };

    /// @brief Интерпретировать первый слой как двумерное отображение
    inline const SqQuad& as2D() const {
        return *reinterpret_cast<const SqQuad*>(verts.data());
    };

    /// @brief Удалить промежуточные вершины и вернуть простой кубоид
    Cube reduce() const;

    /// @brief Разбиение на 8 кубов
    std::array<SqCube, 8> children() const;

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..27)
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..27)
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображает индексы {-1, 0, +1}^3 в [0..27)
    template <int i, int j, int k>
    inline static constexpr int iss() {
        static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1,
                      "Available indices: {-1, 0, +1}");
        return 9*(k + 1) + 3*(j + 1) + i + 1;
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 0, +1}
    /// @tparam j in {-1, 0, +1}
    /// @tparam k in {-1, 0, +1}
    template <int i, int j, int k = -1>
    inline Vector3d &vs() {
        return verts[iss<i, j, k>()];
    }

    /// @brief Оператор доступа по индексам отображения
    /// @tparam i in {-1, 0, +1}
    /// @tparam j in {-1, 0, +1}
    /// @tparam k in {-1, 0, +1}
    template <int i, int j, int k = -1>
    inline const Vector3d &vs() const {
        return verts[iss<i, j, k>()];
    }

    /// @brief Непосредственно отображение
    Vector3d operator()(double x, double y, double z) const;

    /// @brief Неполное отображение (для двумерных ячеек)
    /// Эквивалентно заданию z = -1.
    Vector3d operator()(double x, double y) const;

    /// @brief Непосредственно отображение
    Vector3d get(double x, double y, double z) const;

    Vector3d tangent_x(double x, double y, double z) const;

    Vector3d tangent_y(double x, double y, double z) const;

    Vector3d tangent_z(double x, double y, double z) const;

    /// @brief Якобиан отображения
    double Jacobian(double x, double y, double z) const;

    /// @brief Центр трехмерной ячейки
    Vector3d center() const;

    /// @brief Объем простой трехмерной ячейки
    double volume() const;

    /// @brief Барицентр трехмерной ячейки
    /// @param volume Объем ячейки
    Vector3d centroid(double volume = 0.0) const;
};

} // namespace zephyr::geom