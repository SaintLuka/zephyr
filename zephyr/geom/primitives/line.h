#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @addtogroup Geom-Primitives
/// @brief Геометрические примитивы
/// @{

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

    /// @brief Для отрезков эквивалентно center()
    Vector3d centroid() const;

    /// @brief Барицентр для осесимметричных задач, вращение вокруг оси x.
    Vector3d centroid(bool axial) const;

    /// @brief Нормаль к отрезку, располагается в плоскости с отрезком
    /// и точкой 'c', направлена от точки 'c'.
    Vector3d normal(const Vector3d &c) const;

    /// @brief Длина отрезка
    double length() const;

    /// @brief Площадь в случае осевой симметрии
    double area_as() const;

    /// @brief Нормаль, умноженная на длину отрезка, нормаль располагается
    /// в плоскости с отрезком и точкой 'c', направлена от точки 'c'.
    Vector3d area_n(const Vector3d &c) const;
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

/// @}

} // namespace zephyr::geom