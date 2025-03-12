#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/quad.h>

namespace zephyr::geom {

/// @addtogroup Geom-Primitives
/// @{

/// @struct Cube cube.h
/// @brief Представление шестигранника (топологический куб)
/// @details Линейное отображение куба [-1, 1]^3 на произвольный
/// топологический кубоид.
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

/// @struct SqCube cube.h
/// @brief Представление квадратичного шестигранника (топологический куб)
/// @details Квадратичное отображение куба [-1, 1]^3 на произвольный
/// топологический кубоид.
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

/// @}

} // namespace zephyr::geom