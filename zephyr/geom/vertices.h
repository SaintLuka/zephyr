#pragma once

#include <array>

#include <zephyr/geom/vector.h>

namespace zephyr { namespace geom {

/// @brief Трехмерный вектор-столбец из Eigen
/// с несколькими дополнительными функциями
struct Vertex : public Vector3d {

    template <class... Args>
    Vertex(Args&&... args) : Vector3d(std::forward<Args>(args)...) { }

    template<typename OtherDerived>
    Vertex(const Eigen::MatrixBase<OtherDerived> &other) {
        *reinterpret_cast<zephyr::geom::Vector3d *>(this) = other;
    }

    template<typename OtherDerived>
    Vertex &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        *reinterpret_cast<zephyr::geom::Vector3d *>(this) = other;
        return *this;
    }

    inline bool is_undefined() const {
        return x() != x();
    }

    inline bool is_actual() const {
        return x() == x();
    }

    inline void set_undefined() {
        x() = y() = z() = 0.0 / 0.0;
    }
};

/// @brief Координаты вершин (узлов) ячейки.
/// @details Массив фиксированной длины до max_size элементов. Класс содержит
/// вершины эйлеровой ячейки и может описывать один из следующих примитивов:
//
///    2D элементы:
///    1. Полигоны с числом вершин от 3 до 8. Вершины располагаются в списке
///    против часовой стрелки.
///    2. AMR четырехугольник. Список вершин содержит 9 вершин в виде двумерной
///    таблицы.
///    3D элементы:
///    1. Тетраэдр (4 вершины). Пирамида (5 вершин). Клин (6 вершин).
///    Шестигранник (8 вершин). Нумерация вершин в соответствии с форматом VTK.
///    2. AMR шестигранник. Список вершин содержит 27 вершин в виде трехмерной
///    таблицы.
///
///    Из определения выше следует, что размерность элемента и количество
///    вершин однозначно определяют тип элемента, что в частности используется
///    при сохранении элементов в формате VTK.
class Vertices {
public:

    static const int max_size = 27;      ///< Максимальное число вешин

    /// @brief Конструктор по умолчанию
    Vertices() = default;

    /// @brief Число актуальных вершин.
    int size() const {
        int count = 0;
        for (int i = 0; i < max_size; ++i) {
            if (!list[i].is_undefined()) {
                ++count;
            }
        }
        return count;
    }

    /// @brief Установить неопределенный список
    void set_undefined() {
        for (auto &c: list) {
            c.set_undefined();
        }
    }

    /// @brief Доступ к вершине по индексу
    inline Vertex &operator[](int i) {
        return list[i];
    }

    /// @brief Доступ к вершине по индексу
    inline const Vertex &operator[](int i) const {
        return list[i];
    }

private:
    /// @brief Массив вершин ячейки (до max_size)
    std::array<Vertex, max_size> list;
};

/*
    Содержит набор из нескольких функций, упрощающих индексацию
    в четырехугольных ячейках и шестигранниках, особенно при использовании AMR.
    Отображения строятся из множества трехмерных индексов на одномерный индекс.
    Коротким (s/short) индексом будем называть трехмерный индекс из
    множества {0, 1}^3 или одномерный индекс из [0, 8).
    Длинным (w/wide) индексом будем называть трехмерный индекс из
    множества {0, 1, 2}^3 или одномерный индекс из [0, 27).
 */
namespace topology {

/// @brief Отображение {0, 1}^3 -> [0, 8)
/// @return Индекс [0, 8) вершины ячейки.
inline constexpr int iss(int i, int j, int k = 0) {
    return 4 * k + 2 * j + i;
}

/// @brief Отображение {0, 2}^3 -> [0, 8)
/// @return Индекс [0, 8) вершины ячейки.
inline constexpr int iws(int i, int j, int k = 0) {
    return 2 * k + j + i / 2;
}

/// @brief Отображение {0, 1}^3 -> [0, 27)
/// @return Индекс [0, 27) вершины ячейки.
inline constexpr int isw(int i, int j, int k = 0) {
    return 2 * (9 * k + 3 * j + i);
}

/// @brief Отображение {0, 1, 2}^3 -> [0, 27)
/// @return Индекс [0, 27) вершины ячейки.
inline constexpr int iww(int i, int j, int k = 0) {
    return 9 * k + 3 * j + i;
}

/// @brief Динамическое отображение. Отображает короткий индекс {0, 1}^3
/// на одномерный индекс, соответствующий type::_vertices_::max_size
inline constexpr int iv(int i, int j, int k = 0) {
    return isw(i, j, k);
}

} // topology
} // geom
} // zephyr