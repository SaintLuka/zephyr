#pragma once

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

} // geom
} // zephyr