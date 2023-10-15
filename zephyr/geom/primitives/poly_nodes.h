#pragma once

#include <zephyr/geom/primitives/opt_vertex.h>

namespace zephyr::geom {

/// @brief Координаты вершин (узлов) ячейки.
/// @details Массив фиксированной длины до max_size элементов. Класс содержит
/// узлы подвижной ячейки и может описывать один из следующих примитивов:
/// 2D элементы:
///    Полигоны с числом вершин от 3 до 8. Вершины располагаются в списке
///    против часовой стрелки.
/// 3D элементы:
///    Тетраэдр (4 вершины). Пирамида (5 вершин). Клин (6 вершин).
///    Шестигранник (8 вершин). Нумерация вершин в соответствии с форматом VTK.
///    Произвольный многогранник.
class PolyNodes {
public:

    static const int max_size = 27;      ///< Максимальное число вешин

    /// @brief Конструктор по умолчанию
    PolyNodes() {
        set_undefined();
    }

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
    inline OptVertex &operator[](int i) {
        return list[i];
    }

    /// @brief Доступ к вершине по индексу
    inline const OptVertex &operator[](int i) const {
        return list[i];
    }

private:
    /// @brief Массив вершин ячейки (до max_size)
    std::array<OptVertex, max_size> list;
};

} // namespace zephyr::geom