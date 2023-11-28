#pragma once

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
class BNodes {
public:
    ///@brief Максимальное число вешин
    static const int max_count = 27;

    BNodes() {
        m_nodes.fill(-1);
    }

    /// @brief Доступ к вершине по индексу
    inline int& operator[](int i) {
        assert(i < max_count);
        return m_nodes[i];
    }

    /// @brief Доступ к вершине по индексу
    inline int operator[](int i) const {
        assert(i < max_count);
        return m_nodes[i];
    }

    /*
    struct iterator {
    private:
        int *m_ptr;

        iterator(int *ptr) : m_ptr(ptr) {}

    };
     */

    int count() const {
        int count = 0;
        for (int idx: m_nodes) {
            if (idx >= 0) {
                ++count;
            }
        }
        return count;
    }

private:
    /// @brief Массив индексов вершин ячейки
    std::array<int, max_count> m_nodes;
};

} // namespace zephyr::geom