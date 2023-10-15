#pragma once

#include <zephyr/geom/geom.h>

namespace zephyr::geom {

/// @class Базовый класс для геометрических данных хранилища Storage.
/// На данный момент от него наследуются следующие классы:
///     AmrCell  -- квадратная ячейка адаптивной сетки
///     PolyCell -- многоугольная/многогранная ячейка подвижной сетки
///     PolyNode -- Узел подвижной сетки
class Element {
public:
    int rank;   ///< Ранг процесса (= -1 для неактивных элементов)
    int index;  ///< Индекс ячейки в Storage (>= 0)
    int next;   ///< Новый индекс (в алгоритмах с перестановками)

    /// @brief Конструктор по умолчанию
    Element(int r = -1, int idx = 0) :
        rank(r), index(idx) { }
};

} // namespace zephyr::geom