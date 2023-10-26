#pragma once

namespace zephyr::geom {

/// @class Базовый класс для геометрических данных хранилища AmrStorage.
/// На данный момент от него наследуются следующие классы:
///     AmrCell  -- Квадратная/кубическая ячейка адаптивной сетки
///     PolyCell -- Многоугольная/многогранная ячейка подвижной сетки
///     PolyNode -- Узел подвижной сетки
class Element {
public:
    int rank;   ///< Ранг процесса (= -1 для неактивных элементов)
    int index;  ///< Индекс ячейки в AmrStorage (>= 0)
    int next;   ///< Новый индекс (в алгоритмах с перестановками)

    /// @brief Конструктор по умолчанию
    Element(int r = -1, int idx = 0) :
        rank(r), index(idx), next(-1) { }

    /// @brief Актуальная ячейка? (rank >= 0)
    inline bool is_actual() const {
        return rank >= 0;
    }

    /// @brief Ячейка к удалению (rank < 0)
    inline bool is_undefined() const {
        return rank < 0;
    }

    /// @brief Устанавливает rank = -1 (ячейка вне сетки)
    inline void set_undefined(){
        rank = -1;
    }
};

} // namespace zephyr::geom