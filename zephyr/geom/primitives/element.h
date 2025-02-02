#pragma once

namespace zephyr::geom {

/// @class Базовый класс для геометрических данных хранилища AmrStorage.
/// На данный момент от него наследуются следующие классы:
///     AmrCell -- Квадратная/кубическая ячейка адаптивной сетки, а также
///                многоугольная/многогранная ячейка эйлеровой сетки
///     MovCell -- Многоугольная/многогранная ячейка подвижной сетки
///     MovNode -- Узел подвижной сетки
struct Element {
    int rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    int index;  ///< Индекс элемента в локальном Storage (< 0 для неактивных,
                ///  неопределенных элементов, элементов на удаление)
    int next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)

    /// @brief Конструктор по умолчанию
    explicit Element(int r = -1, int idx = -1) :
        rank(r), index(idx), next(-1) { }

    /// @brief Актуальная ячейка?
    inline bool is_actual() const {
        return index >= 0;
    }

    /// @brief Ячейка к удалению
    inline bool is_undefined() const {
        return index < 0;
    }

    /// @brief Устанавливает index = -1 (ячейка вне сетки)
    inline void set_undefined(){
        index = -1;
    }
};

} // namespace zephyr::geom