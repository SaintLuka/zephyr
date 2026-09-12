#pragma once

#include <functional>

namespace zephyr::mesh {

class Cell;
class Children;

/// @brief Функции для огрубления и распределения (физических) данных в ячейках
/// при адаптации. Задаются пользователем.
/// @ingroup raw-mesh
struct Distributor {
    /// @brief Тип функции, распределяющей данные между дочерними ячейками
    using split_function = std::function<void(const Cell&, Children&)>;

    /// @brief Тип функции, объединяющей данные дочерних ячеек
    using merge_function = std::function<void(const Children&, Cell&)>;

    split_function split;  ///< Распределение данным между дочерними
    merge_function merge;  ///< Объединение данных дочерних ячеек

    /// @brief Конструктор по умолчанию. Определяет Distributor, который
    /// вообще ничего не делает
    Distributor();

    /// @brief Создает Distributor, который вообще ничего не делает.
    static Distributor empty();

    /// @brief Создает Distributor, который определяет функции split и merge
    /// простейшим способом. Для функции split используется простой перенос
    /// данных в дочерние ячейки. Для функции merge используется перенос данных
    /// в родительскую ячейку из первой дочерней.
    static Distributor simple();

    /// @brief Создает дистрибутор, который использует функцию func для
    /// инициализации ячеек. Удобно для задания начальных данных
    static Distributor initializer(std::function<void(Cell&)> func);
};

} // namespace zephyr::mesh