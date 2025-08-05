#pragma once

#include <functional>

namespace zephyr::mesh {

class EuCell;
class Children;

/// @brief Функции для огрубления и распределения (физических) данных в ячейках
/// при адаптации. Задаются пользователем.
/// @ingroup euler-mesh
struct Distributor {
    /// @brief Тип функции, распределяющей данные между дочерними ячейками
    using split_function = std::function<void(EuCell&, Children&)>;

    /// @brief Тип функции, объединяющей данные дочерних ячеек
    using merge_function = std::function<void(Children&, EuCell&)>;

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
    /// в родительсую ячейку из первой дочерней.
    static Distributor simple();
};

} // namespace zephyr::mesh