#pragma once

#include <functional>

#include <zephyr/mesh/storage.h>

namespace zephyr { namespace mesh {

/// @brief Класс содержит набор функций для огрубления и распределения
// (физических) данных в ячейках при адаптации. Определяются пользователем.
struct Distributor {
    using split_function_2D = std::function<void(Storage::Item, const std::array<Storage::Item, 4> &)>;
    using split_function_3D = std::function<void(Storage::Item, const std::array<Storage::Item, 8> &)>;
    using merge_function_2D = std::function<void(const std::array<Storage::Item, 4> &, Storage::Item)>;
    using merge_function_3D = std::function<void(const std::array<Storage::Item, 8> &, Storage::Item)>;

    split_function_2D split2D;
    split_function_3D split3D;

    merge_function_2D merge2D;
    merge_function_3D merge3D;

    /// @brief Конструктор по умолчанию. Определяет Distributor, который
    /// вообще ничего не делает
    Distributor();

    /// @brief Создает Distributor, который вообще ничего не делает.
    static Distributor empty();

    /// @brief Создает Distributor, который определяет функции split2D, split3D,
    /// merge2D, merge3D простейшим способом. Для функций split используется
    /// простой перенос данных в дочерние ячейки. Для функций merge
    /// используется перенос данных в родительсую ячейку из первой дочерней.
    static Distributor simple();

    /// @brief Шаблонная функция split
    template <int dim>
    inline void split(Storage::Item parent, const std::array<Storage::Item, 4*(dim-1)>& children) const;

    /// @brief Шаблонная функция merge
    template <int dim>
    inline void merge(const std::array<Storage::Item, 4*(dim-1)>& children, Storage::Item parent) const;

};

template <>
inline void Distributor::split<2>(Storage::Item parent,
        const std::array<Storage::Item, 4>& children) const {
    return split2D(parent, children);
}

template <>
inline void Distributor::split<3>(Storage::Item parent,
        const std::array<Storage::Item, 8>& children) const {
    return split3D(parent, children);
}

template <>
inline void Distributor::merge<2>(const std::array<Storage::Item, 4>& children,
        Storage::Item parent) const {
    return merge2D(children, parent);
}

template <>
inline void Distributor::merge<3>(const std::array<Storage::Item, 8>& children,
        Storage::Item parent) const {
    return merge3D(children, parent);
}

} // namespace mesh
} // namespace zephyr