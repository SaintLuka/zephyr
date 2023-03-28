#pragma once

#include <functional>

#include <zephyr/data/storage.h>

namespace zephyr { namespace mesh {

using zephyr::data::Storage;

/**
	\brief
        \~russian Класс содержит набор функций для огрубления и распределения
        (физических) данных в ячейках при адаптации. Определяются пользователем.
		\~english .
		\~
*/
struct DataDistributor {
    using split_function_2D = std::function<void(Storage::iterator, const std::array<Storage::iterator, 4> &)>;
    using split_function_3D = std::function<void(Storage::iterator, const std::array<Storage::iterator, 8> &)>;
    using merge_function_2D = std::function<void(const std::array<Storage::iterator, 4> &, Storage::iterator)>;
    using merge_function_3D = std::function<void(const std::array<Storage::iterator, 8> &, Storage::iterator)>;

    split_function_2D split2D = [](Storage::iterator, const std::array<Storage::iterator, 4> &) { };
    split_function_3D split3D = [](Storage::iterator, const std::array<Storage::iterator, 8> &) { };

    merge_function_2D merge2D = [](const std::array<Storage::iterator, 4> &, Storage::iterator) { };
    merge_function_3D merge3D = [](const std::array<Storage::iterator, 8> &, Storage::iterator) { };

    template <unsigned int dim>
    inline void split(Storage::iterator parent, const std::array<Storage::iterator, 4*(dim-1)>& children) const;


    template <unsigned int dim>
    inline void merge(const std::array<Storage::iterator, 4*(dim-1)>& children, Storage::iterator parent) const;

};

template <>
inline void DataDistributor::split<2>(Storage::iterator parent,
        const std::array<Storage::iterator, 4>& children) const {
    return split2D(parent, children);
}

template <>
inline void DataDistributor::split<3>(Storage::iterator parent,
        const std::array<Storage::iterator, 8>& children) const {
    return split3D(parent, children);
}

template <>
inline void DataDistributor::merge<2>(const std::array<Storage::iterator, 4>& children,
        Storage::iterator parent) const {
    return merge2D(children, parent);
}

template <>
inline void DataDistributor::merge<3>(const std::array<Storage::iterator, 8>& children,
        Storage::iterator parent) const {
    return merge3D(children, parent);
}

} // namespace mesh
} // namespace zephyr