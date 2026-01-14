#pragma once

#include <ranges>

/// @brief Основное пространство имен.
namespace zephyr {

template <std::integral I>
using range_t = decltype(std::views::iota(I{}, I{}));

/// @brief Диапазон чисел как в python
template <std::integral Index = int>
range_t<Index> range(Index stop) {
    return std::views::iota(Index{0}, stop);
}

template <std::integral Index = int>
range_t<Index> range(Index start, Index stop) {
    return std::views::iota(start, stop);
}

/// @brief Функция enumerate как в python
template <std::ranges::range Range>
auto enumerate(Range&& range) {
    return std::forward<Range>(range)
           | std::views::transform([](auto&& item) {
               static std::size_t index = 0;
               return std::make_pair(index++, std::forward<decltype(item)>(item));
           });
}

} // namespace zephyr::utils