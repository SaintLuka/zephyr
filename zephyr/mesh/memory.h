#pragma once

#include <vector>

namespace zephyr::mesh {

/// @brief Измеряет расход памяти
struct memory_t {
    size_t needed{0};
    size_t actual{0};

    template <typename T>
    void add(const std::vector<T>& arr) {
        needed += sizeof(T) * arr.size();
        actual += sizeof(T) * arr.capacity();
    }
};

}