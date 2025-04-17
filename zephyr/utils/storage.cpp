#include <zephyr/utils/storage.h>

namespace zephyr::utils {

void Storage::resize(size_t new_size) {
    m_size = new_size;
    resize_basic(new_size);

    for (int i = 0; i < m_values2.size(); ++i) {
        for (auto& vec: m_values2[i]) {
            vec.resize(new_size * custom_type_sizeof(i));
        }
    }
}

void Storage::print() const {
    if (empty()) {
        std::cout << "Storage is empty\n";
    }
    else {
        std::cout << "Storage contains " << size() << " elements\n";
    }
    int count = print_basic();

    for (int i = 0; i < n_custom_types; ++i) {
        if (!m_values2[i].empty()) {
            std::cout << "  " << custom_type_sizeof(i) << " byte type arrays:\n    ";
            for (int j = 0; j < m_names2[i].size() - 1; ++j) {
                std::cout << "'" << m_names2[i][j] << "', ";
            }
            std::cout << "'" << m_names2[i].back() << "'\n";
            ++count;
        }
    }

    if (count < 1) {
        std::cout << "  Has no fields yet\n";
    }
    std::cout << "\n";
}

} // namespace zephyr::utils