#include <zephyr/utils/storage.h>

namespace zephyr::utils {

void SoaStorage::resize(size_t new_size) {
    m_size = new_size;
    resize_basic(new_size);

    for (int i = 0; i < m_values2.size(); ++i) {
        for (auto& vec: m_values2[i]) {
            vec.resize(new_size * custom_type_sizeof(i));
        }
    }
}

void SoaStorage::reserve(size_t new_size) {
    m_size = new_size;
    reserve_basic(new_size);

    for (int i = 0; i < m_values2.size(); ++i) {
        for (auto& vec: m_values2[i]) {
            vec.reserve(new_size * custom_type_sizeof(i));
        }
    }
}

void SoaStorage::print() const {
    if (empty()) {
        std::cout << "SoaStorage is empty\n";
    }
    else {
        std::cout << "SoaStorage contains " << size() << " elements\n";
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


void SoaStorage::copy_data(size_t from, size_t to) {
    copy_data(from, this, to);
}

void SoaStorage::copy_data(size_t from, SoaStorage* dst, size_t to) {
    copy_basics(m_values1, from, dst->m_values1, to);

    for (size_t i = 0; i < n_custom_types; ++i) {
        if (!m_values2[i].empty() && m_values2[i].size() == dst->m_values2[i].size()) {
            size_t datasize = custom_type_sizeof(i);
            size_t count = m_values2[i].size();
            for (size_t k = 0; k < count; ++k) {
                std::memcpy(
                        dst->m_values2[i][k].data() + datasize * to,
                        m_values2[i][k].data() + datasize * from,
                        datasize);
            }
        }
    }
}

} // namespace zephyr::utils