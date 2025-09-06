#include <iostream>

#include <zephyr/mesh/storage.h>

namespace zephyr::mesh {

Storage Storage::same() const {
    Storage new_s;
    new_s.m_size = 0;
    new_s.m_names = m_names;
    for (auto& buf: m_data) {
        new_s.m_data.emplace_back(buf.same());
    }
    return new_s;
}

void Storage::resize(size_t new_size) {
    m_size = new_size;
    for (auto& buf: m_data) {
        buf.resize(new_size);
    }
}

void Storage::reserve(size_t new_size) {
    for (auto& buf: m_data) {
        buf.reserve(new_size);
    }
}

void Storage::shrink_to_fit() {
    for (auto& buf: m_data) {
        buf.shrink_to_fit();
    }
}

void Storage::copy_data(size_t from, size_t to) {
    for (auto& buf: m_data) {
        buf.copy_data(from, to);
    }
}

void Storage::copy_data(size_t from, Storage* dst, size_t to) const {
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i].copy_data(from, dst->m_data[i], to);
    }
}

} // namespace zephyr::mesh