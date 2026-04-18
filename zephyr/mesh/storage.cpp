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

const utils::Buffer& Storage::operator[](const std::string& name) const {
    for (size_t i = 0; i < m_names.size(); ++i) {
        if (m_names[i] == name) {
            return m_data[i];
        }
    }
    throw std::runtime_error("Does not have a buffer named " + name);
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

void Storage::print_names() const {
    if (m_names.empty()) { return; }
    std::cout << "Buffer names: ";
    for (size_t i = 0; i < m_names.size() - 1; ++i) {
        std::cout << m_names[i] << ", ";
    }
    std::cout << m_names.back() << "\n";
}

} // namespace zephyr::mesh