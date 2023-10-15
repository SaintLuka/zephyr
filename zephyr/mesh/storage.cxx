#include <zephyr/mesh/storage.h>

namespace zephyr { namespace mesh {

Storage::Storage(Storage& src)
    : m_size(src.m_size),
      m_itemsize(src.m_itemsize),
      m_data(src.m_data.size()) {

    std::memcpy(m_data.data(), src.m_data.data(),
                m_data.size() * sizeof(Byte));
}

Storage::Storage(Storage&& src)
        : m_size(src.m_size),
          m_itemsize(src.m_itemsize),
          m_data(std::move(src.m_data)) {
}

void Storage::init(int datasize, int size) {
    m_size = std::max(0, size);
    m_itemsize = int(sizeof(AmrCell)) + std::max(0, datasize);
    m_data.resize(m_size * m_itemsize);
}


bool Storage::empty() const {
    return m_size < 1;
}

int Storage::size() const {
    return m_size;
}

int Storage::datasize() const {
    return m_itemsize - int(sizeof(AmrCell));
}

int Storage::itemsize() const {
    return m_itemsize;
}

void Storage::resize(int new_size) {
    m_size = std::max(0, new_size);
    m_data.resize(m_size * m_itemsize);
}


Storage::Item Storage::operator[](int i) {
    return { m_data.data() + m_itemsize * i, m_itemsize };
}

Storage::Item Storage::begin() {
    return { m_data.data(), m_itemsize };
}

Storage::Item Storage::end() {
    return { m_data.data() + m_itemsize * m_size, m_itemsize };
}



} // mesh
} // zephyr
