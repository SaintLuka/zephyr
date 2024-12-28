#include <zephyr/phys/matter/materials.h>

namespace zephyr::phys {


void Materials::clear() {
    m_materials.clear();
}

void Materials::append(Material &mat) {
    m_materials.push_back(mat);
}

void Materials::operator+=(Material &mat) {
    append(mat);
}


} // namespace zephyr::phys