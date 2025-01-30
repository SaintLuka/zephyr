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

void Materials::append(Eos::Ref eos) {
    Material mat(eos);
    append(mat);
}

void Materials::operator+=(Eos::Ref eos) {
    append(eos);
}

MixturePT Materials::mixture_PT() const {
    MixturePT mix;
    for (const auto& mat: m_materials) {
        mix += mat.eos();
    }
    return mix;
}

} // namespace zephyr::phys