#include <zephyr/mesh/decomp/VD3.h>


namespace zephyr::mesh::decomp {

VD3::VD3(const Box &domain, int size)
    : Decomposition(size) {
    std::vector<Vector3d> gs(size);

    auto gen = domain.random2D(0);
    for (int i = 0; i < size; ++i) {
        gs[i] = gen.get();
    }

    m_diagram = VDiagram(domain, gs);
}

int VD3::rank(QCell &elem) const {
    return m_diagram.rank(elem.center());
}

int VD3::rank(AmrStorage::Item &elem) const {
    return m_diagram.rank(elem.center);
}

void VD3::balancing(const std::vector<double> &w) {
    m_diagram.balancing(w);
}

} // namespace zephyr::mesh::decomp