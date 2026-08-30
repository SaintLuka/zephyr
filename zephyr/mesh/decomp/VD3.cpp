#include <zephyr/mesh/decomp/VD3.h>
#include <zephyr/math/random.h>

namespace zephyr::mesh::decomp {

VD3::VD3(const Box &domain, int size)
    : Decomposition(size) {
    std::vector<Vector3d> gs(size);

    auto gen = math::Random::Uniform(domain);
    for (int i = 0; i < size; ++i) {
        gs[i] = gen->next();
    }
    diagram_ = VDiagram(domain, gs);
}

int VD3::rank(const EuCell &elem) const {
    return diagram_.rank(elem.center());
}

void VD3::balancing(const std::vector<double> &w) {
    diagram_.balancing(w);
}

} // namespace zephyr::mesh::decomp