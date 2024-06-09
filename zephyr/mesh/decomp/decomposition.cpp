#include <numeric>

#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/decomp/decomposition.h>

namespace zephyr::mesh::decomp {

using namespace zephyr::utils;

Decomposition::Decomposition() {
    m_size = mpi::size();
}

Decomposition::Decomposition(int size) {
    m_size = std::max(1, size);
}

double Decomposition::imbalance(const std::vector<double>& ws) {
    return *std::max_element(ws.begin(), ws.end()) * ws.size() / std::accumulate(ws.begin(), ws.end(), 0.0)  - 1.0;
}

} // namespace zephyr::mesh::decomp