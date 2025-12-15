#include <map>

#include <zephyr/utils/mpi.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::mesh {

using namespace zephyr::utils;

void EuMesh::set_decomposition(Decomposition::Ref decmp, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    m_decomp = decmp;
    if (update) {
        // вызываю, чтобы инициализировать m_tourists
        m_tourists.update(m_locals);
        redistribute();

        // Вероятно, первый (и единственный) redistribute, почистим память
        m_locals.shrink_to_fit();

        m_tourists.shrink_to_fit();

        m_migrants.clear();
        m_migrants.shrink_to_fit();
    }
#endif
}

void EuMesh::set_decomposition(ORB& orb, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    ORB::Ptr decmp = std::make_shared<ORB>(orb);
    orb = *decmp;

    set_decomposition(decmp, update);
#endif
}

void EuMesh::set_decomposition(const std::string& type, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    auto domain = bbox();
    m_decomp = ORB::create(domain, type, mpi::size());

    set_decomposition(m_decomp, update);
#endif
}

void EuMesh::balancing() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    double load = m_locals.size();
    balancing(load);
#endif
}

void EuMesh::balancing(double load){
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    auto ws = mpi::all_gather(load);
    m_decomp->balancing(ws);
#endif
}

void EuMesh::prebalancing(int n_iters) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    if (mpi::master()) {
        std::cerr << "Warning: prebalancing is not optimal\n";
    }

    for (int i = 0; i < n_iters; ++i) {
        balancing();
        redistribute();
    }
#endif
}

void EuMesh::setup_ranks() {
    // Определим новый rank для всех ячеек из locals
    for_each([this](EuCell cell) {
        cell.set_rank(m_decomp->rank(cell));
    });
}

} // namespace zephyr::mesh
