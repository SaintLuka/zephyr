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

    decomp_ = decmp;
    if (update) {
        // вызываю, чтобы инициализировать tourists_
        tourists_.update(local_cells_, local_nodes_);
        redistribute();

        // Вероятно, первый (и единственный) redistribute, почистим память
        local_cells_.shrink_to_fit();

        tourists_.shrink_to_fit();

        migrants_.clear();
        migrants_.shrink_to_fit();
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
    decomp_ = ORB::create(domain, type, mpi::size());

    set_decomposition(decomp_, update);
#endif
}

const AmrNodes& EuMesh::ghost_nodes() const {
#ifndef ZEPHYR_MPI
    static AmrNodes ghost_nodes;
    return ghost_nodes;
#else
    return tourists_.ghost_nodes();
#endif
}

void EuMesh::balancing() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    bool done = decomp_->exact_balancing(local_cells_.center);
    if (done) return;

    double load = local_cells_.size();
    balancing(load);
#endif
}

void EuMesh::balancing(double load) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    bool done = decomp_->exact_balancing(local_cells_.center);
    if (done) return;

    auto ws = mpi::all_gather(load);
    decomp_->balancing(ws);
#endif
}

void EuMesh::prebalancing(int n_iters) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    bool done = decomp_->exact_balancing(local_cells_.center);
    if (done) {
        redistribute();
        return;
    }
    for (int i = 0; i < n_iters; ++i) {
        balancing();
        redistribute();
    }
#endif
}

void EuMesh::setup_ranks() {
    // Определим новый rank для всех ячеек из locals
    for_each([decomp=decomp_](const EuCell &cell) {
        cell.set_rank(decomp->rank(cell));
    });
}

} // namespace zephyr::mesh
