#include <map>

#include <zephyr/utils/mpi.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/cell.h>
#include <zephyr/mesh/mesh.h>

namespace zephyr::mesh {

using namespace zephyr::utils;

void Mesh::set_decomposition(Decomposition::Ref decmp, bool update) {
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

void Mesh::set_decomposition(ORB& orb, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    ORB::Ptr decmp = std::make_shared<ORB>(orb);
    orb = *decmp;

    set_decomposition(decmp, update);
#endif
}

void Mesh::set_decomposition(const std::string& type, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    auto domain = bbox();
    decomp_ = ORB::create(domain, type, mpi::size());

    set_decomposition(decomp_, update);
#endif
}

RawNodes& Mesh::ghost_nodes() {
#ifndef ZEPHYR_MPI
    static RawNodes ghost_nodes;
    return ghost_nodes;
#else
    return tourists_.ghost_nodes();
#endif
}

const RawNodes& Mesh::ghost_nodes() const {
#ifndef ZEPHYR_MPI
    static RawNodes ghost_nodes;
    return ghost_nodes;
#else
    return tourists_.ghost_nodes();
#endif
}

void Mesh::balancing() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    bool done = decomp_->exact_balancing(local_cells_.center);
    if (done) return;

    double load = local_cells_.n_cells();
    balancing(load);
#endif
}

void Mesh::balancing(double load) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    bool done = decomp_->exact_balancing(local_cells_.center);
    if (done) return;

    auto ws = mpi::all_gather(load);
    decomp_->balancing(ws);
#endif
}

void Mesh::prebalancing(int n_iters) {
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

void Mesh::setup_ranks() {
    // Определим новый rank для всех ячеек из locals
    for_each([decomp=decomp_](const Cell &cell) {
        cell.set_rank(decomp->rank(cell));
    });
}

} // namespace zephyr::mesh
