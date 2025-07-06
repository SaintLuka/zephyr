#include <map>
#include <numeric>
#include <iomanip>

#include <zephyr/utils/mpi.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>


namespace zephyr::mesh {

using namespace zephyr::io;
using namespace zephyr::utils;

void EuMesh::set_decomposition(const std::string& type) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    auto domain = bbox();
    m_decomp = ORB::create(domain, type, mpi::size());

    set_decomposition(m_decomp, true);
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

void EuMesh::set_decomposition(Decomposition::Ref decmp, bool update) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    m_decomp = decmp;
    if (update) {
        // вызываю, чтобы инициализировать m_tourists
    	build_aliens();
    	std::cout << "Has no first redistribute\n";
        //redistribute();

        // Вероятно, первый (и единственный) redistribute, почистим память
    	m_locals.shrink_to_fit();
    	m_aliens.shrink_to_fit();
    	m_tourists.shrink_to_fit();
    	m_migrants.clear();
    }
#endif
}

void EuMesh::prebalancing(int n_iters) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) { return; }

    for (int i = 0; i < n_iters; ++i) {
        balancing();
    }
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


void EuMesh::setup_ranks() {
    // Определим новый rank для всех ячеек из locals
    for (index_t ic = 0; ic < m_locals.size(); ++ic) {
        EuCell cell(&m_locals, ic);
        m_locals.rank[ic] = m_decomp->rank(cell);
    }
}

void EuMesh::build_aliens() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    mpi::barrier();
    mpi::for_each([]() {
        std::cout << "BUILD ALIENS " << mpi::rank() << "\n";
    });
    m_tourists.build_aliens(m_locals, m_aliens);
    mpi::for_each([]() {
        std::cout << "END BUILD ALIENS " << mpi::rank() << "\n";
    });
#endif
}

} // namespace zephyr::mesh
