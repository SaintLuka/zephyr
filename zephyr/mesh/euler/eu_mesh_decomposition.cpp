#include <map>
#include <numeric>

#include <zephyr/utils/mpi.h>

#include <zephyr/geom/grid.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/primitives/bfaces.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_cell.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::mesh {

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
        // вызываю, чтобы инициализировать m_tourism
        build_aliens();
        redistribute();
    }
#endif

}

void EuMesh::balancing(double load){
#ifdef ZEPHYR_MPI
	if (mpi::single()) { return; }

	auto ws = mpi::all_gather(load);
	m_decomp->balancing(ws);
#endif
}

void EuMesh::redistribute() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

	migrate();
	build_aliens();
#endif
}

void EuMesh::send(Post post) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    m_tourism.send(m_locals, post);
#endif
}

void EuMesh::recv(Post post) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    m_tourism.recv(m_aliens, post);
#endif
}

void EuMesh::sync(Post post) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

    m_tourism.send(m_locals, post);
    m_tourism.recv(m_aliens, post);
#endif
}

void EuMesh::migrate() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

	int size = mpi::size();
	int rank = mpi::rank();

	m_migration.reset();

	// По некоторому правилу определяется новый rank для всех ячеек из массива locals
	for (auto& cell: m_locals){
		cell.rank = m_decomp->rank(cell);

		// Подсчитываем число ячеек, которые должны быть перемещены с данного процесса на другие
		++m_migration.m_i[cell.rank];
	}

	m_migration.m = mpi::all_gather_vectors(m_migration.m_i);

	// Переиндексируем локальные ячейки
	for(int i = 0; i < rank; ++i)
		for(int s = 0; s < size; ++s)
			m_migration.m_sum[s] += m_migration.m[size * i + s];

	for (auto& cell: m_locals)
		cell.index = m_migration.m_sum[cell.rank]++;

    sync();

    // Переиндексируем грани
	for (auto& cell: m_locals){
		for(auto& face: cell.faces){
			if (face.is_undefined())
				continue;
			
			if(face.adjacent.alien == -1) {
			    // Зачем этот код??
				face.adjacent.rank  = m_locals[face.adjacent.index].rank;
				face.adjacent.index = m_locals[face.adjacent.index].index;
			} else {
				face.adjacent.rank  = m_aliens[face.adjacent.alien].rank;
				face.adjacent.index = m_aliens[face.adjacent.alien].index;
			}
		}
	}

	// Для сортировки в migrants по ранку за О(n). m_i_sum[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.
	for(int i = 1; i < size; ++i)
		m_migration.m_i_sum[i] = m_migration.m_i_sum[i - 1] + m_migration.m_i[i - 1];

	// Заполняем migrants, сортируем по rank
	auto& migrants = m_migration.m_migrants; 
	// [?!] Просто migrants.resize(m_locals.size()); не получится, т.к. тогда неправильно задан itemsize
	migrants = m_locals;
	for (int i = 0; i < migrants.size(); ++i)
		std::memcpy(migrants[m_migration.m_i_sum[m_locals[i].rank]++].ptr(), m_locals[i].ptr(), migrants.itemsize());

	// Меняем размер m_locals
	int new_size = 0;
	for(int i=0; i<size; ++i)
		new_size += m_migration.m[rank + size * i];
	m_locals.resize(new_size);

	// Считаем вспомогательные для отправки массивы
	for(int i = 0; i < size; ++i)
		m_migration.m_send_counts[i] = m_migration.m_i[i];
	for(int i = 1; i < size; ++i)
		m_migration.m_send_offsets[i] = m_migration.m_send_offsets[i - 1] + m_migration.m_send_counts[i - 1];
	for(int i = 0; i < size; ++i)
		m_migration.m_recv_counts[i] = m_migration.m[rank + size * i];
	for(int i = 1; i < size; ++i)
		m_migration.m_recv_offsets[i] = m_migration.m_recv_offsets[i - 1] + m_migration.m_recv_counts[i - 1];
	
	// Отправляем всем процессам соостветствующие
	// printf("calling...\n");
	MPI_Alltoallv(
		migrants.item(0).ptr(), m_migration.m_send_counts.data(), m_migration.m_send_offsets.data(), m_tourism.get_mpi_type(),
		m_locals.item(0).ptr(), m_migration.m_recv_counts.data(), m_migration.m_recv_offsets.data(), m_tourism.get_mpi_type(), 
		mpi::comm()
	);

#endif
}

/// @brief 
/// Заполняет m_tourism
void EuMesh::build_aliens() {
#ifdef ZEPHYR_MPI
	int size = mpi::size();
	int rank = mpi::rank();

	m_tourism.reset();

	m_tourism.build_border(m_locals, m_aliens);

	// Отправляем 
    m_tourism.send(m_locals);

	for(auto& cell : m_locals){
		for(auto& face : cell.faces){
			if (face.is_undefined()) 
				continue;
				
			if(face.adjacent.rank == rank)
				face.adjacent.alien = -1;
		}
	}

	// Получам в aliens
    m_tourism.recv(m_aliens);

	int al_it = 0;
	for(auto& cell : m_aliens){
		for(auto& face : cell.faces){
			if (face.is_undefined()) 
				continue;

			if(face.adjacent.index != -1 && face.adjacent.rank == rank){
				//printf("I: %d\n", face.adjacent.index);
				auto& curr_cell = m_locals[face.adjacent.index];
				for(auto& l_face : curr_cell.faces){
					if(l_face.adjacent.rank == cell.rank && l_face.adjacent.index == cell.index){
						l_face.adjacent.alien = al_it;
						break;
					}
				}
			}
		}
		++al_it;
	}	

	// [?] если m_border_indices[r].size() == m_count_to_send[r], то зачем второе вообще нужно?
#endif
}

} // namespace zephyr::mesh