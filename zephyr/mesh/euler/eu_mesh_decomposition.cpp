#include <map>
#include <numeric>

#include <zephyr/utils/mpi.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bfaces.h>
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





void EuMesh::exchange() {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;

	m_tourism.exchange_start(m_locals);
	m_tourism.exchange_end(m_aliens);
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

	/* DEBUG
	if(mpi::master()){
		for (auto& cell: m_locals){
			printf("i: %d, r: %d\n", cell.index, cell.rank);
		}
	}*/

	/* DEBUG
	if(mpi::master()){
		printf("rank %d : %d\n", rank, m_locals.size());
		for(int i=0; i<size; ++i)
			printf("m_%d: %d\n", i, m_i[i]);
	}*/


	m_migration.m = mpi::all_gather_vectors(m_migration.m_i);

	// Переиндексируем локальные ячейки
	for(int i = 0; i < rank; ++i)
		for(int s = 0; s < size; ++s)
			m_migration.m_sum[s] += m_migration.m[size * i + s];

	for (auto& cell: m_locals)
		cell.index = m_migration.m_sum[cell.rank]++;

	/* // DEBUG
	if(mpi::master()){
		for(int i=0; i<m_locals.size(); ++i)
			printf("r:%d,i:%d | ", m_locals[i].rank, m_locals[i].index);
		printf("\n");
	}
	*/
	
	exchange();

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

	/*// DEBUG
	if(mpi::master()){
		for(int i=0; i<size; ++i)
			printf("sum: %d\n", m_i_sum[i]);
		printf("\n");
	}
	//*/

	// Заполняем migrants, сортируем по rank
	auto& migrants = m_migration.m_migrants; 
	// [?!] Просто migrants.resize(m_locals.size()); не получится, т.к. тогда неправильно задан itemsize
	migrants = m_locals;
	for (int i = 0; i < migrants.size(); ++i)
		std::memcpy(migrants[m_migration.m_i_sum[m_locals[i].rank]++].ptr(), m_locals[i].ptr(), migrants.itemsize());
				
	/*// DEBUG
	if(mpi::master()){
		for(int i=0; i<m_locals.size(); ++i)
			printf("r:%d,i:%d | ", m_locals[i].rank, m_locals[i].index);
		printf("\n");
		for(int i=0; i<migrants.size(); ++i)
			printf("r:%d,i:%d | ", migrants[i].rank, migrants[i].index);
		printf("\n");
	}
	//*/

	// Меняем размер m_locals
	int new_size = 0;
	for(int i=0; i<size; ++i)
		new_size += m_migration.m[rank + size * i];
	m_locals.resize(new_size);

	// Считаем вспомогательные для отправки массивы
	for(int i=0; i<size; ++i)
		m_migration.m_send_counts[i] = m_migration.m_i[i] * migrants.itemsize();
	for(int i = 1; i < size; ++i)
		m_migration.m_send_offsets[i] = m_migration.m_send_offsets[i - 1] + m_migration.m_send_counts[i - 1];
	for(int i=0; i<size; ++i)
		m_migration.m_recv_counts[i] = m_migration.m[rank + size * i] * migrants.itemsize();
	for(int i=1; i<size; ++i)
		m_migration.m_recv_offsets[i] = m_migration.m_recv_offsets[i - 1] + m_migration.m_recv_counts[i - 1];
	
	// Отправляем всем процессам соостветствующие
	MPI_Alltoallv(
		migrants.item(0).ptr(), m_migration.m_send_counts.data(), m_migration.m_send_offsets.data(), MPI_BYTE,
		m_locals.item(0).ptr(), m_migration.m_recv_counts.data(), m_migration.m_recv_offsets.data(), MPI_BYTE, 
		mpi::comm()
	);

	/*// DEBUG
	if(rank == 3){
		for(int i=0; i<m_locals.size();++i)
		printf("r:%d, i:%d | ", m_locals[i].rank, m_locals[i].index);
		printf("\n");
	}
	//*/

#endif
}

/// @brief 
/// Заполняет m_tourism
void EuMesh::build_aliens() {
#ifdef ZEPHYR_MPI
	int size = mpi::size();
	int rank = mpi::rank();

	m_tourism.reset();

	// Заполняем m_border_indices
	for (auto cell: *this) {
	    // build alien можно вызвать для не совсем нормальной сетки
	    if (cell.geom().is_undefined()) {
            continue;
	    }
		for (auto face: cell.faces()) {
			//if(mpi::master())
			//	printf("face.adjacent().rank: %d\n", face.adjacent().rank);
			auto& border_indices = m_tourism.m_border_indices[face.adjacent().rank];
			if(face.adjacent().rank != rank && (border_indices.empty() || border_indices.back() != cell.geom().index))
				border_indices.push_back(cell.geom().index);
		}
	}

	// Заполняем m_count_to_send
	for(int r = 0; r < size; ++r)
		m_tourism.m_count_to_send[r] = m_tourism.m_border_indices[r].size();
	// Заполняем m_send_offsets
	for(int r = 1; r < size; ++r){
		m_tourism.m_send_offsets[r] = m_tourism.m_send_offsets[r - 1] + m_tourism.m_count_to_send[r - 1];
	}

	// Отправляем m_count_to_send -> получаем m_count_to_recv
	mpi::all_to_all(m_tourism.m_count_to_send, m_tourism.m_count_to_recv);

	// Заполняем m_recv_offsets
	for(int r = 1; r < size; ++r){
		m_tourism.m_recv_offsets[r] = m_tourism.m_recv_offsets[r - 1] + m_tourism.m_count_to_recv[r - 1];
	}
	// Заполняем m_border
	int border_size = m_tourism.m_send_offsets[size - 1] + m_tourism.m_count_to_send[size - 1];
	m_tourism.m_border = m_locals;
	m_tourism.m_border.resize(border_size);

	int temp_border_it = 0;
	for(auto& border_indices : m_tourism.m_border_indices){
		for(auto cell_index : border_indices){
			// [?]
			memcpy(m_tourism.m_border[temp_border_it].ptr(), m_locals[cell_index].ptr(), m_locals.itemsize());

			++temp_border_it;
		}
	}

	m_aliens.resize(m_tourism.m_recv_offsets[size - 1] + m_tourism.m_count_to_recv[size - 1]);

	for(int i = 0; i < size; ++i){
		m_tourism.m_count_to_send[i] *= m_aliens.itemsize();
		m_tourism.m_count_to_recv[i] *= m_aliens.itemsize();
		m_tourism.m_send_offsets[i] *= m_aliens.itemsize();
		m_tourism.m_recv_offsets[i] *= m_aliens.itemsize();
	}

	// Отправляем 
	m_tourism.exchange_start(m_locals);

	for(auto& cell : m_locals){
		for(auto& face : cell.faces){
			if (face.is_undefined()) 
				continue;
				
			if(face.adjacent.rank == rank)
				face.adjacent.alien = -1;
		}
	}

	// Получам в aliens
	m_tourism.exchange_end(m_aliens);

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