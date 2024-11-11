#include <map>
#include <numeric>

#include <zephyr/utils/mpi.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bfaces.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_cell.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <problems/fast.h>

namespace zephyr::mesh {

using namespace zephyr::utils;

void EuMesh::add_decomposition(const decomp::ORB& orb, bool update) {
	m_decomp = std::make_shared<ORB>(orb);
	if (update) {
        redistribute();
	}
}

void EuMesh::redistribute() {
	// Следующий бессмысленный код для демонстрации
	/*
	int r = mpi::rank();

	int s = mpi::size();

	// MPI_COMM_WORLD покороче
	mpi::comm();

	// Указатель на элемент хранилища
	Byte * ptr = m_locals[10].ptr();

	// Размер элемента хранилища в байтах
	int is = m_locals.itemsize();

	for (auto& cell: m_locals) {
		// Важные поля
		cell.rank;
		cell.index;

		for (auto& face: cell.faces) {
			// Часть граней пустые, там резервное место
			if (face.is_undefined()) {
				continue;
			}

			// Противоположный face.is_undefined()
			face.is_actual();

			// Важные поля
			face.adjacent.rank;
			face.adjacent.index;
			face.adjacent.alien;

		}
	}


	// Это вариант обхода как по распределенной сетке,
	// там внутри итераторы сложнее зашиты
	for (auto cell: *this) {
		// К примеру, здесь пропускаются неактуальные грани.
		for (auto face: cell.faces()) {

			// Ещё здесь нельзя просто так получить  геометрические данные,
			// для этого и нужно пользоваться итераторами по m_locals
			cell.geom().rank = 10;

			face.adjacent().rank;
		}
	}
	*/

	migrate();
	build_aliens();
}

void EuMesh::exchange() {
#ifdef ZEPHYR_ENABLE_MPI
	int size = mpi::size();
	int rank = mpi::rank();

	int temp_border_it = 0;
	for(auto& border_indices : m_tourism.m_border_indices){
		for(auto cell_index : border_indices){
			// [?]
			memcpy(m_tourism.m_border[temp_border_it].ptr(), m_locals[cell_index].ptr(), m_locals.itemsize());

			++temp_border_it;
		}
	}

	MPI_Alltoallv(
		m_tourism.m_border.item(0).ptr(), m_tourism.m_count_to_send.data(), m_tourism.m_send_offsets.data(), MPI_BYTE, 
		m_aliens.item(0).ptr(), m_tourism.m_count_to_recv.data(), m_tourism.m_recv_offsets.data(), MPI_BYTE, 
		mpi::comm()
	);
#endif
}


void EuMesh::migrate() {
#ifdef ZEPHYR_ENABLE_MPI
	int size = mpi::size();
	int rank = mpi::rank();

	std::vector<int> m_i(size, 0);
	// По некоторому правилу определяется новый rank для всех ячеек из массива locals
	for (auto& cell: m_locals){
		// m_decomp->rank(cell);
		
		cell.rank = m_decomp->rank(cell);

		// cell.rank = (cell.index / (m_nx*m_nx/2)) * 2 + (cell.index % (m_nx)) / (m_nx/2);
		// Подсчитываем число ячеек, которые должны быть перемещены с данного процесса на другие
		++m_i[cell.rank];
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

	std::vector<int> m = mpi::all_gather_vectors(m_i);

	// Переиндексируем локальные ячейки
	std::vector<int> m_sum(size, 0);
	for(int i = 0; i < rank; ++i)
		for(int s = 0; s < size; ++s)
			m_sum[s] += m[size * i + s];

	for (auto& cell: m_locals)
		cell.index = m_sum[cell.rank]++; 

	/* // DEBUG
	if(mpi::master()){
		for(int i=0; i<m_locals.size(); ++i)
			printf("r:%d,i:%d | ", m_locals[i].rank, m_locals[i].index);
		printf("\n");
	}
	*/

	// Переиндексируем грани
	for (auto& cell: m_locals){
		for(auto& face: cell.faces){
			if(face.adjacent.index >= 0){
				face.adjacent.rank = m_locals[face.adjacent.index].rank;
				face.adjacent.index = m_locals[face.adjacent.index].index;
			}
		}
	}

	// Для сортировки в migrants по ранку за О(n). m_i_sum[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.
	std::vector<int> m_i_sum(size, 0);
	for(int i = 1; i < size; ++i)
		m_i_sum[i] = m_i_sum[i - 1] + m_i[i - 1];

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
		std::memcpy(migrants[m_i_sum[m_locals[i].rank]++].ptr(), m_locals[i].ptr(), migrants.itemsize());

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
		new_size += m[rank + size * i];
	m_locals.resize(new_size);

	// Считаем вспомогательные для отправки массивы
	std::vector<int> send_counts(size, 0);
	for(int i=0; i<size; ++i)
		send_counts[i] = m_i[i] * migrants.itemsize();
	std::vector<int> send_displs(size, 0);
	for(int i = 1; i < size; ++i)
		send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
	std::vector<int> recv_counts(size, 0);
	for(int i=0; i<size; ++i)
		recv_counts[i] = m[rank + size * i] * migrants.itemsize();
	std::vector<int> recv_displs(size, 0);
	for(int i=1; i<size; ++i)
		recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
	
	// Отправляем всем процессам соостветствующие
	MPI_Alltoallv(
		migrants.item(0).ptr(), send_counts.data(), send_displs.data(), MPI_BYTE, 
		m_locals.item(0).ptr(), recv_counts.data(), recv_displs.data(), MPI_BYTE, 
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
#ifdef ZEPHYR_ENABLE_MPI
	int size = mpi::size();
	int rank = mpi::rank();

	// Выделяем память
	// [!] думаю, нужно перенести этот код, чтобы он выполнялся единожды
	m_tourism.m_border_indices.resize(size, std::vector<int>());
	m_tourism.m_count_to_send.resize(size);
	m_tourism.m_count_to_recv.resize(size);
	m_tourism.m_send_offsets.resize(size, 0);
	m_tourism.m_recv_offsets.resize(size, 0);

	// Заполняем m_border_indices
	for (auto cell: *this) {
		for (auto face: cell.faces()) {
			auto& border_indices = m_tourism.m_border_indices[face.adjacent().rank];
			if(face.adjacent().rank != rank && (border_indices.empty() || border_indices.back() != cell.geom().index))
				border_indices.push_back(cell.geom().index);
		}
	}

	// Заполняем m_count_to_send
	for(int r = 0; r < size; ++r)
		m_tourism.m_count_to_send[r] = m_tourism.m_border_indices[r].size();
	// Заполняем m_send_offsets
	for(int r = 1; r < size; ++r)
		m_tourism.m_send_offsets[r] = m_tourism.m_send_offsets[r - 1] + m_tourism.m_count_to_send[r - 1];

	// Отправляем m_count_to_send -> получаем m_count_to_recv
	mpi::all_to_all(m_tourism.m_count_to_send, m_tourism.m_count_to_recv);

	// Заполняем m_recv_offsets
	for(int r = 1; r < size; ++r)
		m_tourism.m_recv_offsets[r] = m_tourism.m_recv_offsets[r - 1] + m_tourism.m_count_to_recv[r - 1];

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

	for(int i=0; i<size; ++i){
		m_tourism.m_count_to_send[i] *= m_aliens.itemsize();
		m_tourism.m_count_to_recv[i] *= m_aliens.itemsize();
		m_tourism.m_send_offsets[i] *= m_aliens.itemsize();
		m_tourism.m_recv_offsets[i] *= m_aliens.itemsize();
	}

	// Отправляем и получам в aliens
	exchange();

	for(auto& cell : m_locals){
		for(auto& face : cell.faces){
			if(face.adjacent.rank == rank)
				face.adjacent.alien = -1;
		}
	}

	int i = 0;
	for(auto& cell : m_aliens){
		for(auto& face : cell.faces){
			if(face.adjacent.rank == rank){
				auto& curr_cell = m_locals[face.adjacent.index];
				for(auto& l_face : curr_cell.faces)
					if(l_face.adjacent.index == cell.index){
						l_face.adjacent.index = -1;
						l_face.adjacent.alien = i++;
						break;
					}
			}
		}
	}

	// [?] если m_border_indices[r].size() == m_count_to_send[r], то зачем второе вообще нужно?
#endif
}

} // namespace zephyr::mesh