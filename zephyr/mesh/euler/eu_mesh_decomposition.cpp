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
}

void EuMesh::redistribute() {
    // Следующий бессмысленный код для демонстрации

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


    migrate();
    build_aliens(); // [?] почему после migrate()?
}

void EuMesh::exchange() {

}

void EuMesh::migrate() {

}

void EuMesh::build_aliens() {
    int rank = mpi::rank();

    // Выделяем память
    // [!] думаю, нужно перенести этот код, чтобы он выполнялся единожды
    m_tourism.m_border_indices.resize(rank);
    m_tourism.m_count_to_send.resize(rank);
    m_tourism.m_send_offsets.resize(rank, 0);

    // Заполняем m_border_indices
    for (auto& cell: m_locals) {
        for (auto& face: cell.faces) {
            auto& border_indices = m_tourism.m_border_indices[face.adjacent.rank];
            if(face.adjacent.rank != rank && (border_indices.empty() || border_indices.back() != cell.index)){
                // [!?] вместо push_back лучше заранее зарезервировать память m_border_indices и использовать присвоение
                m_tourism.m_border_indices[face.adjacent.rank].push_back(cell.index);
            }
        }
    }

    // Заполняем m_count_to_send
    for(int r = 0; r < rank; ++r)
        m_tourism.m_count_to_send[r] = m_tourism.m_border_indices[r].size();
    // Заполняем m_send_offsets
    for(int r = 1; r < rank; ++r)
        m_tourism.m_send_offsets[r] = m_tourism.m_send_offsets[r - 1] + m_tourism.m_count_to_send[r - 1];

    // Заполняем m_border
    int border_size = m_tourism.m_send_offsets[rank - 1] + m_tourism.m_count_to_send[rank - 1];
    m_tourism.m_border.resize(border_size);

    int temp_border_it = 0;
    for(auto& border_indices : m_tourism.m_border_indices)
        for(auto& cell_index : border_indices)
            m_tourism.m_border[temp_border_it++] = m_locals[cell_index];

    // [?] если m_border_indices[r].size() == m_count_to_send[r], то зачем второе вообще нужно?
}

} // namespace zephyr::mesh