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
// [?] почему после migrate()?
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
    build_aliens();
}

void EuMesh::exchange() {

}

// [!] хардкодная рандомная функция реранка
void rerank(zephyr::geom::AmrCell& cell){
    int size = mpi::size();
    cell.rank = cell.index % size;
}

void EuMesh::migrate() {
    int size = mpi::size();
    int rank = mpi::rank();

    std::vector<int> m_i(size, 0);
    // По некоторому правилу определяется новый rank для всех ячеек из массива locals
    for (auto& cell: m_locals){
        rerank(cell);
        // Подсчитываем число ячеек, которые должны быть перемещены с данного процесса на другие
        ++m_i[cell.rank];
    }

    // [?] Надеюсь, это правильно
    std::vector<int> m = mpi::all_gather(m_i);

    // Переиндексируем локальные ячейки
    std::vector<int> m_sum(size, 0);
    for(int i = 0; i < rank - 1; ++i)
        for(int s = 0; s < size; ++s)
            m_sum[s] += m[s][i];

    for (auto& cell: m_locals)
        cell.index = m_sum[cell.rank]++;

    // ... //
}

/// @brief 
/// Заполняет m_tourism
void EuMesh::build_aliens() {
    int size = mpi::size();
    int rank = mpi::rank();
    // Выделяем память
    // [!] думаю, нужно перенести этот код, чтобы он выполнялся единожды
    m_tourism.m_border_indices.resize(size);
    m_tourism.m_count_to_send.resize(size);
    m_tourism.m_send_offsets.resize(size, 0);
    m_tourism.m_recv_offsets.resize(size, 0);

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
    for(int r = 0; r < size; ++r)
        m_tourism.m_count_to_send[r] = m_tourism.m_border_indices[r].size();
    // Заполняем m_send_offsets
    for(int r = 1; r < size; ++r)
        m_tourism.m_send_offsets[r] = m_tourism.m_send_offsets[r - 1] + m_tourism.m_count_to_send[r - 1];

    // Отправляем m_count_to_send -> получаем m_count_to_recv
    // [!] хардкодово... (лишний раз создаю вектор)
    int TAG = 0;
    for(int r = 0; r < size; ++r)
        if(r != rank) mpi::send({m_tourism.m_count_to_send[r]}, r, TAG);
    for(int r = 0; r < size; ++r)
        if(r != rank) mpi::recv({m_tourism.m_count_to_recv[r]}, sizeof(int), r, TAG);

    // Заполняем m_recv_offsets
    for(int r = 1; r < size; ++r)
        m_tourism.m_recv_offsets[r] = m_tourism.m_recv_offsets[r - 1] + m_tourism.m_count_to_recv[r - 1];

    // Заполняем m_border
    int border_size = m_tourism.m_send_offsets[size - 1] + m_tourism.m_count_to_send[size - 1];
    m_tourism.m_border.resize(border_size);

    int temp_border_it = 0;
    for(auto& border_indices : m_tourism.m_border_indices)
        for(auto& cell_index : border_indices)
            m_tourism.m_border[temp_border_it++] = m_locals[cell_index];

    // [?] если m_border_indices[r].size() == m_count_to_send[r], то зачем второе вообще нужно?
}

} // namespace zephyr::mesh