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
    build_aliens();
}

void EuMesh::exchange() {

}

void EuMesh::migrate() {

}

void EuMesh::build_aliens() {

}

} // namespace zephyr::mesh