#include <zephyr/mesh/euler/migration.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>


namespace zephyr::mesh {

using utils::mpi;

#ifdef ZEPHYR_MPI

inline std::ostream &operator<<(std::ostream &os, const std::vector<index_t> &arr) {
    os << "[";
    for (size_t i = 0; i < arr.size() - 1; ++i) {
        os << arr[i] << ", ";
    }
    if (!arr.empty()) {
        os << arr.back() << "]";
    }
    return os;
}

void Migration::init_types(const AmrCells& locals) {
    migrants = locals.same();
}

void Migration::clear() {
    migrants.resize(0, 0, 0);
}

void Migration::shrink_to_fit() {
    migrants.shrink_to_fit();
}

void Migration::fill_router(AmrCells& locals) {
    // Сколько элементов каждого ранга пересылается с данного процесса
    // при миграции, включая пересылки на сам процесс (пересылка r -> r)
    // Строка матрицы пересылок
    std::vector<index_t> cell_send_count(mpi::size(), 0);
    std::vector<index_t> face_send_count(mpi::size(), 0);
    std::vector<index_t> node_send_count(mpi::size(), 0);
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        int new_rank = locals.rank[ic];

        // Число ячеек, граней и вершин, которые нужно переслать с данного процесса на другие
        cell_send_count[new_rank] += 1;
        face_send_count[new_rank] += locals.max_face_count(ic);
        node_send_count[new_rank] += locals.max_node_count(ic);
    }

    // Установить число ячеек на отправку
    cell_route.set_send_count(cell_send_count);
    face_route.set_send_count(face_send_count);
    node_route.set_send_count(node_send_count);

    // Получение полной матрицы пересылок (all to all)
    cell_route.fill_complete();
    face_route.fill_complete();
    node_route.fill_complete();

    /*
    // Количество разных штук на отправку и получение
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "Cell Router\n";
        cell_route.print();
        std::cout << "Face Router\n";
        face_route.print();
        std::cout << "Node Router\n";
        node_route.print();
    });
     */
}

void Migration::reindexing(Tourism& tourists, AmrCells& locals, AmrCells& aliens) {
    // Отправим новые ранги ячеек
    tourists.prepare<MpiTag::RANK>(locals);
    auto send_rnk_1 = tourists.isend<MpiTag::RANK>();
    auto recv_rnk_1 = tourists.irecv<MpiTag::RANK>(aliens);

    // Переиндексируем локальные ячейки

    // Стартовые индексы ячеек, граней и вершин
    std::vector<index_t> cell_index(mpi::size(), 0);
    std::vector<index_t> face_index(mpi::size(), 0);
    std::vector<index_t> node_index(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        for (int i = 0; i < mpi::rank(); ++i) {
            cell_index[r] += cell_route(i, r);
            face_index[r] += face_route(i, r);
            node_index[r] += node_route(i, r);
        }
    }

    // Новые индексы ячеек (получатся после миграции)
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        locals.index[ic] = cell_index[locals.rank[ic]]++;
    }

    // Отправим новые индексы ячеек на alien слой
    tourists.prepare<MpiTag::INDEX>(locals);
    auto send_idx_1 = tourists.isend<MpiTag::INDEX>();
    auto recv_idx_1 = tourists.irecv<MpiTag::INDEX>(aliens);

    // Дождемся получения новых rank и index ячеек
    send_rnk_1.wait();
    send_idx_1.wait();
    recv_rnk_1.wait();
    recv_idx_1.wait();

    // Переиндексируем грани
    auto& faces = locals.faces;
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        for(auto iface: locals.faces_range(ic)){
            if (faces.is_undefined(iface))
                continue;

            if(faces.adjacent.is_local(iface)) {
                faces.adjacent.rank [iface] = locals.rank [faces.adjacent.index[iface]];
                faces.adjacent.index[iface] = locals.index[faces.adjacent.index[iface]];
            } else {
                faces.adjacent.rank [iface] = aliens.rank [faces.adjacent.alien[iface]];
                faces.adjacent.index[iface] = aliens.index[faces.adjacent.alien[iface]];
            }
            faces.adjacent.basic[iface] = locals.index[ic];
        }
    }
}

#endif

}