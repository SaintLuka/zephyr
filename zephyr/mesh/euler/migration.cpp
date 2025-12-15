#include <zephyr/mesh/euler/migration.h>
#include <zephyr/io/vtu_file.h>

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI
using utils::mpi;

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

void Migration::clear() {
    migrants.clear();
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
    m_cell_router.set_send_count(cell_send_count);
    m_face_router.set_send_count(face_send_count);
    m_node_router.set_send_count(node_send_count);

    // Получение полной матрицы пересылок (all to all)
    m_cell_router.fill_complete();
    m_face_router.fill_complete();
    m_node_router.fill_complete();

    /*
    // Количество разных штук на отправку и получение
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "Cell Router\n";
        m_cell_route.print();
        std::cout << "Face Router\n";
        m_face_route.print();
        std::cout << "Node Router\n";
        m_node_route.print();
    });
     */
}

void Migration::reindexing(Tourism& tourism, AmrCells& locals, AmrCells& aliens) {
    // Отправим новые ранги ячеек
    tourism.prepare<MpiTag::RANK>(locals);
    auto send_rnk_1 = tourism.isend<MpiTag::RANK>();
    auto recv_rnk_1 = tourism.irecv<MpiTag::RANK>();

    // Переиндексируем локальные ячейки

    // Стартовые индексы ячеек, граней и вершин
    std::vector<index_t> cell_index(mpi::size(), 0);
    std::vector<index_t> face_index(mpi::size(), 0);
    std::vector<index_t> node_index(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        for (int i = 0; i < mpi::rank(); ++i) {
            cell_index[r] += m_cell_router(i, r);
            face_index[r] += m_face_router(i, r);
            node_index[r] += m_node_router(i, r);
        }
    }

    // Новые индексы ячеек (получатся после миграции)
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        locals.index[ic] = cell_index[locals.rank[ic]]++;
    }

    tourism.prepare<MpiTag::INDEX>(locals);
    auto send_idx_1 = tourism.isend<MpiTag::INDEX>();
    auto recv_idx_1 = tourism.irecv<MpiTag::INDEX>();

    // Дождемся получения новых rank и index ячеек
    send_rnk_1.wait();
    recv_rnk_1.wait();
    send_idx_1.wait();
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