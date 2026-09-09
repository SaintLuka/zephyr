#include <zephyr/mesh/euler/migration.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/mesh/euler/eu_mesh.h>

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

void Migration::init_types(const AmrCells& cells) {
    cell_buffer_ = AmrCells(cells.options());
}

void Migration::clear() {
    cell_buffer_.clear();
    node_buffer_.clear();
}

void Migration::shrink_to_fit() {
    cell_buffer_.shrink_to_fit();
    node_buffer_.shrink_to_fit();
}

void Migration::setup_node_ranks(AmrNodes& nodes,
    const AmrCells& locals, const AmrCells& ghosts) {

    for (index_t in = 0; in < nodes.n_nodes(); ++in) {
        int r = std::numeric_limits<int>::max();
        for (index_t inc: nodes.incident.range(in)) {
            if (nodes.incident.is_undefined(inc)) continue;

            int rnk, idx;
            if (nodes.incident.is_local(inc)) {
                rnk = locals.rank[nodes.incident.index[inc]];
                idx = locals.index[nodes.incident.index[inc]];
            }
            else {
                rnk = ghosts.rank[nodes.incident.ghost[inc]];
                idx = ghosts.index[nodes.incident.ghost[inc]];
            }
            r = std::min(r, rnk);

            // Устанавливаем индекс, который будет после миграции
            nodes.incident.rank [inc] = rnk;
            nodes.incident.index[inc] = idx;
        }
        nodes.rank[in] = r;
    }
}

void Migration::fill_cell_routers(const AmrCells& cells) {
    // Сколько элементов каждого ранга пересылается с данного процесса
    // при миграции, включая пересылки на сам процесс (пересылка r -> r)
    // Строка матрицы пересылок
    std::vector<index_t> cell_send_count(mpi::size(), 0);
    std::vector<index_t> face_send_count(mpi::size(), 0);
    std::vector<index_t> vert_send_count(mpi::size(), 0);
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        int new_rank = cells.rank[ic];

        // Число ячеек, граней и вершин, которые нужно переслать с данного процесса на другие
        cell_send_count[new_rank] += 1;
        face_send_count[new_rank] += cells.faces.max_count(ic);
        vert_send_count[new_rank] += cells.verts.max_count(ic);
    }

    // Установить число ячеек на отправку
    cell_router_.set_send_count(cell_send_count);
    face_router_.set_send_count(face_send_count);
    vert_router_.set_send_count(vert_send_count);

    // Получение полной матрицы пересылок (all to all)
    cell_router_.fill_complete();
    face_router_.fill_complete();
    vert_router_.fill_complete();

    /*
    // Количество разных штук на отправку и получение
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "Cell Router\n";
        cell_router_.print();
        std::cout << "Face Router\n";
        face_router_.print();
        std::cout << "Vert Router\n";
        vert_router_.print();
    });
    */
}

void Migration::fill_node_routers(const AmrNodes& nodes, bool unique_nodes) {
    if (!unique_nodes) {
        node_router_.set_zero_complete();
        inct_router_.set_zero_complete();
        return;
    }

    // Сколько элементов каждого ранга пересылается с данного процесса
    // при миграции, включая пересылки на сам процесс (пересылка r -> r)
    // Строка матрицы пересылок
    std::vector<index_t> node_send_count(mpi::size(), 0);
    std::vector<index_t> inct_send_count(mpi::size(), 0);
    for (index_t in = 0; in < nodes.n_nodes(); ++in) {
        int new_rank = nodes.rank[in];

        // Число узлов и инцидентных ячеек, которые нужно переслать с данного процесса на другие
        node_send_count[new_rank] += 1;
        inct_send_count[new_rank] += nodes.incident.max_count(in);
    }

    // Установить число узлов на отправку
    node_router_.set_send_count(node_send_count);
    inct_router_.set_send_count(inct_send_count);

    // Получение полной матрицы пересылок (all to all)
    node_router_.fill_complete();
    inct_router_.fill_complete();

    // Количество разных штук на отправку и получение
    /*
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "Node Router\n";
        node_router_.print();
        std::cout << "Incident Router\n";
        inct_router_.print();
    });
    */
}

void Migration::cells_reindexing(Tourism& tourism, AmrCells& cells) const {
    // Отправим новые ранги ячеек
    tourism.prepare<MpiTag::RANK>(cells);
    auto send_rnk = tourism.isend<MpiTag::RANK>();
    auto recv_rnk = tourism.irecv<MpiTag::RANK>();

    // Переиндексируем локальные ячейки

    // Стартовые индексы ячеек, граней и вершин
    std::vector<index_t> cell_index(mpi::size(), 0);
    std::vector<index_t> face_index(mpi::size(), 0);
    std::vector<index_t> vert_index(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        for (int i = 0; i < mpi::rank(); ++i) {
            cell_index[r] += cell_router_(i, r);
            face_index[r] += face_router_(i, r);
            vert_index[r] += vert_router_(i, r);
        }
    }

    // Новые индексы ячеек (получатся после миграции)
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        cells.index[ic] = cell_index[cells.rank[ic]]++;
    }

    tourism.prepare<MpiTag::INDEX>(cells);
    auto send_idx = tourism.isend<MpiTag::INDEX>();
    auto recv_idx = tourism.irecv<MpiTag::INDEX>();

    // Дождемся получения новых rank и index ячеек
    send_rnk.wait();
    recv_rnk.wait();
    send_idx.wait();
    recv_idx.wait();
}

void Migration::nodes_reindexing(Tourism& tourism, AmrNodes& nodes) const {
    // Отправим новые ранги ячеек
    tourism.prepare<MpiTag::NODE_RANK>(nodes);
    auto send_rnk = tourism.isend<MpiTag::NODE_RANK>();
    auto recv_rnk = tourism.irecv<MpiTag::NODE_RANK>();

    // Переиндексируем локальные узлы

    // Стартовые индексы узлов, инцидентных ячеек
    std::vector<index_t> node_index(mpi::size(), 0);
    std::vector<index_t> inct_index(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        for (int i = 0; i < mpi::rank(); ++i) {
            node_index[r] += node_router_(i, r);
            inct_index[r] += inct_router_(i, r);
        }
    }

    // Новые индексы узлов (получатся после миграции)
    for (index_t in = 0; in < nodes.n_nodes(); ++in) {
        nodes.index[in] = node_index[nodes.rank[in]]++;
    }

    tourism.prepare<MpiTag::NODE_INDEX>(nodes);
    auto send_idx = tourism.isend<MpiTag::NODE_INDEX>();
    auto recv_idx = tourism.irecv<MpiTag::NODE_INDEX>();

    // Дождемся получения новых rank и index узлов
    send_rnk.wait();
    recv_rnk.wait();
    send_idx.wait();
    recv_idx.wait();
}

void Migration::update_face_adjacent(AmrCells& locals, const AmrCells& ghosts) {
    auto& faces = locals.faces;
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        for(auto iface: locals.faces.range(ic)){
            if (faces.is_undefined(iface))
                continue;

            index_t idx = faces.adjacent.index[iface];
            index_t gst = faces.adjacent.ghost[iface];

            if(faces.adjacent.is_local(iface)) {
                faces.adjacent.rank [iface] = locals.rank [idx];
                faces.adjacent.index[iface] = locals.index[idx];
            } else {
                faces.adjacent.rank [iface] = ghosts.rank [gst];
                faces.adjacent.index[iface] = ghosts.index[gst];
            }
            faces.adjacent.basic[iface] = locals.index[ic];
        }
    }
}
void Migration::update_cell_verts(AmrVerts& verts, const AmrNodes& locals, const AmrNodes& ghosts) {
    for (index_t inode = 0; inode < verts.n_verts(); ++inode) {
        index_t idx = verts.index[inode];
        index_t gst = verts.ghost[inode];

        if (gst < 0) { // local node
            verts.rank [inode] = locals.rank [idx];
            verts.index[inode] = locals.index[idx];
        }
        else {
            verts.rank [inode] = ghosts.rank [idx];
            verts.index[inode] = ghosts.index[idx];
        }
    }
}

#endif

}