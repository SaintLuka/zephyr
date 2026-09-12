#include <bitset>
#include <numeric>
#include <map>
#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/raw/router.h>
#include <zephyr/mesh/raw/tourism.h>
#include <zephyr/utils/threads.h>

#ifdef ZEPHYR_MPI

namespace zephyr::mesh {

using utils::mpi;
using utils::threads;
namespace indexing = geom::indexing;

namespace {

// Вместо смещений записать разности, последнее значение -1.
void pack_offsets(std::span<index_t> offsets) {
    z_assert(!offsets.empty(), "pack_offsets: empty offsets array");

    index_t n_values = offsets.size() - 1;
    for (index_t i = 0; i < n_values; ++i) {
        offsets[i] = offsets[i + 1] - offsets[i];
    }
    offsets.back() = -1;
}

void pack_offsets(std::span<index_t> offsets_1, std::span<index_t> offsets_2) {
    z_assert(!offsets_1.empty(), "pack_offsets: empty offsets_1 or offsets_2");
    z_assert(offsets_1.size() == offsets_2.size(), "pack_offsets: offsets_1 and offsets_2 must have the same size");

    index_t n_values = offsets_1.size() - 1;
    for (index_t i = 0; i < n_values; ++i) {
        offsets_1[i] = offsets_1[i + 1] - offsets_1[i];
        offsets_2[i] = offsets_2[i + 1] - offsets_2[i];
    }
    offsets_1.back() = -1;
    offsets_2.back() = -1;
}

// Распаковать смещения, если они записаны функцией pack_offsets (для border-слоев)
void unpack_offsets_ver1(std::span<index_t> offsets) {
    z_assert(!offsets.empty(), "unpack_offsets_ver1: empty offsets array");

    int prev_count = offsets[0];
    offsets[0] = 0;
    index_t n_values = offsets.size() - 1;
    for (index_t i = 1; i <= n_values; ++i) {
        int temp_count = offsets[i];
        offsets[i] = offsets[i - 1] + prev_count;
        prev_count = temp_count;
    }
}

void unpack_offsets_ver1(std::span<index_t> offsets_1, std::span<index_t> offsets_2) {
    z_assert(!offsets_1.empty(), "unpack_offsets_ver1: empty offsets_1 or offsets_2");
    z_assert(offsets_1.size() == offsets_2.size(), "unpack_offsets_ver1: offsets_1 and offsets_2 must have the same size");

    int prev_count_1 = offsets_1[0];
    int prev_count_2 = offsets_2[0];

    offsets_1[0] = 0;
    offsets_2[0] = 0;

    index_t n_values = offsets_1.size() - 1;
    for (index_t i = 1; i <= n_values; ++i) {
        int temp_count_1 = offsets_1[i];
        int temp_count_2 = offsets_2[i];

        offsets_1[i] = offsets_1[i - 1] + prev_count_1;
        offsets_2[i] = offsets_2[i - 1] + prev_count_2;

        prev_count_1 = temp_count_1;
        prev_count_2 = temp_count_2;
    }
}

// Распаковать смещения, если они записаны функцией pack_offsets,
// но смещены на единицу в массиве (для ghost-слоев).
void unpack_offsets_ver2(std::span<index_t> offsets) {
    z_assert(!offsets.empty(), "unpack_offsets_v2: empty offsets array");

    index_t n_values = offsets.size() - 1;
    offsets[0] = 0;
    for (index_t i = 0; i < n_values; ++i) {
        offsets[i + 1] += offsets[i];
    }
}

void unpack_offsets_ver2(std::span<index_t> offsets_1, std::span<index_t> offsets_2) {
    z_assert(!offsets_1.empty(), "unpack_offsets_ver2: empty offsets_1 or offsets_2");
    z_assert(offsets_1.size() == offsets_2.size(), "unpack_offsets_ver2: offsets_1 and offsets_2 must have the same size");

    index_t n_values = offsets_1.size() - 1;
    offsets_1[0] = 0;
    offsets_2[0] = 0;
    for (index_t i = 0; i < n_values; ++i) {
        offsets_1[i + 1] += offsets_1[i];
        offsets_2[i + 1] += offsets_2[i];
    }
}

} // anonymous namespace

void Tourism::shrink_to_fit() {
    //unique_border_indices_.shrink_to_fit();
    border_cells_indices_.shrink_to_fit();
    border_cells_.shrink_to_fit();
    ghost_cells_.shrink_to_fit();

    border_nodes_indices_.shrink_to_fit();
    border_nodes_.shrink_to_fit();
    ghost_nodes_.shrink_to_fit();
}

void Tourism::init_types(const RawCells& cells) {
    border_cells_ = cells.same();
    ghost_cells_ = cells.same();
}

void Tourism::resize_border_cells() {
    index_t n_border_cells = cell_router_.send_buffer_size();
    index_t n_border_faces = face_router_.send_buffer_size();
    index_t n_border_verts = vert_router_.send_buffer_size();

    border_cells_.resize(n_border_cells, n_border_faces, n_border_verts);
}

void Tourism::extend_border_cells() {
    index_t n_border_cells = cell_router_.send_buffer_size();
    index_t n_border_faces = face_router_.send_buffer_size();
    index_t n_border_verts = vert_router_.send_buffer_size();

    if (n_border_cells > border_cells_.n_cells()) {
        border_cells_.resize(n_border_cells, n_border_faces, n_border_verts);
    }
}

void Tourism::resize_border_nodes() {
    index_t n_border_nodes    = node_router_.send_buffer_size();
    index_t n_border_incident = inct_router_.send_buffer_size();

    border_nodes_.resize(n_border_nodes, n_border_incident);
}

void Tourism::extend_border_nodes() {
    index_t n_border_nodes    = node_router_.send_buffer_size();
    index_t n_border_incident = inct_router_.send_buffer_size();

    if (n_border_nodes    > border_nodes_.n_nodes() ||
        n_border_incident > border_nodes_.incident.n_values()) {
        border_nodes_.resize(n_border_nodes, n_border_incident);
    }
}

void Tourism::resize_ghost_cells() {
    int n_ghost_cells = cell_router_.recv_buffer_size();
    int n_ghost_faces = face_router_.recv_buffer_size();
    int n_ghost_verts = vert_router_.recv_buffer_size();

    ghost_cells_.resize(n_ghost_cells, n_ghost_faces, n_ghost_verts);
}

void Tourism::extend_ghost_cells() {
    int n_ghost_cells = cell_router_.recv_buffer_size();
    int n_ghost_faces = face_router_.recv_buffer_size();
    int n_ghost_verts = vert_router_.recv_buffer_size();

    if (n_ghost_cells > ghost_cells_.n_cells()) {
        ghost_cells_.resize(n_ghost_cells, n_ghost_faces, n_ghost_verts);
    }
}

void Tourism::resize_ghost_nodes() {
    index_t n_ghost_nodes    = node_router_.recv_buffer_size();
    index_t n_ghost_incident = inct_router_.recv_buffer_size();

    ghost_nodes_.resize(n_ghost_nodes, n_ghost_incident);
}

void Tourism::extend_ghost_nodes() {
    index_t n_ghost_nodes    = node_router_.recv_buffer_size();
    index_t n_ghost_incident = inct_router_.recv_buffer_size();

    if (n_ghost_nodes    > ghost_nodes_.n_nodes() ||
        n_ghost_incident > ghost_nodes_.incident.n_values()) {
        ghost_nodes_.resize(n_ghost_nodes, n_ghost_incident);
    }
}

void Tourism::fill_nodes_send_count(const RawNodes& nodes) {
    if (nodes.empty()) {
        // Сетка без уникальных узлов
        node_router_.set_zero_send_count();
        inct_router_.set_zero_send_count();
        return;
    }

    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последнего учтенного узла
    std::vector<index_t> last_append(size, -1);

    std::vector<index_t> node_send_count(size, 0);
    std::vector<index_t> inct_send_count(size, 0);

    for (index_t in = 0; in < nodes.n_nodes(); ++in) {
        // Для сетки с неактуальными узлами
        if (nodes.is_undefined(in)) continue;

        for (index_t inc: nodes.incident.range(in)) {
            if (nodes.incident.is_undefined(inc)) continue;

            int neib_rank = nodes.incident.rank[inc];
            if (neib_rank != rank && last_append[neib_rank] != in) {
                last_append[neib_rank] = in;
                node_send_count[neib_rank] += 1;
                inct_send_count[neib_rank] += nodes.incident.max_count(in);
            }
        }
    }

    // Установить число на обмены
    node_router_.set_send_count(node_send_count);
    inct_router_.set_send_count(inct_send_count);
}

void Tourism::fill_border_nodes_indices(const RawNodes& nodes) {
    if (nodes.empty()) {
        border_nodes_indices_.clear();
        return;
    }

    const int rank = mpi::rank();

    // Индекс последнего учтенного узла
    std::vector<index_t> last_append(mpi::size(), -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> node_index = node_router_.send_offset();

    border_nodes_indices_.resize(node_router_.send_buffer_size());

    for (index_t in = 0; in < nodes.n_nodes(); ++in) {
        // Для сетки с неактуальными узлами
        if (nodes.is_undefined(in)) continue;

        for (index_t inc: nodes.incident.range(in)) {
            if (nodes.incident.is_undefined(inc)) {
                continue;
            }

            int neib_rank = nodes.incident.rank[inc];
            if (neib_rank != rank) {
                if (last_append[neib_rank] != in) {
                    last_append[neib_rank] = in;
                    border_nodes_indices_[node_index[neib_rank]++] = in;
                }
            }
        }
    }
}

void Tourism::fill_cells_send_count(const RawCells& cells) {
    if (cells.empty()) {
        cell_router_.set_zero_send_count();
        face_router_.set_zero_send_count();
        vert_router_.set_zero_send_count();
        return;
    }

    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    std::vector<index_t> last_append(size, -1);

    std::vector<index_t> cell_send_count(size, 0);
    std::vector<index_t> face_send_count(size, 0);
    std::vector<index_t> vert_send_count(size, 0);

    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (cells.is_undefined(ic)) continue;

        for (index_t iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = cells.faces.adjacent.rank[iface];
            if (neib_rank != rank && last_append[neib_rank] != ic) {
                last_append[neib_rank] = ic;
                cell_send_count[neib_rank] += 1;
                face_send_count[neib_rank] += cells.faces.max_count(ic);
                vert_send_count[neib_rank] += cells.verts.max_count(ic);
            }
        }
    }

    // Установить число на обмены
    cell_router_.set_send_count(cell_send_count);
    face_router_.set_send_count(face_send_count);
    vert_router_.set_send_count(vert_send_count);
}

void Tourism::fill_border_cells_indices(const RawCells& cells) {
    if (cells.empty()) {
        border_cells_indices_.clear();
        return;
    }

    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    index_t last_unique_append = -1;
    std::vector<index_t> last_append(mpi::size(), -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> cell_index = cell_router_.send_offset();

    border_cells_indices_.resize(cell_router_.send_buffer_size());
    // unique_border_indices_.reserve(m_cell_route.send_buffer_size());

    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (cells.is_undefined(ic)) continue;

        for (index_t iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = cells.faces.adjacent.rank[iface];
            if (neib_rank != rank) {
                if (last_append[neib_rank] != ic) {
                    last_append[neib_rank] = ic;
                    border_cells_indices_[cell_index[neib_rank]++] = ic;
                }
                if (last_unique_append != ic) {
                    last_unique_append = ic;
                    // unique_border_indices_.push_back(ic);
                }
            }
        }
    }
}

void Tourism::fill_cells_send_count(const RawCells& cells, const RawNodes& nodes) {
    if (cells.empty()) {
        cell_router_.set_zero_send_count();
        face_router_.set_zero_send_count();
        vert_router_.set_zero_send_count();
        return;
    }

    z_assert(!nodes.empty(), "fill_cells_send_count: empty nodes");

    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    std::vector<index_t> last_append(size, -1);

    std::vector<index_t> cell_send_count(size, 0);
    std::vector<index_t> face_send_count(size, 0);
    std::vector<index_t> vert_send_count(size, 0);

    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (cells.is_undefined(ic)) continue;

        for (index_t iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = cells.faces.adjacent.rank[iface];
            if (neib_rank != rank && last_append[neib_rank] != ic) {
                last_append[neib_rank] = ic;
                cell_send_count[neib_rank] += 1;
                face_send_count[neib_rank] += cells.faces.max_count(ic);
                vert_send_count[neib_rank] += cells.verts.max_count(ic);
            }
        }

        for (index_t i: cells.verts.range(ic)) {
            index_t gst = cells.verts.ghost[i];
            index_t idx = cells.verts.index[i];

            index_t inode = gst < 0 ? idx : gst;
            const RawNodes& neibs = gst < 0 ? nodes : ghost_nodes_;

            if (inode < 0 || inode >= neibs.n_nodes()) {
                std::cout << "Rank " << mpi::rank() << "\n";
                std::cout << "\t" << ic << ", " << i << ", " << inode << ", " << neibs.n_nodes() << ", " << gst << ", " << idx << "\n";
            }

            z_assert(0 <= inode && inode < neibs.n_nodes(), "assrt 1");
            z_assert(!neibs.rank.empty() && !neibs.index.empty(), "assrt 2");

            for (index_t inc: neibs.incident.range(inode)) {
                if (neibs.incident.is_undefined(inc)) {
                    continue;
                }

                z_assert(0 <= inc && inc < neibs.incident.n_values(), "asrt 3");

                int neib_rank = neibs.incident.rank[inc];
                if (neib_rank != rank && last_append[neib_rank] != ic) {
                    last_append[neib_rank] = ic;
                    cell_send_count[neib_rank] += 1;
                    face_send_count[neib_rank] += cells.faces.max_count(ic);
                    vert_send_count[neib_rank] += cells.verts.max_count(ic);
                }
            }
        }
    }

    // Установить число на обмены
    cell_router_.set_send_count(cell_send_count);
    face_router_.set_send_count(face_send_count);
    vert_router_.set_send_count(vert_send_count);
}

void Tourism::fill_border_cells_indices(const RawCells& cells, const RawNodes& nodes) {
    if (cells.empty()) {
        border_cells_indices_.clear();
        return;
    }

    z_assert(!nodes.empty(), "fill_border_cells_indices: empty nodes");

    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    index_t last_unique_append = -1;
    std::vector<index_t> last_append(mpi::size(), -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> cell_index = cell_router_.send_offset();

    border_cells_indices_.resize(cell_router_.send_buffer_size());
    // unique_border_indices_.reserve(m_cell_route.send_buffer_size());

    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (cells.is_undefined(ic)) continue;

        for (index_t iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = cells.faces.adjacent.rank[iface];
            if (neib_rank != rank) {
                if (last_append[neib_rank] != ic) {
                    last_append[neib_rank] = ic;
                    border_cells_indices_[cell_index[neib_rank]++] = ic;
                }
                if (last_unique_append != ic) {
                    last_unique_append = ic;
                    // unique_border_indices_.push_back(ic);
                }
            }
        }

        for (index_t i: cells.verts.range(ic)) {
            index_t gst = cells.verts.ghost[i];
            index_t idx = cells.verts.index[i];

            index_t inode = gst < 0 ? idx : gst;
            const RawNodes& neibs = gst < 0 ? nodes : ghost_nodes_;

            for (index_t inc: neibs.incident.range(inode)) {
                if (neibs.incident.is_undefined(inc)) {
                    continue;
                }

                int neib_rank = neibs.incident.rank[inc];
                if (neib_rank != rank) {
                    if (last_append[neib_rank] != ic) {
                        last_append[neib_rank] = ic;
                        border_cells_indices_[cell_index[neib_rank]++] = ic;
                    }
                    if (last_unique_append != ic) {
                        last_unique_append = ic;
                        // unique_border_indices_.push_back(ic);
                    }
                }
            }
        }
    }
}

void Tourism::prepare_cells_geometry(const RawCells& cells) {
    z_assert(border_cells_.has_nodes() == cells.has_nodes(), "Has no nodes");
    z_assert(ghost_cells_. has_nodes() == cells.has_nodes(), "Has no nodes");

    index_t face_idx = 0;
    index_t vert_idx = 0;
    for (index_t ic = 0; ic < border_cells_indices_.size(); ++ic) {
        cells.copy_geom(border_cells_indices_[ic], border_cells_, ic, face_idx, vert_idx);

        face_idx += cells.faces.max_count(border_cells_indices_[ic]);
        vert_idx += cells.verts.max_count(border_cells_indices_[ic]);
    }
}

void Tourism::prepare_nodes_geometry(const RawNodes& nodes) {
    index_t inc_offset = 0;
    for (index_t in = 0; in < border_nodes_indices_.size(); ++in) {
        nodes.copy_geom(border_nodes_indices_[in], border_nodes_, in, inc_offset);

        inc_offset += nodes.incident.max_count(border_nodes_indices_[in]);
    }
}

void Tourism::build_border_nodes(const RawNodes& nodes) {
    // Заполнить router.send_count
    fill_nodes_send_count(nodes);

    // Заполнить индексы узлов для отправки
    fill_border_nodes_indices(nodes);

    // Подготовить массив border_nodes
    resize_border_nodes();

    // Скопировать геометрию из nodes в border_nodes
    prepare_nodes_geometry(nodes);
}

void Tourism::build_border_cells(const RawCells& cells, const RawNodes& nodes) {
    if (!cells.has_nodes()) {
        // Нет уникальных узлов - окрестность Неймана

        // Заполнить router.send_count
        fill_cells_send_count(cells);

        // Заполнить индексы ячеек для отправки
        fill_border_cells_indices(cells);
    }
    else {
        // Есть уникальные узлы - окрестность Мура

        // Заполнить router.send_count
        fill_cells_send_count(cells, nodes);

        // Заполнить индексы ячеек для отправки
        fill_border_cells_indices(cells, nodes);
    }

    // Подготовить массив border_cells
    resize_border_cells();

    // Перенести геометрию из cells в border_cells
    prepare_cells_geometry(cells);
}

// Инициализация индекса ghost = -1 для большинства граней
void set_undef_ghosts(RawCells& cells, int rank) {
    threads::parallel_for(
        index_t{0}, cells.n_cells(),
        [&cells, rank](index_t ic) {
            for (index_t iface: cells.faces.range(ic)) {
                if (cells.faces.is_actual(iface) &&
                    cells.faces.adjacent.rank[iface] == rank) {
                    cells.faces.adjacent.ghost[iface] = -1;
                }
            }
        });
}

index_t Tourism::find_ghost_node(int rank, index_t index) const {
    const auto beg = ghost_nodes_.index.begin() + node_router_.recv_offset(rank);
    const auto end = beg + node_router_.recv_count(rank);
    auto it = std::lower_bound(beg, end, index);
    if (it != end && *it == index) {
        return std::distance(ghost_nodes_.index.begin(), it);
    }
    return -1;
}

index_t Tourism::find_ghost_cell(int rank, index_t index) const {
    const auto beg = ghost_cells_.index.begin() + cell_router_.recv_offset(rank);
    const auto end = beg + cell_router_.recv_count(rank);
    auto it = std::lower_bound(beg, end, index);
    if (it != end && *it == index) {
        return std::distance(ghost_cells_.index.begin(), it);
    }
    return -1;

}

// выставить index/ghost в verts
void Tourism::find_connections_verts(RawVerts &verts, int rank) const {
    for (index_t iv = 0; iv < verts.n_verts(); ++iv) {
        index_t rnk = verts.rank[iv];
        if (rnk < 0 || rnk >= mpi::size()) {
            throw std::runtime_error("Bad rank");
        }
        if (rnk == rank) {
            verts.ghost[iv] = -1;
        }
        else {
            index_t idx = verts.index[iv];
            index_t gst = find_ghost_node(rnk, idx);
            verts.ghost[iv] = gst;
            if (gst >= 0) {
                z_assert(gst < ghost_nodes_.n_nodes(), "bad ghost");

                geom::Vector3d v1 = verts.coord[iv];
                geom::Vector3d v2 = ghost_nodes_.coord[gst];

                if ((v2 - v1).norm() > 1.0e-15) {
                    std::cout << "Bad rank " << mpi::rank() << ", " << rnk << ", " << iv << ", " << idx << ", " << gst << "\n";
                    std::cout << "  v1: " << v1.transpose() << "; v2: " << v2.transpose() << "\n";

                    for (int i = 0; i < ghost_nodes_.n_nodes(); ++i) {
                        std::cout << "\t" << ghost_nodes_.rank[i] << "; " << ghost_nodes_.index[i] << "; " << ghost_nodes_.coord[i].transpose() << "\n";
                    }
                }
                z_assert((v2 - v1).norm() < 1.0e-15, "diff vert #2");
            }
        }
    }
}

// Обходим ячейки в ghost и ищем связи
void Tourism::find_connections_neumann(RawCells& cells, int rank) const {
    // Инициализация индекса ghost = -1 для большинства граней
    set_undef_ghosts(cells, mpi::rank());

    for (index_t ic = 0; ic < ghost_cells_.n_cells(); ++ic) {
        for (index_t iface: ghost_cells_.faces.range(ic)) {
            if (ghost_cells_.faces.is_undefined(iface)) {
                continue;
            }

            if (ghost_cells_.faces.adjacent.index[iface] >= 0 &&
                ghost_cells_.faces.adjacent.rank[iface] == rank) {
                // Индекс соседа
                index_t jc = ghost_cells_.faces.adjacent.index[iface];

                for (index_t l_face: cells.faces.range(jc)) {
                    if (cells.faces.adjacent.rank [l_face] == ghost_cells_.rank [ic] &&
                        cells.faces.adjacent.index[l_face] == ghost_cells_.index[ic]) {

                        cells.faces.adjacent.ghost[l_face] = ic;
                        break;
                    }
                }
            }
        }
    }
}

// Обходим ячейки в ghost и ищем связи
void Tourism::find_connections_moore(RawCells& cells, RawNodes &nodes, int rank) {
    auto mark_incident = [this, rank](RawNodes& nodes) {
        for (index_t inc = 0; inc < nodes.incident.n_values(); ++inc) {
            if (nodes.incident.is_undefined(inc)) continue;

            int rnk = nodes.incident.rank[inc];
            if (rnk == rank) {
                nodes.incident.ghost[inc] = -1;
            }
            else {
                index_t idx = nodes.incident.index[inc];
                nodes.incident.ghost[inc] = find_ghost_cell(rnk, idx);
            }
        }
    };

    mark_incident(nodes);
    mark_incident(ghost_nodes_);

    auto mark_adjacent = [this, rank](RawCells& cells) {
        for (index_t iface = 0; iface < cells.n_faces(); ++iface) {
            if (cells.faces.is_undefined(iface)) continue;

            int rnk = cells.faces.adjacent.rank[iface];
            if (rnk == rank) {
                cells.faces.adjacent.ghost[iface] = -1;
            }
            else {
                index_t idx = cells.faces.adjacent.index[iface];
                cells.faces.adjacent.ghost[iface] = find_ghost_cell(rnk, idx);
            }
        }
    };

    mark_adjacent(cells);
    mark_adjacent(ghost_cells_);
}

void Tourism::build_ghost_cells(const RawCells& cells, const RawNodes& nodes) {
    // Построить border-слой
    build_border_cells(cells, nodes);

    // Заполнить recv массивы
    cell_router_.fill_partial();
    face_router_.fill_partial();
    vert_router_.fill_partial();

    // Расширить массив ghosts для получения геометрии
    resize_ghost_cells();

    // Отправить и получить геометрию
    sync_cells_geometry();
}

void Tourism::build_ghost_nodes(const RawNodes &nodes) {
    // Построить border-слой узлов
    build_border_nodes(nodes);

    // Заполнить recv массивы
    node_router_.fill_partial();
    inct_router_.fill_partial();

    // Расширить массив ghosts для получения геометрии
    resize_ghost_nodes();

    // Отправить и получить геометрию
    sync_nodes_geometry();
}

void Tourism::update(RawCells& cells, RawNodes& nodes) {
    if (cells.has_nodes()) {
        build_ghost_nodes(nodes);
        find_connections_verts(cells.verts, mpi::rank());
    }
    build_ghost_cells(cells, nodes);

    if (!cells.has_nodes()) {
        find_connections_neumann(cells, mpi::rank());
    }
    else {
        find_connections_moore(cells, nodes, mpi::rank());
    }
}

// Определить дочерние ячейки, которые прилегают к граням с рангом rank.
// То есть дочерние ячейки, которые окажутся в border-блоке ранга rank.
// @param faces Список граней ячеек из border-слоя
// @param ic Индекс родительской ячейки
// @param rank Ранг процесса, к которому ищется прилегание.
// @return bitset<8> - true/false, прилегает дочерняя ячейка или нет.
template<int dim>
std::bitset<8> border_children(const RawFaces& faces, index_t ic, int rank) {
    std::bitset<8> children; children.reset();
    index_t face_beg = Side<dim>::n_subfaces() * ic;
    for (Side<dim> side: Side<dim>::items()) {
        if (faces.is_undefined(face_beg + side[1])) {
            // Simple Face
            index_t iface = face_beg + side;
            if (faces.adjacent.rank[iface] == rank) {
                for (int i: indexing::children(side)) {
                    children[i] = true;
                }
            }
        }
        else { // Complex Face
            for (auto subface: side.subfaces()) {
                index_t iface = face_beg + subface;
                if (faces.adjacent.rank[iface] == rank) {
                    children[indexing::child(subface)] = true;
                }
            }
        }
    }
    return children;
}

template<int dim>
std::vector<index_t> Tourism::setup_border_next() {
    // Уникальные ячейки на огрубление, ключ (b_idx, level, z_idx)
    std::map<std::tuple<index_t, int, index_t>, index_t> coarse_cells;

    std::vector<index_t> n_border_cells(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == mpi::rank()) continue;

        // Проходим по border-блоку и выставляем индексы
        // Исправим NEXT и INDEX только внутри border блоков!
        index_t next_index = 0; // Локальный новый индекс border-ячейки

        coarse_cells.clear();
        // i - индекс в border_cells_, border_cells_indices_
        for (index_t i: cell_router_.send_indices(r)) {
            if (border_cells_.flag[i] == 0) {
                border_cells_.next[i] = next_index;
                next_index += 1;
            }
            else if (border_cells_.flag[i] == 1) {
                // bitset<8> для дочерних ячеек
                auto children = dim == 2 ? border_children<2>(border_cells_.faces, i, r) :
                                           border_children<3>(border_cells_.faces, i, r);
                // Кодируем список дочерних ячеек
                border_cells_.next[i] = amr::pack_children(next_index, children);
                next_index += static_cast<index_t>(children.count());
            }
            else {
                // Самый неприятный случай, border-ячейка огрубляется

                // Полный индекс родительской ячейки
                std::tuple<index_t, index_t, index_t> parent = {
                    border_cells_.b_idx[i],
                    border_cells_.level[i],
                    border_cells_.z_idx[i] / indexing::CpC(dim),
                };

                auto parent_it = coarse_cells.find(parent);
                if (parent_it != coarse_cells.end()) {
                    border_cells_.next[i] = parent_it->second;
                }
                else {
                    coarse_cells[parent] = next_index;
                    border_cells_.next[i] = next_index;
                    next_index += 1;
                }
            }
        }
        n_border_cells[r] = next_index;
    }
    return n_border_cells;
}

template std::vector<index_t> Tourism::setup_border_next<2>();
template std::vector<index_t> Tourism::setup_border_next<3>();

template<int dim>
void Tourism::update_border_indices(const std::vector<index_t>& locals_next) {
    std::vector<index_t> prev_border_indices = border_cells_indices_;
    border_cells_indices_.resize(cell_router_.send_buffer_size());

    index_t last_border_next = 0;
    for (index_t i = 0; i < prev_border_indices.size(); ++i) {
        z_assert(i < border_cells_.flag.size(), "out of range #1521");
        z_assert(i < border_cells_.next.size(), "out of range #1522");
        z_assert(i < prev_border_indices.size(), "out of range #1523");

        if (border_cells_.flag[i] == 0) {
            index_t border_next = border_cells_.next[i];

            z_assert(border_next < border_cells_indices_.size(), "out of range #1524");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1525");

            border_cells_indices_[border_next] = locals_next[prev_border_indices[i]];
            last_border_next = std::max(last_border_next, border_next);
        }
        else if (border_cells_.flag[i] < 0) {
            index_t border_next = border_cells_.next[i];

            z_assert(border_next < border_cells_indices_.size(), "out of range #1526");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1527");

            index_t parent_index = locals_next[prev_border_indices[i]];

            z_assert(parent_index < locals_next.size(), "out of range #1528");

            border_cells_indices_[border_next] = locals_next[parent_index];
            last_border_next = std::max(last_border_next, border_next);
        }
        else {
            auto [border_next, children] = amr::unpack_children(border_cells_.next[i]);
            index_t main_child = locals_next[prev_border_indices[i]];
            for (int c = 0; c < indexing::CpC(dim); ++c) {
                if (children[c]) {
                    if (border_next >= border_cells_indices_.size()) {
                        std::cout << border_cells_.next[i] << "; " << border_next << "; " << children << "; " << border_cells_indices_.size() << "\n";
                    }
                    z_assert(border_next < border_cells_indices_.size(), "out of range #1529");
                    z_assert(main_child + c < locals_next.size(), "out of range #1530");
                    border_cells_indices_[border_next] = locals_next[main_child + c];
                    ++border_next;
                }
            }
            last_border_next = std::max(last_border_next, border_next);
        }
    }
}

template void Tourism::update_border_indices<2>(const std::vector<index_t>& locals_next);
template void Tourism::update_border_indices<3>(const std::vector<index_t>& locals_next);

#define WRITE_DBG 0

template<int dim>
void Tourism::setup_positions(const std::vector<index_t>& locals_next) {
    int rank = mpi::rank();

    // Размеры border-блоков / новое число ячеек на отправку
    auto n_block_cells = setup_border_next<dim>();

    // Отправим значения NEXT, последнее использование старого роутера
    auto send_next = isend<MpiTag::NEXT>();
    auto recv_next = irecv<MpiTag::NEXT>();

    // ========================================================================
    //              Посчитаем смещения для новых border и ghosts
    // ========================================================================

    std::vector<index_t> n_block_faces(mpi::size(), 0);
    std::vector<index_t> n_block_verts(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        n_block_faces[r] = (dim == 2 ? 8 : 24) * n_block_cells[r];
        n_block_verts[r] = (dim == 2 ? 9 : 27) * n_block_cells[r];
    }

    // Получим значения NEXT, далее можем менять роутеры
    send_next.wait();
    recv_next.wait();

#if WRITE_DBG
    static size_t pvd_counter = 0;
    static io::Variables vars = {"flag", "next", "rank", "level", "index", "b_idx", "z_idx"};
    static io::PvdFile bef_border("sp_border_bef", "debug");
    static io::PvdFile aft_border("sp_border_aft", "debug");
    static io::PvdFile bef_ghosts("sp_ghosts_bef", "debug");
    static io::PvdFile aft_ghosts("sp_ghosts_aft", "debug");

    if (pvd_counter == 0) {
        bef_border.variables = vars;
        bef_ghosts.variables = vars;
        aft_border.variables = vars;
        aft_ghosts.variables = vars;
    }

    bef_border.save(border_cells_, pvd_counter);
    bef_ghosts.save(ghost_cells_, pvd_counter);
#endif

    Router prev_router = cell_router_;

    // Установить число на отправку
    cell_router_.set_send_count(n_block_cells);
    face_router_.set_send_count(n_block_faces);
    vert_router_.set_send_count(n_block_verts);

    // Заполнить recv массивы
    cell_router_.fill_partial();
    face_router_.fill_partial();
    vert_router_.fill_partial();

    // ========================================================================
    //          Сделаем глобальную индексацию next в border и ghosts
    // ========================================================================

    // Добавляем смещения, теперь индексы NEXT в border идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // border-ячейки указана следующая позиция внутри нового border-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.send_indices(r)) {
            border_cells_.next[i] += cell_router_.send_offset(r);
        }
    }

    // Добавляем смещения, теперь индексы NEXT в ghost идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // ghost-ячейки указана следующая позиция внутри нового ghost-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.recv_indices(r)) {
            ghost_cells_.next[i] += cell_router_.recv_offset(r);
        }
    }

#if WRITE_DBG
    aft_border.save(border_cells_, pvd_counter);
    aft_ghosts.save(ghost_cells_, pvd_counter);
    ++pvd_counter;
#endif

    // ========================================================================
    //           Подготовим новые массивы, без актуальных данных
    // ========================================================================

    // Подготовить border массив
    extend_border_cells();

    // Подготовить ghosts массив
    extend_ghost_cells();

    // Выставить корректные индексы в border_cells_indices_
    update_border_indices<dim>(locals_next);
}

template <>
void Tourism::setup_positions<0>(const std::vector<index_t>&) {
    // Вызывается для пустых

    // Установить число на отправку
    std::vector<index_t> send_count(mpi::size(), 0);
    cell_router_.set_send_count(send_count);
    face_router_.set_send_count(send_count);
    vert_router_.set_send_count(send_count);

    // Заполнить recv массивы
    cell_router_.fill_partial();
    face_router_.fill_partial();
    vert_router_.fill_partial();
}

template void Tourism::setup_positions<2>(const std::vector<index_t>&);
template void Tourism::setup_positions<3>(const std::vector<index_t>&);


template<int dim>
void set_amr_indices(std::vector<index_t>& faces_beg, std::vector<index_t>& verts_beg) {
    z_assert(faces_beg.size() == verts_beg.size(), "restore amr sizes mismatch");

    constexpr int n_faces = Side<dim>::n_subfaces();
    constexpr int n_verts = dim == 2 ? 9 : 27;

    threads::parallel_for(
        index_t{0}, index_t(faces_beg.size()),
        [&faces_beg, &verts_beg](index_t ic) {
            faces_beg[ic] = n_faces * ic;
            verts_beg[ic] = n_verts * ic;
        });
}

template <int dim>
void set_amr_incident(std::vector<index_t>& offsets) {
    constexpr int n_inc = RawIncident::max_incident_amr(dim);
    threads::parallel_for(
        index_t{0}, index_t(offsets.size()),
        [&offsets](index_t ic) {
            offsets[ic] = n_inc * ic;
        });
}

void Tourism::send_geometry(const RawCells& cells) {
    prepare_cells_geometry(cells);
    sync_cells_geometry();
}

void Tourism::restore_indices(RawCells& cells) const {
    for (index_t ic: border_cells_indices_) {
        for (index_t iface: cells.faces.range(ic)) {
            index_t ghost_index = cells.faces.adjacent.ghost[iface];
            if (ghost_index >= 0) {
                cells.faces.adjacent.index[iface] = ghost_cells_.index[ghost_index];
            }
        }
    }
}

void Tourism::sync_nodes_geometry() {
    bool amr = border_cells_.adaptive();

    if (!amr) {
        // Оптимизируем пересылку инцидентных ячеек
        pack_offsets(border_nodes_.incident.offsets);
    }

    // ============================= ISEND ====================================

    // Отправить данные узлов
    RequestsList nodes_send; nodes_send.reserve(5);
    nodes_send += node_router_.isend(border_nodes_.rank,  MpiTag::NODE_RANK);
    nodes_send += node_router_.isend(border_nodes_.next,  MpiTag::NODE_NEXT);
    nodes_send += node_router_.isend(border_nodes_.index, MpiTag::NODE_INDEX);
    nodes_send += node_router_.isend(border_nodes_.coord, MpiTag::NODE_COORD);
    if (!amr) {
        nodes_send += node_router_.isend(border_nodes_.incident.offsets, MpiTag::INCT_BEG);
    }

    // Отправить списки инцидентных ячеек
    RequestsList incident_send; incident_send.reserve(4);
    incident_send += inct_router_.isend(border_nodes_.incident.role,  MpiTag::INCT_ROLE);
    incident_send += inct_router_.isend(border_nodes_.incident.rank,  MpiTag::INCT_RANK);
    incident_send += inct_router_.isend(border_nodes_.incident.index, MpiTag::INCT_INDEX);
    incident_send += inct_router_.isend(border_nodes_.incident.ghost, MpiTag::INCT_GHOST);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    RequestsList nodes_recv; nodes_recv.reserve(5);
    nodes_recv += node_router_.irecv(ghost_nodes_.rank,  MpiTag::NODE_RANK);
    nodes_recv += node_router_.irecv(ghost_nodes_.next,  MpiTag::NODE_NEXT);
    nodes_recv += node_router_.irecv(ghost_nodes_.index, MpiTag::NODE_INDEX);
    nodes_recv += node_router_.irecv(ghost_nodes_.coord, MpiTag::NODE_COORD);

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    if (!amr) {
        nodes_recv += node_router_.irecv(ghost_nodes_.incident.offsets.data() + 1, MpiTag::INCT_BEG);
    }

    // Получить данные граней
    RequestsList incident_recv; incident_recv.reserve(4);
    incident_recv += inct_router_.irecv(ghost_nodes_.incident.role,  MpiTag::INCT_ROLE);
    incident_recv += inct_router_.irecv(ghost_nodes_.incident.rank,  MpiTag::INCT_RANK);
    incident_recv += inct_router_.irecv(ghost_nodes_.incident.index, MpiTag::INCT_INDEX);
    incident_recv += inct_router_.irecv(ghost_nodes_.incident.ghost, MpiTag::INCT_GHOST);

    // =========================== WAIT ISEND =================================

    nodes_send.wait();     // Завершить отправку узлов
    incident_send.wait();  // Завершить отправку инцидентных ячеек

    // =========================== WAIT IRECV =================================

    nodes_recv.wait();     // Завершить получение узлов
    incident_recv.wait();  // Завершить получение инцидентных ячеек

    if (!amr) {
        // Восстановить смещения в списке инцидентных ячеек
        unpack_offsets_ver2(ghost_nodes_.incident.offsets);

        // Поддерживать смещения в border массиве не обязательно
        // unpack_offsets_ver1(border_nodes_.incident.offsets);
    }
    else {
        if (border_cells_.dim() == 2) {
            set_amr_incident<2>(ghost_nodes_.incident.offsets);
        }
        else {
            set_amr_incident<3>(ghost_nodes_.incident.offsets);
        }
    }
}

void Tourism::sync_cells_geometry() {
    bool amr = border_cells_.adaptive();
    bool axial = border_cells_.axial();

    if (!amr) {
        // Оптимизируем пересылку индексов граней/вершин
        pack_offsets(
            border_cells_.faces.offsets,
            border_cells_.verts.offsets);
    }

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    RequestsList cells_send; cells_send.reserve(16);
    cells_send += cell_router_.isend(border_cells_.rank, MpiTag::RANK);
    cells_send += cell_router_.isend(border_cells_.next, MpiTag::NEXT);
    cells_send += cell_router_.isend(border_cells_.index, MpiTag::INDEX);
    cells_send += cell_router_.isend(border_cells_.flag, MpiTag::FLAG);
    cells_send += cell_router_.isend(border_cells_.level, MpiTag::LEVEL);
    cells_send += cell_router_.isend(border_cells_.b_idx, MpiTag::B_IDX);
    cells_send += cell_router_.isend(border_cells_.z_idx, MpiTag::Z_IDX);
    cells_send += cell_router_.isend(border_cells_.center, MpiTag::CENTER);
    cells_send += cell_router_.isend(border_cells_.volume, MpiTag::VOLUME);
    if (axial) {
        cells_send += cell_router_.isend(border_cells_.volume_alt, MpiTag::VOLUME_ALT);
    }
    if (!amr) {
        cells_send += cell_router_.isend(border_cells_.faces.offsets, MpiTag::FACE_BEG);
        cells_send += cell_router_.isend(border_cells_.verts.offsets, MpiTag::VERT_BEG);
    }

    // Отправить данные граней
    RequestsList faces_send; faces_send.reserve(16);
    faces_send += face_router_.isend(border_cells_.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_send += face_router_.isend(border_cells_.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_send += face_router_.isend(border_cells_.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_send += face_router_.isend(border_cells_.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_send += face_router_.isend(border_cells_.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_send += face_router_.isend(border_cells_.faces.boundary, MpiTag::BOUNDARY);
    faces_send += face_router_.isend(border_cells_.faces.normal, MpiTag::NORMAL);
    faces_send += face_router_.isend(border_cells_.faces.center, MpiTag::FACE_CENTER);
    faces_send += face_router_.isend(border_cells_.faces.area, MpiTag::AREA);
    if (axial) {
        faces_send += face_router_.isend(border_cells_.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_send += face_router_.isend(border_cells_.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    RequestsList verts_send; verts_send.reserve(3);
    verts_send += vert_router_.isend(border_cells_.verts.coord, MpiTag::VERT_COORD);
    if (border_cells_.verts.has_nodes()) {
        verts_send += vert_router_.isend(border_cells_.verts.rank,  MpiTag::VERT_RANK);
        verts_send += vert_router_.isend(border_cells_.verts.index, MpiTag::VERT_INDEX);
        verts_send += vert_router_.isend(border_cells_.verts.ghost, MpiTag::VERT_GHOST);
    }

    // ============================= IRECV ====================================

    // Получить данные ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += cell_router_.irecv(ghost_cells_.rank, MpiTag::RANK);
    cells_recv += cell_router_.irecv(ghost_cells_.next, MpiTag::NEXT);
    cells_recv += cell_router_.irecv(ghost_cells_.index, MpiTag::INDEX);
    cells_recv += cell_router_.irecv(ghost_cells_.flag, MpiTag::FLAG);
    cells_recv += cell_router_.irecv(ghost_cells_.level, MpiTag::LEVEL);
    cells_recv += cell_router_.irecv(ghost_cells_.b_idx, MpiTag::B_IDX);
    cells_recv += cell_router_.irecv(ghost_cells_.z_idx, MpiTag::Z_IDX);
    cells_recv += cell_router_.irecv(ghost_cells_.center, MpiTag::CENTER);
    cells_recv += cell_router_.irecv(ghost_cells_.volume, MpiTag::VOLUME);
    if (axial) {
        cells_recv += cell_router_.irecv(ghost_cells_.volume_alt, MpiTag::VOLUME_ALT);
    }

    if (!amr) {
        // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
        cells_recv += cell_router_.irecv(ghost_cells_.faces.offsets.data() + 1, MpiTag::FACE_BEG);
        cells_recv += cell_router_.irecv(ghost_cells_.verts.offsets.data() + 1, MpiTag::VERT_BEG);
    }

    // Получить данные граней
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += face_router_.irecv(ghost_cells_.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += face_router_.irecv(ghost_cells_.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += face_router_.irecv(ghost_cells_.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_recv += face_router_.irecv(ghost_cells_.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += face_router_.irecv(ghost_cells_.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_recv += face_router_.irecv(ghost_cells_.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += face_router_.irecv(ghost_cells_.faces.normal, MpiTag::NORMAL);
    faces_recv += face_router_.irecv(ghost_cells_.faces.center, MpiTag::FACE_CENTER);
    faces_recv += face_router_.irecv(ghost_cells_.faces.area, MpiTag::AREA);
    if (axial) {
        faces_recv += face_router_.irecv(ghost_cells_.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_recv += face_router_.irecv(ghost_cells_.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    RequestsList verts_recv; verts_recv.reserve(3);
    verts_recv += vert_router_.irecv(ghost_cells_.verts.coord, MpiTag::VERT_COORD);
    if (ghost_cells_.verts.has_nodes()) {
        verts_recv += vert_router_.irecv(ghost_cells_.verts.rank,  MpiTag::VERT_RANK);
        verts_recv += vert_router_.irecv(ghost_cells_.verts.index, MpiTag::VERT_INDEX);
        verts_recv += vert_router_.irecv(ghost_cells_.verts.ghost, MpiTag::VERT_GHOST);
    }

    // =========================== WAIT ISEND =================================

    cells_send.wait();  // Завершить отправку ячеек
    faces_send.wait();  // Завершить отправку граней
    verts_send.wait();  // Завершить отправку вершин

    // =========================== WAIT IRECV =================================

    cells_recv.wait();  // Завершить получение ячеек
    faces_recv.wait();  // Завершить получение граней
    verts_recv.wait();  // Завершить получение вершин

    if (!amr) {
        // Восстановить индексацию граней
        unpack_offsets_ver2(ghost_cells_.faces.offsets, ghost_cells_.verts.offsets);

        // Поддерживать индексацию граней в border массиве не обязательно
        // unpack_offsets_ver1(border_cells_.faces.offsets, ghost_cells.verts.offsets);
    }
    else {
        if (border_cells_.dim() == 2) {
            set_amr_indices<2>(ghost_cells_.faces.offsets, ghost_cells_.verts.offsets);
        }
        else {
            set_amr_indices<3>(ghost_cells_.faces.offsets, ghost_cells_.verts.offsets);
        }
    }
}

} // namespace zephyr::mesh

#endif