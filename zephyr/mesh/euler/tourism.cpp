#include <bitset>
#include <numeric>
#include <map>
#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/euler/router.h>
#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/utils/threads.h>

#ifdef ZEPHYR_MPI

namespace zephyr::mesh {

using utils::mpi;
using utils::threads;
namespace indexing = geom::indexing;

void Tourism::shrink_to_fit() {
    //m_unique_border_indices.shrink_to_fit();
    m_border_indices.shrink_to_fit();
    m_border.shrink_to_fit();
    m_aliens.shrink_to_fit();
}

void Tourism::init_types(const AmrCells& locals) {
    m_border = locals.same();
    m_aliens = locals.same();
}

void Tourism::resize_border() {
    index_t n_border_cells = m_cell_router.send_buffer_size();
    index_t n_border_faces = m_face_router.send_buffer_size();
    index_t n_border_nodes = m_node_router.send_buffer_size();

    m_border.resize(n_border_cells, n_border_faces, n_border_nodes);
}

void Tourism::extend_border() {
    index_t n_border_cells = m_cell_router.send_buffer_size();
    index_t n_border_faces = m_face_router.send_buffer_size();
    index_t n_border_nodes = m_node_router.send_buffer_size();

    if (n_border_cells > m_border.size()) {
        m_border.resize(n_border_cells, n_border_faces, n_border_nodes);
    }
}

void Tourism::resize_aliens() {
    int n_alien_cells = m_cell_router.recv_buffer_size();
    int n_alien_faces = m_face_router.recv_buffer_size();
    int n_alien_nodes = m_node_router.recv_buffer_size();

    m_aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);
}

void Tourism::extend_aliens() {
    int n_alien_cells = m_cell_router.recv_buffer_size();
    int n_alien_faces = m_face_router.recv_buffer_size();
    int n_alien_nodes = m_node_router.recv_buffer_size();

    if (n_alien_cells > m_aliens.size()) {
        m_aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);
    }
}

void Tourism::fill_send_count(const AmrCells& locals) {
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    std::vector<index_t> last_append(size, -1);

    std::vector<index_t> cell_send_count(size, 0);
    std::vector<index_t> face_send_count(size, 0);
    std::vector<index_t> node_send_count(size, 0);

    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (locals.is_undefined(ic)) continue;

        for (index_t iface: locals.faces_range(ic)) {
            if (locals.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = locals.faces.adjacent.rank[iface];
            if (neib_rank != rank && last_append[neib_rank] != ic) {
                last_append[neib_rank] = ic;
                cell_send_count[neib_rank] += 1;
                face_send_count[neib_rank] += locals.max_face_count(ic);
                node_send_count[neib_rank] += locals.max_node_count(ic);
            }
        }
    }

    // Установить число на обмены
    m_cell_router.set_send_count(cell_send_count);
    m_face_router.set_send_count(face_send_count);
    m_node_router.set_send_count(node_send_count);
}

void Tourism::fill_indices(const AmrCells& locals) {
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    index_t last_unique_append = -1;
    std::vector<index_t> last_append(mpi::size(), -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> cell_index = m_cell_router.send_offset();

    m_border_indices.resize(m_cell_router.send_buffer_size());
    // m_unique_border_indices.reserve(m_cell_route.send_buffer_size());

    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (locals.is_undefined(ic)) continue;

        for (index_t iface: locals.faces_range(ic)) {
            if (locals.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = locals.faces.adjacent.rank[iface];
            if (neib_rank != rank) {
                if (last_append[neib_rank] != ic) {
                    last_append[neib_rank] = ic;
                    m_border_indices[cell_index[neib_rank]++] = ic;
                }
                if (last_unique_append != ic) {
                    last_unique_append = ic;
                    // m_unique_border_indices.push_back(ic);
                }
            }
        }
    }
}

void Tourism::prepare_geometry(const AmrCells& locals) {
    index_t face_idx = 0;
    index_t node_idx = 0;
    for (index_t ic = 0; ic < m_border_indices.size(); ++ic) {
        locals.copy_geom(m_border_indices[ic], m_border, ic, face_idx, node_idx);

        face_idx += locals.max_face_count(m_border_indices[ic]);
        node_idx += locals.max_node_count(m_border_indices[ic]);
    }
}

void Tourism::build_border(const AmrCells& locals) {
    // Заполняем router.send_count
    fill_send_count(locals);

    // Заполняем индексы
    fill_indices(locals);

    // Подготовим массив border
    resize_border();

    // Перенести геометрию из locals в border
    prepare_geometry(locals);
}

// Инициализация индекса alien = -1 для большинства граней
void set_undef_aliens(AmrCells& locals, int rank) {
    threads::parallel_for(
        index_t{0}, locals.n_cells(),
        [&locals, rank](index_t ic) {
            for (index_t iface: locals.faces_range(ic)) {
                if (locals.faces.is_actual(iface) &&
                    locals.faces.adjacent.rank[iface] == rank) {
                    locals.faces.adjacent.alien[iface] = -1;
                }
            }
        });
}

// Обходим ячейки в alien и ищем связи
void Tourism::find_connections(AmrCells& locals, int rank) const {
    for (index_t ic = 0; ic < m_aliens.n_cells(); ++ic) {
        for (index_t iface: m_aliens.faces_range(ic)) {
            if (m_aliens.faces.is_undefined(iface)) {
                continue;
            }

            if (m_aliens.faces.adjacent.index[iface] >= 0 &&
                m_aliens.faces.adjacent.rank[iface] == rank) {
                // Индекс соседа
                index_t jc = m_aliens.faces.adjacent.index[iface];

                for (index_t l_face: locals.faces_range(jc)) {
                    if (locals.faces.adjacent.rank [l_face] == m_aliens.rank [ic] &&
                        locals.faces.adjacent.index[l_face] == m_aliens.index[ic]) {

                        locals.faces.adjacent.alien[l_face] = ic;
                        break;
                        }
                }
            }
        }
    }
}

void Tourism::update(AmrCells& locals) {
    // Построить border-слой
    build_border(locals);

    // Заполнить recv массивы
    m_cell_router.fill_partial();
    m_face_router.fill_partial();
    m_node_router.fill_partial();

    // Расширить массив aliens для получения геометрии
    resize_aliens();

    // Отправить и получить геометрию
    sync_geometry();

    // Инициализация индекса alien = -1 для большинства граней
    set_undef_aliens(locals, mpi::rank());

    // Обходим ячейки в alien и ищем связи
    find_connections(locals, mpi::rank());
}

/// @brief Определить дочерние ячейки, которые прилегают к граням с рангом rank.
/// То есть дочерние ячейки, которые окажутся в border-блоке ранга rank.
/// @param faces Список граней ячеек из border-слоя
/// @param ic Индекс родительской ячейки
/// @param rank Ранг процесса, к которому ищется прилегание.
/// @return bitset<8> - true/false, прилегает дочерняя ячейка или нет.
template<int dim>
std::bitset<8> border_children(const AmrFaces& faces, index_t ic, int rank) {
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
        // i - индекс в m_border, m_border_indices
        for (index_t i: m_cell_router.send_indices(r)) {
            if (m_border.flag[i] == 0) {
                m_border.next[i] = next_index;
                next_index += 1;
            }
            else if (m_border.flag[i] == 1) {
                // bitset<8> для дочерних ячеек
                auto children = dim == 2 ? border_children<2>(m_border.faces, i, r) :
                                           border_children<3>(m_border.faces, i, r);
                // Кодируем список дочерних ячеек
                m_border.next[i] = amr::pack_children(next_index, children);
                next_index += static_cast<index_t>(children.count());
            }
            else {
                // Самый неприятный случай, border-ячейка огрубляется

                // Полный индекс родительской ячейки
                std::tuple<index_t, index_t, index_t> parent = {
                    m_border.b_idx[i],
                    m_border.level[i],
                    m_border.z_idx[i] / indexing::CpC(dim),
                };

                auto parent_it = coarse_cells.find(parent);
                if (parent_it != coarse_cells.end()) {
                    m_border.next[i] = parent_it->second;
                }
                else {
                    coarse_cells[parent] = next_index;
                    m_border.next[i] = next_index;
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
    std::vector<index_t> prev_border_indices = m_border_indices;
    m_border_indices.resize(m_cell_router.send_buffer_size());

    index_t last_border_next = 0;
    for (index_t i = 0; i < prev_border_indices.size(); ++i) {
        z_assert(i < m_border.flag.size(), "out of range #1521");
        z_assert(i < m_border.next.size(), "out of range #1522");
        z_assert(i < prev_border_indices.size(), "out of range #1523");

        if (m_border.flag[i] == 0) {
            index_t border_next = m_border.next[i];

            z_assert(border_next < m_border_indices.size(), "out of range #1524");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1525");

            m_border_indices[border_next] = locals_next[prev_border_indices[i]];
            last_border_next = std::max(last_border_next, border_next);
        }
        else if (m_border.flag[i] < 0) {
            index_t border_next = m_border.next[i];

            z_assert(border_next < m_border_indices.size(), "out of range #1526");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1527");

            index_t parent_index = locals_next[prev_border_indices[i]];

            z_assert(parent_index < locals_next.size(), "out of range #1528");

            m_border_indices[border_next] = locals_next[parent_index];
            last_border_next = std::max(last_border_next, border_next);
        }
        else {
            auto [border_next, children] = amr::unpack_children(m_border.next[i]);
            index_t main_child = locals_next[prev_border_indices[i]];
            for (int c = 0; c < indexing::CpC(dim); ++c) {
                if (children[c]) {
                    if (border_next >= m_border_indices.size()) {
                        std::cout << m_border.next[i] << "; " << border_next << "; " << children << "; " << m_border_indices.size() << "\n";
                    }
                    z_assert(border_next < m_border_indices.size(), "out of range #1529");
                    z_assert(main_child + c < locals_next.size(), "out of range #1530");
                    m_border_indices[border_next] = locals_next[main_child + c];
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
    //              Посчитаем смещения для новых border и aliens
    // ========================================================================

    std::vector<index_t> n_block_faces(mpi::size(), 0);
    std::vector<index_t> n_block_nodes(mpi::size(), 0);
    for (int r = 0; r < mpi::size(); ++r) {
        n_block_faces[r] = (dim == 2 ? 8 : 24) * n_block_cells[r];
        n_block_nodes[r] = (dim == 2 ? 9 : 27) * n_block_cells[r];
    }

    // Получим значения NEXT, далее можем менять роутеры
    send_next.wait();
    recv_next.wait();

#if WRITE_DBG
    static size_t pvd_counter = 0;
    static io::Variables vars = {"flag", "next", "rank", "level", "index", "b_idx", "z_idx"};
    static io::PvdFile bef_border("sp_border_bef", "debug");
    static io::PvdFile aft_border("sp_border_aft", "debug");
    static io::PvdFile bef_aliens("sp_aliens_bef", "debug");
    static io::PvdFile aft_aliens("sp_aliens_aft", "debug");

    if (pvd_counter == 0) {
        bef_border.variables = vars;
        bef_aliens.variables = vars;
        aft_border.variables = vars;
        aft_aliens.variables = vars;
    }

    bef_border.save(m_border, pvd_counter);
    bef_aliens.save(m_aliens, pvd_counter);
#endif

    Router prev_router = m_cell_router;

    // Установить число на отправку
    m_cell_router.set_send_count(n_block_cells);
    m_face_router.set_send_count(n_block_faces);
    m_node_router.set_send_count(n_block_nodes);

    // Заполнить recv массивы
    m_cell_router.fill_partial();
    m_face_router.fill_partial();
    m_node_router.fill_partial();

    // ========================================================================
    //          Сделаем глобальную индексацию next в border и aliens
    // ========================================================================

    // Добавляем смещения, теперь индексы NEXT в border идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // border-ячейки указана следующая позиция внутри нового border-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.send_indices(r)) {
            m_border.next[i] += m_cell_router.send_offset(r);
        }
    }

    // Добавляем смещения, теперь индексы NEXT в alien идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // alien-ячейки указана следующая позиция внутри нового alien-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.recv_indices(r)) {
            m_aliens.next[i] += m_cell_router.recv_offset(r);
        }
    }

#if WRITE_DBG
    aft_border.save(m_border, pvd_counter);
    aft_aliens.save(m_aliens, pvd_counter);
    ++pvd_counter;
#endif

    // ========================================================================
    //           Подготовим новые массивы, без актуальных данных
    // ========================================================================

    // Подготовить border массив
    extend_border();

    // Подготовить aliens массив
    extend_aliens();

    // Выставить корректные индексы в m_border_indices
    update_border_indices<dim>(locals_next);
}

template void Tourism::setup_positions<2>(const std::vector<index_t>&);
template void Tourism::setup_positions<3>(const std::vector<index_t>&);

void Tourism::pack_border_indices() {
    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в face_begin и node_begin количество элементов на ячейку
    for (index_t ic = 0; ic < m_border.size(); ++ic) {
        m_border.face_begin[ic] = m_border.face_begin[ic + 1] - m_border.face_begin[ic];
        m_border.node_begin[ic] = m_border.node_begin[ic + 1] - m_border.node_begin[ic];
    }
    m_border.face_begin.back() = -1;
    m_border.node_begin.back() = -1;
}

void Tourism::unpack_border_indices() {
    int prev_n_face = m_border.face_begin[0];
    int prev_n_node = m_border.node_begin[0];

    m_border.face_begin[0] = 0;
    m_border.node_begin[0] = 0;

    for (index_t ic = 1; ic <= m_border.n_cells(); ++ic) {
        int temp_n_face = m_border.face_begin[ic];
        int temp_n_node = m_border.node_begin[ic];

        m_border.face_begin[ic] = m_border.face_begin[ic - 1] + prev_n_face;
        m_border.node_begin[ic] = m_border.node_begin[ic - 1] + prev_n_node;

        prev_n_face = temp_n_face;
        prev_n_node = temp_n_node;
    }
}

void Tourism::unpack_aliens_indices() {
    // Восстанавливаем индексацию, сейчас в массивах хранится число граней или вершин ячейки
    m_aliens.face_begin[0] = 0;
    m_aliens.node_begin[0] = 0;
    for (index_t ic = 0; ic < m_aliens.n_cells(); ++ic) {
        m_aliens.face_begin[ic + 1] += m_aliens.face_begin[ic];
        m_aliens.node_begin[ic + 1] += m_aliens.node_begin[ic];
    }
}

template<int dim>
void set_amr_indices(std::vector<index_t>& faces_beg, std::vector<index_t>& nodes_beg) {
    z_assert(faces_beg.size() == nodes_beg.size(), "restore amr sizes mismatch");

    constexpr int n_faces = Side<dim>::n_subfaces();
    constexpr int n_nodes = dim == 2 ? 9 : 27;

    threads::parallel_for(
        index_t{0}, index_t(faces_beg.size()),
        [&faces_beg, &nodes_beg](index_t ic) {
            faces_beg[ic] = n_faces * ic;
            nodes_beg[ic] = n_nodes * ic;
        });
}

void Tourism::send_geometry(const AmrCells& locals) {
    prepare_geometry(locals);
    sync_geometry();
}

void Tourism::restore_indices(AmrCells& locals) const {
    for (index_t ic: m_border_indices) {
        for (index_t iface: locals.faces_range(ic)) {
            index_t alien_index = locals.faces.adjacent.alien[iface];
            if (alien_index >= 0) {
                locals.faces.adjacent.index[iface] = m_aliens.index[alien_index];
            }
        }
    }
}

void Tourism::sync_geometry() {
    bool amr = m_border.adaptive();
    bool axial = m_border.axial();

    if (!amr) {
        // Оптимизируем пересылку индексов граней/вершин
        pack_border_indices();
    }

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    RequestsList cells_send; cells_send.reserve(16);
    cells_send += m_cell_router.isend(m_border.rank, MpiTag::RANK);
    cells_send += m_cell_router.isend(m_border.next, MpiTag::NEXT);
    cells_send += m_cell_router.isend(m_border.index, MpiTag::INDEX);
    cells_send += m_cell_router.isend(m_border.flag, MpiTag::FLAG);
    cells_send += m_cell_router.isend(m_border.level, MpiTag::LEVEL);
    cells_send += m_cell_router.isend(m_border.b_idx, MpiTag::B_IDX);
    cells_send += m_cell_router.isend(m_border.z_idx, MpiTag::Z_IDX);
    cells_send += m_cell_router.isend(m_border.center, MpiTag::CENTER);
    cells_send += m_cell_router.isend(m_border.volume, MpiTag::VOLUME);
    if (axial) {
        cells_send += m_cell_router.isend(m_border.volume_alt, MpiTag::VOLUME_ALT);
    }
    if (!amr) {
        cells_send += m_cell_router.isend(m_border.face_begin, MpiTag::FACE_BEG);
        cells_send += m_cell_router.isend(m_border.node_begin, MpiTag::NODE_BEG);
    }

    // Отправить данные граней
    RequestsList faces_send; faces_send.reserve(16);
    faces_send += m_face_router.isend(m_border.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_send += m_face_router.isend(m_border.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_send += m_face_router.isend(m_border.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    faces_send += m_face_router.isend(m_border.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_send += m_face_router.isend(m_border.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_send += m_face_router.isend(m_border.faces.boundary, MpiTag::BOUNDARY);
    faces_send += m_face_router.isend(m_border.faces.normal, MpiTag::NORMAL);
    faces_send += m_face_router.isend(m_border.faces.center, MpiTag::FACE_CENTER);
    faces_send += m_face_router.isend(m_border.faces.area, MpiTag::AREA);
    if (axial) {
        faces_send += m_face_router.isend(m_border.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_send += m_face_router.isend(m_border.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    auto nodes_send = m_node_router.isend(m_border.verts, MpiTag::VERTICES);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += m_cell_router.irecv(m_aliens.rank, MpiTag::RANK);
    cells_recv += m_cell_router.irecv(m_aliens.next, MpiTag::NEXT);
    cells_recv += m_cell_router.irecv(m_aliens.index, MpiTag::INDEX);
    cells_recv += m_cell_router.irecv(m_aliens.flag, MpiTag::FLAG);
    cells_recv += m_cell_router.irecv(m_aliens.level, MpiTag::LEVEL);
    cells_recv += m_cell_router.irecv(m_aliens.b_idx, MpiTag::B_IDX);
    cells_recv += m_cell_router.irecv(m_aliens.z_idx, MpiTag::Z_IDX);
    cells_recv += m_cell_router.irecv(m_aliens.center, MpiTag::CENTER);
    cells_recv += m_cell_router.irecv(m_aliens.volume, MpiTag::VOLUME);
    if (axial) {
        cells_recv += m_cell_router.irecv(m_aliens.volume_alt, MpiTag::VOLUME_ALT);
    }

    if (!amr) {
        // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
        cells_recv += m_cell_router.irecv(m_aliens.face_begin.data() + 1, MpiTag::FACE_BEG);
        cells_recv += m_cell_router.irecv(m_aliens.node_begin.data() + 1, MpiTag::NODE_BEG);
    }

    // Получить данные граней
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += m_face_router.irecv(m_aliens.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += m_face_router.irecv(m_aliens.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += m_face_router.irecv(m_aliens.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    faces_recv += m_face_router.irecv(m_aliens.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += m_face_router.irecv(m_aliens.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_recv += m_face_router.irecv(m_aliens.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += m_face_router.irecv(m_aliens.faces.normal, MpiTag::NORMAL);
    faces_recv += m_face_router.irecv(m_aliens.faces.center, MpiTag::FACE_CENTER);
    faces_recv += m_face_router.irecv(m_aliens.faces.area, MpiTag::AREA);
    if (axial) {
        faces_recv += m_face_router.irecv(m_aliens.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_recv += m_face_router.irecv(m_aliens.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    auto nodes_recv = m_node_router.irecv(m_aliens.verts, MpiTag::VERTICES);

    // =========================== WAIT ISEND =================================

    cells_send.wait();  // Завершить отправку ячеек
    faces_send.wait();  // Завершить отправку граней
    nodes_send.wait();  // Завершить отправку вершин

    // =========================== WAIT IRECV =================================

    cells_recv.wait();  // Завершить получение ячеек
    faces_recv.wait();  // Завершить получение граней
    nodes_recv.wait();  // Завершить получение вершин

    if (!amr) {
        // Восстановить индексацию граней
        unpack_aliens_indices();

        // Поддерживать индексацию граней в border массиве не обязательно
        // unpack_border_indices();
    }
    else {
        if (m_border.dim() == 2) {
            set_amr_indices<2>(m_aliens.face_begin, m_aliens.node_begin);
        }
        else {
            set_amr_indices<3>(m_aliens.face_begin, m_aliens.node_begin);
        }
    }
}

} // namespace zephyr::mesh

#endif