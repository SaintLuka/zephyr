#include <bitset>
#include <numeric>
#include <map>
#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/mesh/euler/router.h>

#include <zephyr/io/pvd_file.h>

#ifdef ZEPHYR_MPI

namespace zephyr::mesh {

using utils::mpi;

void Tourism::shrink_to_fit() {
    //m_unique_border_indices.shrink_to_fit();
    m_border_indices.shrink_to_fit();
    m_border.shrink_to_fit();
}

void Tourism::init_types(const AmrCells& locals) {
    m_border = locals.same();
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
    m_cell_route.set_send_count(cell_send_count);
    m_face_route.set_send_count(face_send_count);
    m_node_route.set_send_count(node_send_count);
}

void Tourism::fill_indices(const AmrCells& locals) {
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    index_t last_unique_append = -1;
    std::vector<index_t> last_append(size, -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> cell_index = m_cell_route.send_offset();

    m_border_indices.resize(m_cell_route.send_buffer_size());
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

void Tourism::pack_geometry(const AmrCells& locals) {
    index_t face_idx = 0;
    index_t node_idx = 0;
    for (index_t ic = 0; ic < m_border_indices.size(); ++ic) {
        locals.copy_geom(m_border_indices[ic], m_border, ic, face_idx, node_idx);

        face_idx += locals.max_face_count(m_border_indices[ic]);
        node_idx += locals.max_node_count(m_border_indices[ic]);
    }

    if (mpi::rank() == 1) {
        //std::cout << "PACK GEOM #0: " << m_border.size() << ", " << m_border.face_begin.size() << "; bs size: " << m_border_indices.size() << "\n";
        //for (auto i: m_border.face_begin) std::cout << i << " ";
        //std::cout << "\n";
    }
}

void Tourism::build_border(const AmrCells& locals) {
    // Заполняем router.send_count
    fill_send_count(locals);

    // Заполняем индексы
    fill_indices(locals);

    // Подготовим массив border
    index_t n_border_cells = m_cell_route.send_buffer_size();
    index_t n_border_faces = m_face_route.send_buffer_size();
    index_t n_border_nodes = m_node_route.send_buffer_size();

    m_border.resize(n_border_cells, n_border_faces, n_border_nodes);

    pack_geometry(locals);
}

void Tourism::build_border_basic(const AmrCells& locals) {
    // Заполняем router.send_count
    fill_send_count(locals);

    // Заполняем индексы
    fill_indices(locals);

    // Подготовим массив border
    index_t n_border_cells = m_cell_route.send_buffer_size();
    index_t n_border_faces = m_face_route.send_buffer_size();
    index_t n_border_nodes = m_node_route.send_buffer_size();

    m_border.resize(n_border_cells, n_border_faces, n_border_nodes);

    pack_geometry(locals);
}

void Tourism::build_aliens(AmrCells& locals, AmrCells& aliens) {
    // Построить border-слой
    build_border(locals);

    // Заполнить recv массивы
    m_cell_route.fill_partial();
    m_face_route.fill_partial();
    m_node_route.fill_partial();

    // Подготовим массив aliens для получения геометрии
    int n_alien_cells = m_cell_route.recv_buffer_size();
    int n_alien_faces = m_face_route.recv_buffer_size();
    int n_alien_nodes = m_node_route.recv_buffer_size();
    
    aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);

    send_geometry(aliens);

    // ========================================================================

    const int rank = mpi::rank();

    // Инициализация индекса alien = -1 для большинства граней
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        for (index_t iface: locals.faces_range(ic)) {
            if (locals.faces.is_actual(iface) &&
                locals.faces.adjacent.rank[iface] == rank) {
                locals.faces.adjacent.alien[iface] = -1;
            }
        }
    }

    // Обходим ячейки в alien и ищем связи
    for (index_t ic = 0; ic < aliens.n_cells(); ++ic) {
        for (index_t iface: aliens.faces_range(ic)) {
            if (aliens.faces.is_undefined(iface)) {
                continue;
            }

            if (aliens.faces.adjacent.index[iface] >= 0 &&
                aliens.faces.adjacent.rank[iface] == rank) {

                // Индекс соседа
                index_t jc = aliens.faces.adjacent.index[iface];

                for (index_t l_face: locals.faces_range(jc)) {
                    if (locals.faces.adjacent.rank [l_face] == aliens.rank [ic] &&
                        locals.faces.adjacent.index[l_face] == aliens.index[ic]) {

                        locals.faces.adjacent.alien[l_face] = ic;
                        break;
                    }
                }
            }
        }
    }
}

void Tourism::build_aliens_basic(AmrCells& locals, AmrCells& aliens) {
    // Построить border-слой
    build_border_basic(locals);

    // Заполнить recv массивы
    m_cell_route.fill_partial();
    m_face_route.fill_partial();
    m_node_route.fill_partial();

    // Подготовим массив aliens для получения геометрии
    int n_alien_cells = m_cell_route.recv_buffer_size();
    int n_alien_faces = m_face_route.recv_buffer_size();
    int n_alien_nodes = m_node_route.recv_buffer_size();

    aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);

    // ========================== Пересылка ячеек =============================

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в face_begin и node_begin количество элементов на ячейку
    for (index_t ic = 0; ic < m_border.size(); ++ic) {
        m_border.face_begin[ic] = m_border.face_begin[ic + 1] - m_border.face_begin[ic];
        m_border.node_begin[ic] = m_border.node_begin[ic + 1] - m_border.node_begin[ic];
    }

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    auto send_rank     = m_cell_route.isend(m_border.rank, MpiTag::RANK);
    auto send_index    = m_cell_route.isend(m_border.index, MpiTag::INDEX);
    auto send_face_beg = m_cell_route.isend(m_border.face_begin, MpiTag::FACE_BEG);
    auto send_node_beg = m_cell_route.isend(m_border.node_begin, MpiTag::NODE_BEG);

    // Отправить данные граней
    auto send_adj_rank  = m_face_route.isend(m_border.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto send_adj_index = m_face_route.isend(m_border.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto send_adj_alien = m_face_route.isend(m_border.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto send_boundary  = m_face_route.isend(m_border.faces.boundary, MpiTag::BOUNDARY);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    auto recv_rank     = m_cell_route.irecv(aliens.rank, MpiTag::RANK);
    auto recv_index    = m_cell_route.irecv(aliens.index, MpiTag::INDEX);

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    auto recv_face_beg = m_cell_route.irecv(aliens.face_begin.data() + 1, MpiTag::FACE_BEG);
    auto recv_node_beg = m_cell_route.irecv(aliens.node_begin.data() + 1, MpiTag::NODE_BEG);

    // Получить данные граней
    auto recv_adj_rank  = m_face_route.irecv(aliens.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto recv_adj_index = m_face_route.irecv(aliens.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto recv_adj_alien = m_face_route.irecv(aliens.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto recv_boundary  = m_face_route.irecv(aliens.faces.boundary, MpiTag::BOUNDARY);

    // =========================== WAIT ISEND =================================

    // Завершить отправку данных ячеек
    send_rank.wait();
    send_index.wait();
    send_face_beg.wait();
    send_node_beg.wait();

    // Завершить отправку данных граней
    send_adj_rank.wait();
    send_adj_index.wait();
    send_adj_alien.wait();
    send_boundary.wait();

    // =========================== WAIT IRECV =================================

    // Завершить получение данных ячеек
    recv_rank.wait();
    recv_index.wait();
    recv_face_beg.wait();
    recv_node_beg.wait();

    // Завершить получение данных граней
    recv_adj_rank.wait();
    recv_adj_index.wait();
    recv_adj_alien.wait();
    recv_boundary.wait();

    // Восстанавливаем индексацию, сейчас в массивах хранится число граней или вершин ячейки
    aliens.face_begin[0] = 0;
    aliens.node_begin[0] = 0;
    for (index_t ic = 0; ic < aliens.n_cells(); ++ic) {
        aliens.face_begin[ic + 1] += aliens.face_begin[ic];
        aliens.node_begin[ic + 1] += aliens.node_begin[ic];
    }

    m_border.face_begin[0] = 0;
    m_border.node_begin[0] = 0;
    for (index_t ic = 0; ic < m_border.n_cells(); ++ic) {
        m_border.face_begin[ic + 1] += m_border.face_begin[ic];
        m_border.node_begin[ic + 1] += m_border.node_begin[ic];
    }

    // ========================================================================

    const int rank = mpi::rank();

    // Инициализация индекса alien = -1 для большинства граней
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        for (index_t iface: locals.faces_range(ic)) {
            if (locals.faces.is_actual(iface) &&
                locals.faces.adjacent.rank[iface] == rank) {
                locals.faces.adjacent.alien[iface] = -1;
            }
        }
    }

    // Обходим ячейки в alien и ищем связи
    for (index_t ic = 0; ic < aliens.n_cells(); ++ic) {
        for (index_t iface: aliens.faces_range(ic)) {
            if (aliens.faces.is_undefined(iface)) {
                continue;
            }

            if (aliens.faces.adjacent.index[iface] >= 0 &&
                aliens.faces.adjacent.rank[iface] == rank) {

                // Индекс соседа
                index_t jc = aliens.faces.adjacent.index[iface];

                for (index_t l_face: locals.faces_range(jc)) {
                    if (locals.faces.adjacent.rank [l_face] == aliens.rank [ic] &&
                        locals.faces.adjacent.index[l_face] == aliens.index[ic]) {

                        locals.faces.adjacent.alien[l_face] = ic;
                        break;
                    }
                }
            }
        }
    }
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
                for (int i: side.children()) {
                    children[i] = true;
                }
            }
        }
        else { // Complex Face
            for (auto subface: side.subfaces()) {
                index_t iface = face_beg + subface;
                if (faces.adjacent.rank[iface] == rank) {
                    children[subface.child()] = true;
                }
            }
        }
    }
    return children;
}

index_t pack_children(index_t next, std::bitset<8> children) {
    z_assert(next < (1 << 23), "Too large next");
    return static_cast<index_t>((children.to_ulong() << 24) + next);
}

std::tuple<index_t, std::bitset<8>> unpack_children(index_t next) {
    return {static_cast<std::make_unsigned_t<index_t>>(next) % (1 << 24), std::bitset<8>(next >> 24) };
}

void Tourism::resize_to_router(AmrCells& aliens) {
    // Подготовим массив border
    index_t n_border_cells = m_cell_route.send_buffer_size();
    index_t n_border_faces = m_face_route.send_buffer_size();
    index_t n_border_nodes = m_node_route.send_buffer_size();

    m_border.resize(n_border_cells, n_border_faces, n_border_nodes);

    // Подготовим массив aliens для получения геометрии
    index_t n_alien_cells = m_cell_route.recv_buffer_size();
    index_t n_alien_faces = m_face_route.recv_buffer_size();
    index_t n_alien_nodes = m_node_route.recv_buffer_size();

    aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);
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
        for (index_t i: m_cell_route.send_indices(r)) {
            if (m_border.flag[i] == 0) {
                m_border.next[i] = next_index;
                next_index += 1;
            }
            else if (m_border.flag[i] == 1) {
                // bitset<8> для дочерних ячеек
                auto children = dim == 2 ? border_children<2>(m_border.faces, i, r) :
                                           border_children<3>(m_border.faces, i, r);
                // Кодируем список дочерних ячеек
                m_border.next[i] = pack_children(next_index, children);
                next_index += static_cast<index_t>(children.count());
            }
            else {
                // Самый неприятный случай, border-ячейка огрубляется

                // Полный индекс родительской ячейки
                std::tuple<index_t, index_t, index_t> parent = {
                    m_border.b_idx[i],
                    m_border.level[i],
                    m_border.z_idx[i] / CpC(dim),
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

#define WRITE_DBG 0

template<int dim>
void Tourism::update_border_indices(const std::vector<index_t>& locals_next) {
    std::vector<index_t> prev_border_indices = m_border_indices;
    m_border_indices.resize(m_cell_route.send_buffer_size());

    index_t last_border_next = 0;
    for (index_t i = 0; i < prev_border_indices.size(); ++i) {
        z_assert(i < m_border.flag.size(), "out of range #1521");
        z_assert(i < m_border.next.size(), "out of range #1522");
        z_assert(i < old_border_indices.size(), "out of range #1523");

        if (m_border.flag[i] == 0) {
            index_t border_next = m_border.next[i];

            z_assert(border_next < m_border_indices.size(), "out of range #1524");
            z_assert(old_border_indices[i] < locals.next.size(), "out of range #1525");

            m_border_indices[border_next] = locals_next[prev_border_indices[i]];
            last_border_next = std::max(last_border_next, border_next);
        }
        else if (m_border.flag[i] < 0) {
            index_t border_next = m_border.next[i];

            z_assert(border_next < m_border_indices.size(), "out of range #1526");
            z_assert(old_border_indices[i] < locals.next.size(), "out of range #1527");

            index_t parent_index = locals_next[prev_border_indices[i]];

            z_assert(parent_index < locals_next.size(), "out of range #1528");

            m_border_indices[border_next] = locals_next[parent_index];
            last_border_next = std::max(last_border_next, border_next);
        }
        else {
            auto [border_next, children] = unpack_children(m_border.next[i]);
            index_t main_child = locals_next[prev_border_indices[i]];
            for (int c = 0; c < CpC(dim); ++c) {
                if (children[c]) {
                    if (border_next >= m_border_indices.size()) {
                        std::cout << m_border.next[i] << "; " << border_next << "; " << children << "; " << m_border_indices.size() << "\n";
                    }
                    z_assert(border_next < m_border_indices.size(), "out of range #1529");
                    z_assert(main_child + c < locals.next.size(), "out of range #1530");
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

template<int dim>
void Tourism::setup_positions(AmrCells& aliens) {
    int rank = mpi::rank();

    // Размеры border-блоков / новое число ячеек на отправку
    auto n_block_cells = setup_border_next<dim>();

    // Отправим значения NEXT, последнее использование старого роутера
    auto send_next = isend<MpiTag::NEXT>();
    auto recv_next = irecv<MpiTag::NEXT>(aliens);

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

#ifdef WRITE_DBG
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
    bef_aliens.save(aliens, pvd_counter);
#endif

    Router prev_router = m_cell_route;

    // Установить число на отправку
    m_cell_route.set_send_count(n_block_cells);
    m_face_route.set_send_count(n_block_faces);
    m_node_route.set_send_count(n_block_nodes);

    // Заполнить recv массивы
    m_cell_route.fill_partial();
    m_face_route.fill_partial();
    m_node_route.fill_partial();

    // ========================================================================
    //          Сделаем глобальную индексацию next в border и aliens
    // ========================================================================

    // Добавляем смещения, теперь индексы NEXT в border идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // border-ячейки указана следующая позиция внутри нового border-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.send_indices(r)) {
            m_border.next[i] += m_cell_route.send_offset(r);
        }
    }

    // Добавляем смещения, теперь индексы NEXT в alien идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // alien-ячейки указана следующая позиция внутри нового alien-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.recv_indices(r)) {
            aliens.next[i] += m_cell_route.recv_offset(r);
        }
    }

#if WRITE_DBG
    aft_border.save(m_border, pvd_counter);
    aft_aliens.save(aliens, pvd_counter);
    ++pvd_counter;
#endif

    // ========================================================================
    //           Подготовим новые массивы, без актуальных данных
    // ========================================================================

    // Подготовим массив border
    index_t n_border_cells = m_cell_route.send_buffer_size();
    index_t n_border_faces = m_face_route.send_buffer_size();
    index_t n_border_nodes = m_node_route.send_buffer_size();

    if (n_border_cells > prev_router.send_buffer_size()) {
        m_border.resize(n_border_cells, n_border_faces, n_border_nodes);
    }

    // Подготовим массив aliens для получения геометрии
    index_t n_alien_cells = m_cell_route.recv_buffer_size();
    index_t n_alien_faces = m_face_route.recv_buffer_size();
    index_t n_alien_nodes = m_node_route.recv_buffer_size();

    if (n_alien_cells > prev_router.recv_buffer_size()) {
        aliens.resize(n_alien_cells, n_alien_faces, n_alien_nodes);
    }
}

template void Tourism::setup_positions<2>(AmrCells& aliens);
template void Tourism::setup_positions<3>(AmrCells& aliens);

void Tourism::send_geometry(AmrCells& aliens) {
    // ========================== Пересылка ячеек =============================

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в face_begin и node_begin количество элементов на ячейку
    for (index_t ic = 0; ic < m_border.size(); ++ic) {
        m_border.face_begin[ic] = m_border.face_begin[ic + 1] - m_border.face_begin[ic];
        m_border.node_begin[ic] = m_border.node_begin[ic + 1] - m_border.node_begin[ic];
    }
    m_border.face_begin.back() = -1;
    m_border.node_begin.back() = -1;

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    auto send_rank     = m_cell_route.isend(m_border.rank, MpiTag::RANK);
    auto send_next     = m_cell_route.isend(m_border.next, MpiTag::NEXT);
    auto send_index    = m_cell_route.isend(m_border.index, MpiTag::INDEX);
    auto send_flag     = m_cell_route.isend(m_border.flag, MpiTag::FLAG);
    auto send_level    = m_cell_route.isend(m_border.level, MpiTag::LEVEL);
    auto send_b_idx    = m_cell_route.isend(m_border.b_idx, MpiTag::B_IDX);
    auto send_z_idx    = m_cell_route.isend(m_border.z_idx, MpiTag::Z_IDX);
    auto send_center   = m_cell_route.isend(m_border.center, MpiTag::CENTER);
    auto send_volume   = m_cell_route.isend(m_border.volume, MpiTag::VOLUME);
    auto send_volume_a = m_cell_route.isend(m_border.volume_alt, MpiTag::VOLUME_ALT);
    auto send_face_beg = m_cell_route.isend(m_border.face_begin, MpiTag::FACE_BEG);
    auto send_node_beg = m_cell_route.isend(m_border.node_begin, MpiTag::NODE_BEG);

    // Отправить данные граней
    auto send_adj_rank  = m_face_route.isend(m_border.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto send_adj_index = m_face_route.isend(m_border.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto send_adj_alien = m_face_route.isend(m_border.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto send_adj_basic = m_face_route.isend(m_border.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto send_boundary  = m_face_route.isend(m_border.faces.boundary, MpiTag::BOUNDARY);
    auto send_normal    = m_face_route.isend(m_border.faces.normal, MpiTag::NORMAL);
    auto send_f_center  = m_face_route.isend(m_border.faces.center, MpiTag::FACE_CENTER);
    auto send_area      = m_face_route.isend(m_border.faces.area, MpiTag::AREA);
    auto send_area_alt  = m_face_route.isend(m_border.faces.area_alt, MpiTag::AREA_ALT);
    auto send_face_vs   = m_face_route.isend(m_border.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    auto send_vertices  = m_node_route.isend(m_border.verts, MpiTag::VERTICES);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    auto recv_rank     = m_cell_route.irecv(aliens.rank, MpiTag::RANK);
    auto recv_next     = m_cell_route.irecv(aliens.next, MpiTag::NEXT);
    auto recv_index    = m_cell_route.irecv(aliens.index, MpiTag::INDEX);
    auto recv_flag     = m_cell_route.irecv(aliens.flag, MpiTag::FLAG);
    auto recv_level    = m_cell_route.irecv(aliens.level, MpiTag::LEVEL);
    auto recv_b_idx    = m_cell_route.irecv(aliens.b_idx, MpiTag::B_IDX);
    auto recv_z_idx    = m_cell_route.irecv(aliens.z_idx, MpiTag::Z_IDX);
    auto recv_center   = m_cell_route.irecv(aliens.center, MpiTag::CENTER);
    auto recv_volume   = m_cell_route.irecv(aliens.volume, MpiTag::VOLUME);
    auto recv_volume_a = m_cell_route.irecv(aliens.volume_alt, MpiTag::VOLUME_ALT);

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    auto recv_face_beg = m_cell_route.irecv(aliens.face_begin.data() + 1, MpiTag::FACE_BEG);
    auto recv_node_beg = m_cell_route.irecv(aliens.node_begin.data() + 1, MpiTag::NODE_BEG);

    // Получить данные граней
    auto recv_adj_rank  = m_face_route.irecv(aliens.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto recv_adj_index = m_face_route.irecv(aliens.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto recv_adj_alien = m_face_route.irecv(aliens.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto recv_adj_basic = m_face_route.irecv(aliens.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto recv_boundary  = m_face_route.irecv(aliens.faces.boundary, MpiTag::BOUNDARY);
    auto recv_normal    = m_face_route.irecv(aliens.faces.normal, MpiTag::NORMAL);
    auto recv_f_center  = m_face_route.irecv(aliens.faces.center, MpiTag::FACE_CENTER);
    auto recv_area      = m_face_route.irecv(aliens.faces.area, MpiTag::AREA);
    auto recv_area_alt  = m_face_route.irecv(aliens.faces.area_alt, MpiTag::AREA_ALT);
    auto recv_face_vs   = m_face_route.irecv(aliens.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    auto recv_vertices  = m_node_route.irecv(aliens.verts, MpiTag::VERTICES);

    // =========================== WAIT ISEND =================================

    // Завершить отправку данных ячеек
    send_rank.wait();
    send_next.wait();
    send_index.wait();
    send_flag.wait();
    send_level.wait();
    send_b_idx.wait();
    send_z_idx.wait();
    send_center.wait();
    send_volume.wait();
    send_volume_a.wait();
    send_face_beg.wait();
    send_node_beg.wait();

    // Завершить отправку данных граней
    send_adj_rank.wait();
    send_adj_index.wait();
    send_adj_alien.wait();
    send_adj_basic.wait();
    send_boundary.wait();
    send_normal.wait();
    send_f_center.wait();
    send_area.wait();
    send_area_alt.wait();
    send_face_vs.wait();

    // Завершить отправку вершин
    send_vertices.wait();

    // =========================== WAIT IRECV =================================

    // Завершить получение данных ячеек
    recv_rank.wait();
    recv_next.wait();
    recv_index.wait();
    recv_flag.wait();
    recv_level.wait();
    recv_b_idx.wait();
    recv_z_idx.wait();
    recv_center.wait();
    recv_volume.wait();
    recv_volume_a.wait();
    recv_face_beg.wait();
    recv_node_beg.wait();

    // Завершить получение данных граней
    recv_adj_rank.wait();
    recv_adj_index.wait();
    recv_adj_alien.wait();
    recv_adj_basic.wait();
    recv_boundary.wait();
    recv_normal.wait();
    recv_f_center.wait();
    recv_area.wait();
    recv_area_alt.wait();
    recv_face_vs.wait();

    // Завершить получение вершин
    recv_vertices.wait();

    // Восстанавливаем индексацию, сейчас в массивах хранится число граней или вершин ячейки
    aliens.face_begin[0] = 0;
    aliens.node_begin[0] = 0;
    for (index_t ic = 0; ic < aliens.n_cells(); ++ic) {
        aliens.face_begin[ic + 1] += aliens.face_begin[ic];
        aliens.node_begin[ic + 1] += aliens.node_begin[ic];
    }

    // Восстанавливаем индексацию, сейчас в массивах хранится число граней или вершин ячейки
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

} // namespace zephyr::mesh

#endif