#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/mesh/euler/router.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>

using namespace zephyr::io;

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI

Tourism::Tourism() {
    m_border_indices.resize(mpi::size());
}

void Tourism::shrink_to_fit() {
    for (auto& indices: m_border_indices) {
        indices.shrink_to_fit();
    }
    m_border.shrink_to_fit();
}

void Tourism::init_types(const AmrCells& locals) {
    m_border = locals.same();
}

void Tourism::build_border(AmrCells& locals) {
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Очистить списки и смещения
    for (auto &indices: m_border_indices) {
        indices.clear();
    }

    // Заполняем m_border_indices
    auto &faces = locals.faces;
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        // build alien можно вызвать для не совсем нормальной сетки
        if (locals.is_undefined(ic)) {
            continue;
        }

        for (index_t iface: locals.faces_range(ic)) {
            if (faces.is_undefined(iface)) {
                continue;
            }

            auto &indices = m_border_indices[faces.adjacent.rank[iface]];
            if (faces.adjacent.rank[iface] != rank &&
                (indices.empty() || indices.back() != locals.index[ic])) {
                indices.push_back(locals.index[ic]);
            }
        }
    }

    {
        // Заполняем cell_send_count
        std::vector<index_t> cell_send_count(size);
        for (int r = 0; r < size; ++r) {
            cell_send_count[r] = m_border_indices[r].size();
        }

        // Заполняем face_send_count и node_send_count
        std::vector<index_t> face_send_count(size);
        std::vector<index_t> node_send_count(size);
        for (int r = 0; r < size; ++r) {
            auto &indices = m_border_indices[r];
            face_send_count[r] = 0;
            node_send_count[r] = 0;
            for (index_t ic: indices) {
                face_send_count[r] += locals.max_face_count(ic);
                node_send_count[r] += locals.max_node_count(ic);
            }
        }

        // Установить число на обмены
        cell_route.set_send_count(cell_send_count);
        face_route.set_send_count(face_send_count);
        node_route.set_send_count(node_send_count);
    }

    // Заполняем m_border
    index_t n_border_cells = cell_route.send_buffer_size();
    index_t n_border_faces = face_route.send_buffer_size();
    index_t n_border_nodes = node_route.send_buffer_size();

    // Подготовить массив для border
    m_border.resize(n_border_cells, n_border_faces, n_border_nodes);

    // Копируем геометрию в border слой
    index_t cell_idx = 0;
    index_t face_idx = 0;
    index_t node_idx = 0;
    for (auto &indices: m_border_indices) {
        for (index_t ic: indices) {
            locals.copy_geom(ic, m_border, cell_idx, face_idx, node_idx);

            cell_idx += 1;
            face_idx += locals.max_face_count(ic);
            node_idx += locals.max_node_count(ic);
        }
    }

    /*
    static int counter = 0;
    static PvdFile pvd("border");
    pvd.variables = {"rank", "index", "faces2D"};
    pvd.save(m_border, counter);
    ++counter;
     */
}

void Tourism::build_aliens(AmrCells& locals, AmrCells& aliens) {
    const int rank = mpi::rank();

    build_border(locals);

    // Заполнить recv массивы
    cell_route.fill_partial();
    face_route.fill_partial();
    node_route.fill_partial();

    // Подготовим массив aliens для получения геометрии
    int n_alien_cells = cell_route.recv_buffer_size();
    int n_alien_faces = face_route.recv_buffer_size();
    int n_alien_nodes = node_route.recv_buffer_size();
    
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
    auto send_rank     = cell_route.isend(m_border.rank, MpiTag::RANK);
    auto send_next     = cell_route.isend(m_border.next, MpiTag::NEXT);
    auto send_index    = cell_route.isend(m_border.index, MpiTag::INDEX);
    auto send_flag     = cell_route.isend(m_border.flag, MpiTag::FLAG);
    auto send_level    = cell_route.isend(m_border.level, MpiTag::LEVEL);
    auto send_b_idx    = cell_route.isend(m_border.b_idx, MpiTag::B_IDX);
    auto send_z_idx    = cell_route.isend(m_border.z_idx, MpiTag::Z_IDX);
    auto send_center   = cell_route.isend(m_border.center, MpiTag::CENTER);
    auto send_volume   = cell_route.isend(m_border.volume, MpiTag::VOLUME);
    auto send_volume_a = cell_route.isend(m_border.volume_alt, MpiTag::VOLUME_ALT);
    auto send_face_beg = cell_route.isend(m_border.face_begin, MpiTag::FACE_BEG);
    auto send_node_beg = cell_route.isend(m_border.node_begin, MpiTag::NODE_BEG);

    // Отправить данные граней
    auto send_adj_rank  = face_route.isend(m_border.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto send_adj_index = face_route.isend(m_border.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto send_adj_alien = face_route.isend(m_border.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto send_adj_basic = face_route.isend(m_border.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto send_boundary  = face_route.isend(m_border.faces.boundary, MpiTag::BOUNDARY);
    auto send_normal    = face_route.isend(m_border.faces.normal, MpiTag::NORMAL);
    auto send_f_center  = face_route.isend(m_border.faces.center, MpiTag::FACE_CENTER);
    auto send_area      = face_route.isend(m_border.faces.area, MpiTag::AREA);
    auto send_area_alt  = face_route.isend(m_border.faces.area_alt, MpiTag::AREA_ALT);
    auto send_face_vs   = face_route.isend(m_border.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    auto send_vertices  = node_route.isend(m_border.verts, MpiTag::VERTICES);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    auto recv_rank     = cell_route.irecv(aliens.rank,  MpiTag::RANK);
    auto recv_next     = cell_route.irecv(aliens.next,  MpiTag::NEXT);
    auto recv_index    = cell_route.irecv(aliens.index, MpiTag::INDEX);
    auto recv_flag     = cell_route.irecv(aliens.flag,  MpiTag::FLAG);
    auto recv_level    = cell_route.irecv(aliens.level, MpiTag::LEVEL);
    auto recv_b_idx    = cell_route.irecv(aliens.b_idx, MpiTag::B_IDX);
    auto recv_z_idx    = cell_route.irecv(aliens.z_idx, MpiTag::Z_IDX);
    auto recv_center   = cell_route.irecv(aliens.center, MpiTag::CENTER);
    auto recv_volume   = cell_route.irecv(aliens.volume, MpiTag::VOLUME);
    auto recv_volume_a = cell_route.irecv(aliens.volume_alt, MpiTag::VOLUME_ALT);

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    auto recv_face_beg = cell_route.irecv(aliens.face_begin.data() + 1, MpiTag::FACE_BEG);
    auto recv_node_beg = cell_route.irecv(aliens.node_begin.data() + 1, MpiTag::NODE_BEG);
    
    // Получить данные граней
    auto recv_adj_rank  = face_route.irecv(aliens.faces.adjacent.rank,  MpiTag::ADJ_RANK);
    auto recv_adj_index = face_route.irecv(aliens.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto recv_adj_alien = face_route.irecv(aliens.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto recv_adj_basic = face_route.irecv(aliens.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto recv_boundary  = face_route.irecv(aliens.faces.boundary, MpiTag::BOUNDARY);
    auto recv_normal    = face_route.irecv(aliens.faces.normal,   MpiTag::NORMAL);
    auto recv_f_center  = face_route.irecv(aliens.faces.center,   MpiTag::FACE_CENTER);
    auto recv_area      = face_route.irecv(aliens.faces.area,     MpiTag::AREA);
    auto recv_area_alt  = face_route.irecv(aliens.faces.area_alt, MpiTag::AREA_ALT);
    auto recv_face_vs   = face_route.irecv(aliens.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    auto recv_vertices  = node_route.irecv(aliens.verts, MpiTag::VERTICES);

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

    // ========================================================================

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

    /*
    static int counter = 0;
    static PvdFile pvd("alien");
    pvd.variables = {"rank", "index", "faces2D"};
    pvd.save(aliens, counter);
    ++counter;
     */
}

#endif

} // namespace zephyr::mesh