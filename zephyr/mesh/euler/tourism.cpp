#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/mesh/euler/router.h>

#ifdef ZEPHYR_MPI

namespace zephyr::mesh {

using zephyr::utils::mpi;

void Tourism::shrink_to_fit() {
    m_unique_border_indices.shrink_to_fit();
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
    m_unique_border_indices.reserve(m_cell_route.send_buffer_size());

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
                    m_unique_border_indices.push_back(ic);
                }
            }
        }
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

    // Копируем геометрию в border-слой
    index_t face_idx = 0;
    index_t node_idx = 0;
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        locals.copy_geom(m_border_indices[ic], m_border, ic, face_idx, node_idx);

        face_idx += locals.max_face_count(m_border_indices[ic]);
        node_idx += locals.max_node_count(m_border_indices[ic]);
    }
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

    // Копируем геометрию в border-слой
    index_t face_idx = 0;
    index_t node_idx = 0;
    for (size_t ic = 0; ic < m_border_indices.size(); ++ic) {
        locals.copy_geom_basic(m_border_indices[ic], m_border, ic, face_idx, node_idx);

        face_idx += locals.max_face_count(m_border_indices[ic]);
        node_idx += locals.max_node_count(m_border_indices[ic]);
    }
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

} // namespace zephyr::mesh

#endif