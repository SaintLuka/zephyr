#pragma once

#include <vector>

#include <zephyr/mesh/raw/tourism.h>

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI

class Migration {
public:
    /// @brief Инициализирует буфер с теми же опциями
    void init_types(const RawCells& cells);

    /// @brief Очистить буферы
    void clear();

    /// @brief Привести буферы к актуальным размерам
    void shrink_to_fit();

    /// @brief Основная функция перераспределения ячеек
    /// @param tourism Актуальные слои для пересылок
    /// @param cells Локальные ячейки
    /// @param vars Список переменных типа Storable<T>, которые переносятся
    /// при перераспределении ячеек.
    /// @details Входной параметр tourists требуется для связывания соседних
    /// ячеек при построении глобальной индексации примитивов.
    /// На данный момент не занимается построением обменного слоя, после
    /// миграции следует вызвать build_ghosts().
    template <typename... Args>
    void migrate(Tourism& tourism, RawCells& cells, RawNodes& nodes, Args&&... vars);

protected:
    /// @brief Выставить новые ранги узлов исходя из рангов ячеек.
    /// Владельцем узла считается процесс с минимальным рангом среди
    /// инцидентных ячеек.
    static void setup_node_ranks(RawNodes& nodes, const RawCells& locals, const RawCells& ghosts);

    /// @brief Составить матрицу пересылок перед миграцей.
    /// Проверяет новые ранги ячеек и подсчитывает число пересылок ячеек,
    /// граней и вершин. Заполняет cell_router, face_router, vert_route.
    /// Требует коллективной MPI-операции, по типу all-to-all.
    void fill_cell_routers(const RawCells& cells);

    void fill_node_routers(const RawNodes& nodes, bool unique_nodes);

    // Новая индексация ячеек (какая будет после миграции), пересылка
    // и получение новых index и rank в ghost-слой.
    void cells_reindexing(Tourism& tourism, RawCells& cells) const;

    // Новая индексация узлов (какая будет после миграции), пересылка
    // и получение новых index и rank в ghost-слой.
    void nodes_reindexing(Tourism& tourism, RawNodes& nodes) const;

    // Обновить face.adjacent.index в соответствии с индексами в массивах
    static void update_face_adjacent(RawCells& locals, const RawCells& ghosts);

    static void update_cell_verts(RawVerts& verts, const RawNodes& locals, const RawNodes& ghosts);

    // Копировать геометрию и поля данных в хранилище migrants
    // Все аргументы Vars должны иметь тип Storable<T>
    template <typename... Vars>
    void fill_migrant_cells(RawCells& cells,
        const std::tuple<Vars...>& loc_vars,
        const std::tuple<Vars...>& mig_vars);

    // Копировать геометрию и поля данных в хранилище migrants
    // Все аргументы Vars должны иметь тип Storable<T>
    template <typename... Vars>
    void fill_migrant_nodes(RawNodes& nodes,
        const std::tuple<Vars...>& loc_vars,
        const std::tuple<Vars...>& mig_vars);

    /// @brief Вспомогательный буфер ячеек
    RawCells cell_buffer_;

    // Маршрутизаторы для отправки примитивов из cells
    Router cell_router_;
    Router face_router_;
    Router vert_router_;

    /// @brief Вспомогательный буфер узлов
    RawNodes node_buffer_;

    // Маршрутизаторы для отправки примитивов из cells
    Router node_router_;
    Router inct_router_;
};

template <typename... Vars>
void Migration::fill_migrant_cells(
    RawCells& cells,
    const std::tuple<Vars...>& loc_vars,
    const std::tuple<Vars...>& mig_vars) {
    using utils::Buffer;

    // Все переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Для сортировки в migrants по ранку за О(n).
    // cell_offsets[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.

    // Заполняем migrants, сортируем по rank
    z_assert(cells.has_nodes() == cell_buffer_.has_nodes(), "Migration: buffer and cells options mismatch");
    cell_buffer_.resize(
            cells.n_cells(),
            cells.n_faces(),
            cells.n_verts());

    // Сортировка migrants по rank, получается нормальный RawCells
    // стартовые индексы (?)
    auto cell_index = cell_router_.send_offset();
    auto face_index = face_router_.send_offset();
    auto vert_index = vert_router_.send_offset();

    // Кортеж указателей на буферы данных
    std::array<Buffer*, n_vars> data_src = cells.data[loc_vars];
    std::array<Buffer*, n_vars> data_dst = cell_buffer_.data[mig_vars];

    for (index_t ic = 0; ic < cell_buffer_.n_cells(); ++ic) {
        int r = cells.rank[ic];

        index_t jc = cell_index[r];

        cells.copy_geom(ic, cell_buffer_, jc, face_index[r], vert_index[r]);

        // Копирование данных
        //data_dst[jc] = data_src[ic];
        for (int v = 0; v < n_vars; ++v) {
            data_src[v]->copy_data(ic, *data_dst[v], jc);
        }

        cell_index[r] += 1;
        face_index[r] += cells.faces.max_count(ic);
        vert_index[r] += cells.verts.max_count(ic);
    }

    /*
    if (mpi::master()) {
        std::cout << "cells.fb: " << cells.faces.offsets << "\n";
        std::cout << "cells.nb: " << cells.verts.offsets << "\n";

        std::cout << "migrants.rank: " << migrants.rank << "\n";
        std::cout << "migrants.b_idx: " << migrants.b_idx << "\n";
        std::cout << "migrants.index: " << migrants.index << "\n";
        std::cout << "migrants.fb: " << migrants.faces.offsets << "\n";
        std::cout << "migrants.nb: " << migrants.verts.offsets << "\n";
    }
    */

    //PvdFile pvd("migrants");
    //pvd.variables = {"rank", "index", "faces2D"};
    //pvd.save(migrants, 0.0);
}

template <typename... Vars>
void Migration::fill_migrant_nodes(
    RawNodes& nodes,
    const std::tuple<Vars...>& loc_vars,
    const std::tuple<Vars...>& mig_vars) {
    using utils::Buffer;

    // Все переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Для сортировки в migrants по ранку за О(n).
    // cell_offsets[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.

    // Заполняем migrants, сортируем по rank
    node_buffer_.resize(
            nodes.n_nodes(),
            nodes.n_incident());

    // Сортировка migrants по rank, получается нормальный RawCells
    // стартовые индексы (?)
    auto node_index = node_router_.send_offset();
    auto inc_offset = inct_router_.send_offset();

    for (index_t in = 0; in < node_buffer_.n_nodes(); ++in) {
        int r = nodes.rank[in];
        index_t jn = node_index[r];

        nodes.copy_geom(in, node_buffer_, jn, inc_offset[r]);

        node_index[r] += 1;
        inc_offset[r] += nodes.incident.max_count(in);
    }
}

template <typename... Vars>
void Migration::migrate(Tourism& tourism, RawCells& cells, RawNodes& nodes, Vars&&... vars) {
    using utils::Buffer;
    RawCells& ghost_cells = tourism.ghost_cells();
    RawNodes& ghost_nodes = tourism.ghost_nodes();

    // Все дополнительные переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Кортежи переменных std::tuple{Storable<T1>, Storable<T2>, ...}
    // loc_vars - переменные в cells
    // mig_vars - переменные в migrants
    auto loc_vars = std::tuple{std::forward<Vars>(vars)...};
    auto mig_vars = cell_buffer_.data.add_replace(cells.data, loc_vars);

    // Составим полные матрицы пересылок (совпадают на всех процессах).
    // То есть заполним cell_router, face_router и vert_router.
    fill_cell_routers(cells);

    // Переиндексировать ячейки, отправить новые индексы и ранги в ghost-ячейки
    cells_reindexing(tourism, cells);

    // Обновить значения faces.adjacent.index
    update_face_adjacent(cells, ghost_cells);

    // Выставить ранги узлов, обновить nodes.incident.index
    setup_node_ranks(nodes, cells, ghost_cells);

    // Составим полные матрицы пересылок (совпадают на всех процессах).
    // То есть заполним node_router, inct_router.
    fill_node_routers(nodes, cells.has_nodes());

    if (cells.has_nodes()) {
        // Переиндексировать узлы, отправить новые индексы и ранги в ghost-узлы
        nodes_reindexing(tourism, nodes);

        // Обновить значения cells.verts.rank, cells.verts.index
        update_cell_verts(cells.verts, nodes, ghost_nodes);
    }

    // Перенести ячейки из cells в migrants в нужном порядке
    fill_migrant_cells(cells, loc_vars, mig_vars);

    if (cells.has_nodes()) {
        fill_migrant_nodes(nodes, {}, {});
    }

    // Меняем размеры cells
    cells.resize(cell_router_.recv_buffer_size(),
                 face_router_.recv_buffer_size(),
                 vert_router_.recv_buffer_size());

    nodes.resize(node_router_.recv_buffer_size(),
                 inct_router_.recv_buffer_size());

    // ========================== Пересылка ячеек =============================

    // Указатели на буферы данных
    std::array<Buffer*, n_vars> data_src = cell_buffer_.data[mig_vars];
    std::array<Buffer*, n_vars> data_dst = cells.data[loc_vars];

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в faces.offsets и verts.offsets количество элементов на ячейку
    for (index_t ic = 0; ic < cell_buffer_.n_cells(); ++ic) {
        cell_buffer_.faces.offsets[ic] = cell_buffer_.faces.offsets[ic + 1] - cell_buffer_.faces.offsets[ic];
        cell_buffer_.verts.offsets[ic] = cell_buffer_.verts.offsets[ic + 1] - cell_buffer_.verts.offsets[ic];
    }

    if (cells.has_nodes()) {
        for (index_t in = 0; in < node_buffer_.n_nodes(); ++in) {
            node_buffer_.incident.offsets[in] = node_buffer_.incident.offsets[in + 1] - node_buffer_.incident.offsets[in];
        }
    }

    // ============================= ISEND ====================================

    // Отправить геометрию ячеек
    RequestsList cells_send; cells_send.reserve(16);
    cells_send += cell_router_.isend(cell_buffer_.rank,  MpiTag::RANK);
    cells_send += cell_router_.isend(cell_buffer_.next,  MpiTag::NEXT);
    cells_send += cell_router_.isend(cell_buffer_.index, MpiTag::INDEX);
    cells_send += cell_router_.isend(cell_buffer_.flag,  MpiTag::FLAG);
    cells_send += cell_router_.isend(cell_buffer_.level, MpiTag::LEVEL);
    cells_send += cell_router_.isend(cell_buffer_.b_idx, MpiTag::B_IDX);
    cells_send += cell_router_.isend(cell_buffer_.z_idx, MpiTag::Z_IDX);
    cells_send += cell_router_.isend(cell_buffer_.center, MpiTag::CENTER);
    cells_send += cell_router_.isend(cell_buffer_.volume, MpiTag::VOLUME);
    cells_send += cell_router_.isend(cell_buffer_.volume_alt, MpiTag::VOLUME_ALT);
    cells_send += cell_router_.isend(cell_buffer_.faces.offsets, MpiTag::FACE_BEG);
    cells_send += cell_router_.isend(cell_buffer_.verts.offsets, MpiTag::VERT_BEG);

    // Отправить данные ячеек
    RequestsList data_send; data_send.reserve(n_vars);
    for (int v = 0; v < n_vars; ++v) {
        data_send += cell_router_.isend(*data_src[v], static_cast<MpiTag>(12345 + v));
    }

    // Отправить данные граней
    RequestsList faces_send; faces_send.reserve(16);
    faces_send += face_router_.isend(cell_buffer_.faces.adjacent.rank,  MpiTag::ADJ_RANK);
    faces_send += face_router_.isend(cell_buffer_.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_send += face_router_.isend(cell_buffer_.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_send += face_router_.isend(cell_buffer_.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_send += face_router_.isend(cell_buffer_.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_send += face_router_.isend(cell_buffer_.faces.boundary, MpiTag::BOUNDARY);
    faces_send += face_router_.isend(cell_buffer_.faces.normal,   MpiTag::NORMAL);
    faces_send += face_router_.isend(cell_buffer_.faces.center,   MpiTag::FACE_CENTER);
    faces_send += face_router_.isend(cell_buffer_.faces.area,     MpiTag::AREA);
    faces_send += face_router_.isend(cell_buffer_.faces.area_alt, MpiTag::AREA_ALT);
    faces_send += face_router_.isend(cell_buffer_.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    RequestsList verts_send; verts_send.reserve(3);
    verts_send += vert_router_.isend(cell_buffer_.verts.coord, MpiTag::VERT_COORD);
    if (cell_buffer_.verts.has_nodes()) {
        verts_send += vert_router_.isend(cell_buffer_.verts.rank,  MpiTag::VERT_RANK);
        verts_send += vert_router_.isend(cell_buffer_.verts.index, MpiTag::VERT_INDEX);
        verts_send += vert_router_.isend(cell_buffer_.verts.ghost, MpiTag::VERT_GHOST);
    }

    // Отправить данные узлов
    RequestsList nodes_send;
    RequestsList incident_send;
    if (cells.has_nodes()) {
        nodes_send.reserve(5);
        nodes_send += node_router_.isend(node_buffer_.rank,  MpiTag::NODE_RANK);
        nodes_send += node_router_.isend(node_buffer_.next,  MpiTag::NODE_NEXT);
        nodes_send += node_router_.isend(node_buffer_.index, MpiTag::NODE_INDEX);
        nodes_send += node_router_.isend(node_buffer_.coord, MpiTag::NODE_COORD);
        nodes_send += node_router_.isend(node_buffer_.incident.offsets, MpiTag::INCT_BEG);

        // Отправить списки инцидентных ячеек
        incident_send.reserve(4);
        incident_send += inct_router_.isend(node_buffer_.incident.role,  MpiTag::INCT_ROLE);
        incident_send += inct_router_.isend(node_buffer_.incident.rank,  MpiTag::INCT_RANK);
        incident_send += inct_router_.isend(node_buffer_.incident.index, MpiTag::INCT_INDEX);
        incident_send += inct_router_.isend(node_buffer_.incident.ghost, MpiTag::INCT_GHOST);
    }

    // ============================= IRECV ====================================

    // Получить геометрию ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += cell_router_.irecv(cells.rank, MpiTag::RANK);
    cells_recv += cell_router_.irecv(cells.next, MpiTag::NEXT);
    cells_recv += cell_router_.irecv(cells.index, MpiTag::INDEX);
    cells_recv += cell_router_.irecv(cells.flag, MpiTag::FLAG);
    cells_recv += cell_router_.irecv(cells.level, MpiTag::LEVEL);
    cells_recv += cell_router_.irecv(cells.b_idx, MpiTag::B_IDX);
    cells_recv += cell_router_.irecv(cells.z_idx, MpiTag::Z_IDX);
    cells_recv += cell_router_.irecv(cells.center, MpiTag::CENTER);
    cells_recv += cell_router_.irecv(cells.volume, MpiTag::VOLUME);
    cells_recv += cell_router_.irecv(cells.volume_alt, MpiTag::VOLUME_ALT);

    // Получить данные ячеек
    RequestsList data_recv; data_recv.reserve(n_vars);
    for (int v = 0; v < n_vars; ++v) {
        data_recv += cell_router_.irecv(*data_dst[v], static_cast<MpiTag>(12345 + v));
    };

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += cell_router_.irecv(cells.faces.offsets.data() + 1, MpiTag::FACE_BEG);
    faces_recv += cell_router_.irecv(cells.verts.offsets.data() + 1, MpiTag::VERT_BEG);

    // Получить данные граней
    faces_recv += face_router_.irecv(cells.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += face_router_.irecv(cells.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += face_router_.irecv(cells.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_recv += face_router_.irecv(cells.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += face_router_.irecv(cells.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_recv += face_router_.irecv(cells.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += face_router_.irecv(cells.faces.normal, MpiTag::NORMAL);
    faces_recv += face_router_.irecv(cells.faces.center, MpiTag::FACE_CENTER);
    faces_recv += face_router_.irecv(cells.faces.area, MpiTag::AREA);
    faces_recv += face_router_.irecv(cells.faces.area_alt, MpiTag::AREA_ALT);
    faces_recv += face_router_.irecv(cells.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    RequestsList verts_recv; verts_recv.reserve(3);
    verts_recv += vert_router_.irecv(cells.verts.coord, MpiTag::VERT_COORD);
    if (cell_buffer_.verts.has_nodes()) {
        verts_recv += vert_router_.irecv(cells.verts.rank,  MpiTag::VERT_RANK);
        verts_recv += vert_router_.irecv(cells.verts.index, MpiTag::VERT_INDEX);
        verts_recv += vert_router_.irecv(cells.verts.ghost, MpiTag::VERT_GHOST);
    }

    // Получить данные узлов
    RequestsList nodes_recv;
    RequestsList incident_recv;
    if (cells.has_nodes()) {
        nodes_recv.reserve(5);
        nodes_recv += node_router_.irecv(nodes.rank,  MpiTag::NODE_RANK);
        nodes_recv += node_router_.irecv(nodes.next,  MpiTag::NODE_NEXT);
        nodes_recv += node_router_.irecv(nodes.index, MpiTag::NODE_INDEX);
        nodes_recv += node_router_.irecv(nodes.coord, MpiTag::NODE_COORD);

        // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
        nodes_recv += node_router_.irecv(nodes.incident.offsets.data() + 1, MpiTag::INCT_BEG);

        // Получить данные инцидентных ячеек
        incident_recv.reserve(4);
        incident_recv += inct_router_.irecv(nodes.incident.role,  MpiTag::INCT_ROLE);
        incident_recv += inct_router_.irecv(nodes.incident.rank,  MpiTag::INCT_RANK);
        incident_recv += inct_router_.irecv(nodes.incident.index, MpiTag::INCT_INDEX);
        incident_recv += inct_router_.irecv(nodes.incident.ghost, MpiTag::INCT_GHOST);
    }

    // =========================== WAIT ISEND =================================

    cells_send.wait();  // Завершить отправку ячеек
    data_send.wait();   // Завершить отправку данных ячеек
    faces_send.wait();  // Завершить отправку граней
    verts_send.wait();  // Завершить отправку вершин
    if (cells.has_nodes()) {
        nodes_send.wait();
        incident_send.wait();
    }

    // =========================== WAIT IRECV =================================

    cells_recv.wait();  // Завершить получение ячеек
    data_recv.wait();   // Завершить получение данных ячеек
    faces_recv.wait();  // Завершить получение граней
    verts_recv.wait();  // Завершить получение вершин
    if (cells.has_nodes()) {
        nodes_recv.wait();
        incident_recv.wait();
    }

    // Восстановить индексацию граней
    cells.faces.offsets[0] = 0;
    cells.verts.offsets[0] = 0;
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        cells.faces.offsets[ic + 1] += cells.faces.offsets[ic];
        cells.verts.offsets[ic + 1] += cells.verts.offsets[ic];
    }

    if (cells.has_nodes()) {
        for (index_t in = 0; in < nodes.n_nodes(); ++in) {
            nodes.incident.offsets[in + 1] += nodes.incident.offsets[in];
        }
    }

    // ========================================================================

    /*
    {
        static int counter = 0;
        static io::PvdFile pvd("migrants");
        pvd.variables = {"rank", "index", "faces2D", "verts2D"};
        pvd.save(cells, counter);
        ++counter;
        save_markers(nodes, "nodes");
    }
    utils::mpi::barrier();
    utils::mpi::for_each([]() {
        std::cout << "SUCCESS MIGRATION " << utils::mpi::rank() << "\n";
    });
    */

    // Построить обменные слои
    tourism.update(cells, nodes);
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh