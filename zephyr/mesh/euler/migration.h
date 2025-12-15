#pragma once

#include <vector>

#include <zephyr/mesh/euler/tourism.h>

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI

class Migration {
public:
    /// @brief Очистить буферы
    void clear();

    /// @brief Привести буферы к актуальным размерам
    void shrink_to_fit();

    /// @brief Основная функция перераспределения ячеек
    /// @param tourists Актуальные слои для пересылок
    /// @param locals Локальные ячейки
    /// @param vars Список переменных типа Storable<T>, которые переносятся
    /// при перераспределении ячеек.
    /// @details Входной параметр tourists требуется для связывания соседних
    /// ячеек при построении глобальной индексации примитивов.
    /// На данный момент не занимается построением обменного слоя, после
    /// миграции следует вызвать build_aliens().
    template <typename... Args>
    void migrate(Tourism& tourists, AmrCells& locals, Args&&... vars);

protected:
    /// @brief Составить матрицу пересылок перед миграцей.
    /// Проверяет новые ранги ячеек и подсчитывает число пересылок ячеек,
    /// граней и вершин. Заполняет cell_route, face_route, node_route.
    /// Требует коллективной MPI-операции, по типу all-to-all.
    void fill_router(AmrCells& locals);

    // Новая индексация ячеек (какая будет после миграции), пересылка
    // и получение новых index и rank в alien-слой.
    void reindexing(Tourism& tourism,
                    AmrCells& locals,
                    AmrCells& aliens);

    // Копировать геометрию и поля данных в хранилище migrants
    // Все аргументы Vars должны иметь тип Storable<T>
    template <typename... Vars>
    void fill_migrants(AmrCells& locals,
            const std::tuple<Vars...>& loc_vars,
            const std::tuple<Vars...>& mig_vars);

    /// @brief Вспомогательный буфер
    AmrCells migrants;

    Router m_cell_router;
    Router m_face_router;
    Router m_node_router;
};

template <typename... Vars>
void Migration::fill_migrants(AmrCells& locals,
        const std::tuple<Vars...>& loc_vars,
        const std::tuple<Vars...>& mig_vars
        ) {
    using utils::Buffer;

    // Все переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Для сортировки в migrants по ранку за О(n).
    // cell_offsets[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.

    // Заполняем migrants, сортируем по rank
    migrants.resize(
            locals.n_cells(),
            locals.n_faces(),
            locals.n_nodes());

    // Сортировка migrants по rank, получается нормальный AmrCells
    // стартовые индексы (?)
    auto cell_index = m_cell_router.send_offset();
    auto face_index = m_face_router.send_offset();
    auto node_index = m_node_router.send_offset();

    // Кортеж указателей на буферы данных
    std::array<Buffer*, n_vars> data_src = locals  .data[loc_vars];
    std::array<Buffer*, n_vars> data_dst = migrants.data[mig_vars];

    for (index_t ic = 0; ic < migrants.size(); ++ic) {
        int r = locals.rank[ic];

        index_t jc = cell_index[r];

        locals.copy_geom(ic, migrants, jc, face_index[r], node_index[r]);

        // Копирование данных
        //data_dst[jc] = data_src[ic];
        for (int v = 0; v < n_vars; ++v) {
            data_src[v]->copy_data(ic, *data_dst[v], jc);
        }

        cell_index[r] += 1;
        face_index[r] += locals.max_face_count(ic);
        node_index[r] += locals.max_node_count(ic);
    }

    /*
    if (mpi::master()) {
        std::cout << "locals.fb: " << locals.face_begin << "\n";
        std::cout << "locals.nb: " << locals.node_begin << "\n";

        std::cout << "migrants.rank: " << migrants.rank << "\n";
        std::cout << "migrants.b_idx: " << migrants.b_idx << "\n";
        std::cout << "migrants.index: " << migrants.index << "\n";
        std::cout << "migrants.fb: " << migrants.face_begin << "\n";
        std::cout << "migrants.nb: " << migrants.node_begin << "\n";
    }
    */

    //PvdFile pvd("migrants");
    //pvd.variables = {"rank", "index", "faces2D"};
    //pvd.save(migrants, 0.0);
}

template <typename... Vars>
void Migration::migrate(Tourism& tourists, AmrCells& locals, Vars&&... vars) {
    using utils::Buffer;
    AmrCells& aliens = tourists.aliens();

    // Все дополнительные переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Кортежи переменных std::tuple{Storable<T1>, Storable<T2>, ...}
    // loc_vars - переменные в locals
    // mig_vars - переменные в migrants
    auto loc_vars = std::tuple{std::forward<Vars>(vars)...};
    auto mig_vars = migrants.data.add_replace(locals.data, loc_vars);

    // Составим полные матрицы пересылок (совпадают на всех процессах).
    // То есть заполним m_cell_route, m_face_route и m_node_route.
    fill_router(locals);

    // Переиндексировать ячейки, отправить новые индексы и ранги в alien-ячейки
    reindexing(tourists, locals, aliens);

    // Перенести ячейки из locals в migrants в нужном порядке
    fill_migrants(locals, loc_vars, mig_vars);

    // Меняем размеры locals
    locals.resize(m_cell_router.recv_buffer_size(),
                  m_face_router.recv_buffer_size(),
                  m_node_router.recv_buffer_size());

    // ========================== Пересылка ячеек =============================

    // Указатели на буферы данных
    std::array<Buffer*, n_vars> data_src = migrants.data[mig_vars];
    std::array<Buffer*, n_vars> data_dst = locals.data[loc_vars];

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в face_begin и node_begin количество элементов на ячейку
    for (index_t ic = 0; ic < migrants.size(); ++ic) {
        migrants.face_begin[ic] = migrants.face_begin[ic + 1] - migrants.face_begin[ic];
        migrants.node_begin[ic] = migrants.node_begin[ic + 1] - migrants.node_begin[ic];
    }

    // ============================= ISEND ====================================

    // Отправить геометрию ячеек
    RequestsList cells_send; cells_send.reserve(16);
    cells_send += m_cell_router.isend(migrants.rank,  MpiTag::RANK);
    cells_send += m_cell_router.isend(migrants.next,  MpiTag::NEXT);
    cells_send += m_cell_router.isend(migrants.index, MpiTag::INDEX);
    cells_send += m_cell_router.isend(migrants.flag,  MpiTag::FLAG);
    cells_send += m_cell_router.isend(migrants.level, MpiTag::LEVEL);
    cells_send += m_cell_router.isend(migrants.b_idx, MpiTag::B_IDX);
    cells_send += m_cell_router.isend(migrants.z_idx, MpiTag::Z_IDX);
    cells_send += m_cell_router.isend(migrants.center, MpiTag::CENTER);
    cells_send += m_cell_router.isend(migrants.volume, MpiTag::VOLUME);
    cells_send += m_cell_router.isend(migrants.volume_alt, MpiTag::VOLUME_ALT);
    cells_send += m_cell_router.isend(migrants.face_begin, MpiTag::FACE_BEG);
    cells_send += m_cell_router.isend(migrants.node_begin, MpiTag::NODE_BEG);

    // Отправить данные ячеек
    RequestsList data_send; data_send.reserve(n_vars);
    for (int v = 0; v < n_vars; ++v) {
        data_send += m_cell_router.isend(*data_src[v], static_cast<MpiTag>(12345 + v));
    }

    // Отправить данные граней
    RequestsList faces_send; faces_send.reserve(16);
    faces_send += m_face_router.isend(migrants.faces.adjacent.rank,  MpiTag::ADJ_RANK);
    faces_send += m_face_router.isend(migrants.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_send += m_face_router.isend(migrants.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    faces_send += m_face_router.isend(migrants.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_send += m_face_router.isend(migrants.faces.boundary, MpiTag::BOUNDARY);
    faces_send += m_face_router.isend(migrants.faces.normal,   MpiTag::NORMAL);
    faces_send += m_face_router.isend(migrants.faces.center,   MpiTag::FACE_CENTER);
    faces_send += m_face_router.isend(migrants.faces.area,     MpiTag::AREA);
    faces_send += m_face_router.isend(migrants.faces.area_alt, MpiTag::AREA_ALT);
    faces_send += m_face_router.isend(migrants.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    auto nodes_send = m_node_router.isend(migrants.verts, MpiTag::VERTICES);

    // ============================= IRECV ====================================

    // Получить геометрию ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += m_cell_router.irecv(locals.rank, MpiTag::RANK);
    cells_recv += m_cell_router.irecv(locals.next, MpiTag::NEXT);
    cells_recv += m_cell_router.irecv(locals.index, MpiTag::INDEX);
    cells_recv += m_cell_router.irecv(locals.flag, MpiTag::FLAG);
    cells_recv += m_cell_router.irecv(locals.level, MpiTag::LEVEL);
    cells_recv += m_cell_router.irecv(locals.b_idx, MpiTag::B_IDX);
    cells_recv += m_cell_router.irecv(locals.z_idx, MpiTag::Z_IDX);
    cells_recv += m_cell_router.irecv(locals.center, MpiTag::CENTER);
    cells_recv += m_cell_router.irecv(locals.volume, MpiTag::VOLUME);
    cells_recv += m_cell_router.irecv(locals.volume_alt, MpiTag::VOLUME_ALT);

    // Получить данные ячеек
    RequestsList data_recv; data_recv.reserve(n_vars);
    for (int v = 0; v < n_vars; ++v) {
        data_recv += m_cell_router.irecv(*data_dst[v], static_cast<MpiTag>(12345 + v));
    };

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += m_cell_router.irecv(locals.face_begin.data() + 1, MpiTag::FACE_BEG);
    faces_recv += m_cell_router.irecv(locals.node_begin.data() + 1, MpiTag::NODE_BEG);

    // Получить данные граней
    faces_recv += m_face_router.irecv(locals.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += m_face_router.irecv(locals.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += m_face_router.irecv(locals.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    faces_recv += m_face_router.irecv(locals.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += m_face_router.irecv(locals.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += m_face_router.irecv(locals.faces.normal, MpiTag::NORMAL);
    faces_recv += m_face_router.irecv(locals.faces.center, MpiTag::FACE_CENTER);
    faces_recv += m_face_router.irecv(locals.faces.area, MpiTag::AREA);
    faces_recv += m_face_router.irecv(locals.faces.area_alt, MpiTag::AREA_ALT);
    faces_recv += m_face_router.irecv(locals.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    auto nodes_recv = m_node_router.irecv(locals.verts, MpiTag::VERTICES);

    // =========================== WAIT ISEND =================================

    cells_send.wait();  // Завершить отправку ячеек
    data_send.wait();   // Завершить отправку данных ячеек
    faces_send.wait();  // Завершить отправку граней
    nodes_send.wait();  // Завершить отправку вершин

    // =========================== WAIT IRECV =================================

    cells_recv.wait();  // Завершить получение ячеек
    data_recv.wait();   // Завершить получение данных ячеек
    faces_recv.wait();  // Завершить получение граней
    nodes_recv.wait();  // Завершить получение вершин

    // Восстановить индексацию граней
    locals.face_begin[0] = 0;
    locals.node_begin[0] = 0;
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        locals.face_begin[ic + 1] += locals.face_begin[ic];
        locals.node_begin[ic + 1] += locals.node_begin[ic];
    }

    // ========================================================================

    /*
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "  fb: " << locals.face_begin << "\n";
        std::cout << "  nb: " << locals.node_begin << "\n";
    });
    */
    /*
    static int counter = 0;
    static PvdFile pvd("migrants");
    pvd.variables = {"rank", "index", "faces2D"};
    pvd.save(locals, counter);
    ++counter;

    mpi::barrier();
    mpi::for_each([]() {
        std::cout << "SUCCESS MIGRATION " << mpi::rank() << "\n";
    });
     */

    // Построить обменные слои
    tourists.update(locals);
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh