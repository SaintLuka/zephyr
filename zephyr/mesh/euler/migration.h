#pragma once

#include <vector>

#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/io/pvd_file.h>

using zephyr::io::PvdFile;

namespace zephyr::mesh {

#ifdef ZEPHYR_MPI

class Migration {
public:
    /// @brief Очистить буфферы
    void clear();

    /// @brief Привести буфферы к актуальным размерам
    void shrink_to_fit();

    /// @brief Основная функция перераспределения ячеек
    /// @param tourists Актуальные слои для пересылок
    /// @param locals Локальные ячейки
    /// @param aliens Слой удаленных ячеек
    /// @param vars Список переменных типа Storable<T>, которые переносятся
    /// при перераспределении ячеек.
    /// @details Входной параметр tourists требуется для связывании соседних
    /// ячеек при построении глобальной индексации примитивов.
    /// На данный момент не занимается построением обменного слоя, после
    /// миграции следует вызвать build_aliens().
    template <typename... Args>
    void migrate(
            Tourism & tourists,
            AmrCells& locals,
            AmrCells& aliens,
            Args&&... vars);

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

    /// @brief Вспомогательный буффер
    AmrCells migrants;

    Router cell_route;
    Router face_route;
    Router node_route;
};

template <typename Src, typename Dst, size_t... Is>
void data_assign_impl(Src& src, Dst& dst, index_t from, index_t to, std::index_sequence<Is...>) {
    ((std::get<Is>(dst)[to] = std::get<Is>(src)[from]), ...);
}

template <typename... Vars>
void Migration::fill_migrants(AmrCells& locals,
        const std::tuple<Vars...>& loc_vars,
        const std::tuple<Vars...>& mig_vars
        ) {
    // Для сортировки в migrants по ранку за О(n).
    // cell_offsets[i] показывает с какого индекса в migrants ставить i-ранковую ячейку.

    // Заполняем migrants, сортируем по rank
    migrants.resize(
            locals.n_cells(),
            locals.n_faces(),
            locals.n_nodes());

    // Сортировка migrants по rank, получается нормальный AmrCells
    // стартовые индексы (?)
    auto cell_index = cell_route.send_offset();
    auto face_index = face_route.send_offset();
    auto node_index = node_route.send_offset();

    // Кортеж указателей на буферы данных
    auto data_src = locals  .data.data_tuple(loc_vars);
    auto data_dst = migrants.data.data_tuple(mig_vars);

    for (index_t ic = 0; ic < migrants.size(); ++ic) {
        int r = locals.rank[ic];

        index_t jc = cell_index[r];

        locals.copy_geom(ic, migrants, jc, face_index[r], node_index[r]);
        /*
        for (int iface: migrants.faces_range(jc)) {
            if (migrants.faces.is_undefined(iface)) {
                continue;
            }
        }
         */

        // Копирование данных
        //data_dst[jc] = data_src[ic];
        data_assign_impl(data_src, data_dst, ic, jc, std::index_sequence_for<Vars...>{});

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

template <typename... Args>
void Migration::migrate(Tourism& tourists, AmrCells& locals, AmrCells& aliens, Args&&... vars) {
    // Все дополнительные переменные имеют тип Storable<T>
    utils::assert_all_storable<Args...>();

    // Кортежи переменных std::tuple{Storable<T1>, Storable<T2>, ...}
    // loc_vars - переменные в locals
    // mig_vars - переменные в migrants
    auto loc_vars = std::tuple{std::forward<Args>(vars)...};
    auto mig_vars = migrants.data.add_replace(std::forward<Args>(vars)...);

    // Составим полные матрицы пересылок (совпадают на всех процессах).
    // То есть заполним m_cell_route, m_face_route и m_node_route.
    fill_router(locals);

    // Переиндексировать ячейки, отправить новые индексы и ранги в alien-ячейки
    reindexing(tourists, locals, aliens);

    // Перенести ячейки из locals в migrants в нужном порядке
    fill_migrants(locals, loc_vars, mig_vars);

    // Меняем размеры locals
    locals.resize(cell_route.recv_buffer_size(),
                  face_route.recv_buffer_size(),
                  node_route.recv_buffer_size());

    // ========================== Пересылка ячеек =============================

    // Указатели на буферы данных
    auto data_src = migrants.data.data_tuple(mig_vars);
    auto data_dst = locals  .data.data_tuple(loc_vars);

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в face_begin и node_begin количество элементов на ячейку
    for (index_t ic = 0; ic < migrants.size(); ++ic) {
        migrants.face_begin[ic] = migrants.face_begin[ic + 1] - migrants.face_begin[ic];
        migrants.node_begin[ic] = migrants.node_begin[ic + 1] - migrants.node_begin[ic];
    }

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    auto send_rank     = cell_route.isend(migrants.rank,  MpiTag::RANK);
    auto send_next     = cell_route.isend(migrants.next,  MpiTag::NEXT);
    auto send_index    = cell_route.isend(migrants.index, MpiTag::INDEX);
    auto send_flag     = cell_route.isend(migrants.flag,  MpiTag::FLAG);
    auto send_level    = cell_route.isend(migrants.level, MpiTag::LEVEL);
    auto send_b_idx    = cell_route.isend(migrants.b_idx, MpiTag::B_IDX);
    auto send_z_idx    = cell_route.isend(migrants.z_idx, MpiTag::Z_IDX);
    auto send_center   = cell_route.isend(migrants.center, MpiTag::CENTER);
    auto send_volume   = cell_route.isend(migrants.volume, MpiTag::VOLUME);
    auto send_volume_a = cell_route.isend(migrants.volume_alt, MpiTag::VOLUME_ALT);
    auto send_face_beg = cell_route.isend(migrants.face_begin, MpiTag::FACE_BEG);
    auto send_node_beg = cell_route.isend(migrants.node_begin, MpiTag::NODE_BEG);

    // Отправить данные ячеек
    auto send_data     = [&]<size_t... Is>(std::index_sequence<Is...> ) {
        return std::tuple{ cell_route.isend(std::get<Is>(data_src), MpiTag(12345 + Is))...  };
    }(std::index_sequence_for<Args...>());

    // Отправить данные граней
    auto send_adj_rank  = face_route.isend(migrants.faces.adjacent.rank,  MpiTag::ADJ_RANK);
    auto send_adj_index = face_route.isend(migrants.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto send_adj_alien = face_route.isend(migrants.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto send_adj_basic = face_route.isend(migrants.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto send_boundary  = face_route.isend(migrants.faces.boundary, MpiTag::BOUNDARY);
    auto send_normal    = face_route.isend(migrants.faces.normal,   MpiTag::NORMAL);
    auto send_f_center  = face_route.isend(migrants.faces.center,   MpiTag::FACE_CENTER);
    auto send_area      = face_route.isend(migrants.faces.area,     MpiTag::AREA);
    auto send_area_alt  = face_route.isend(migrants.faces.area_alt, MpiTag::AREA_ALT);
    auto send_face_vs   = face_route.isend(migrants.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    auto send_vertices  = node_route.isend(migrants.verts, MpiTag::VERTICES);

    // ============================= IRECV ====================================

    // Получить данные ячеек
    auto recv_rank     = cell_route.irecv(locals.rank, MpiTag::RANK);
    auto recv_next     = cell_route.irecv(locals.next, MpiTag::NEXT);
    auto recv_index    = cell_route.irecv(locals.index, MpiTag::INDEX);
    auto recv_flag     = cell_route.irecv(locals.flag, MpiTag::FLAG);
    auto recv_level    = cell_route.irecv(locals.level, MpiTag::LEVEL);
    auto recv_b_idx    = cell_route.irecv(locals.b_idx, MpiTag::B_IDX);
    auto recv_z_idx    = cell_route.irecv(locals.z_idx, MpiTag::Z_IDX);
    auto recv_center   = cell_route.irecv(locals.center, MpiTag::CENTER);
    auto recv_volume   = cell_route.irecv(locals.volume, MpiTag::VOLUME);
    auto recv_volume_a = cell_route.irecv(locals.volume_alt, MpiTag::VOLUME_ALT);

    // Получить данные ячеек
    auto recv_data     = [&]<size_t... Is>(std::index_sequence<Is...> ) {
        return std::tuple{ cell_route.irecv(std::get<Is>(data_dst), MpiTag(12345 + Is))...  };
    }(std::index_sequence_for<Args...>());

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    auto recv_face_beg = cell_route.irecv(locals.face_begin.data() + 1, MpiTag::FACE_BEG);
    auto recv_node_beg = cell_route.irecv(locals.node_begin.data() + 1, MpiTag::NODE_BEG);

    // Получить данные граней
    auto recv_adj_rank  = face_route.irecv(locals.faces.adjacent.rank, MpiTag::ADJ_RANK);
    auto recv_adj_index = face_route.irecv(locals.faces.adjacent.index, MpiTag::ADJ_INDEX);
    auto recv_adj_alien = face_route.irecv(locals.faces.adjacent.alien, MpiTag::ADJ_ALIEN);
    auto recv_adj_basic = face_route.irecv(locals.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    auto recv_boundary  = face_route.irecv(locals.faces.boundary, MpiTag::BOUNDARY);
    auto recv_normal    = face_route.irecv(locals.faces.normal, MpiTag::NORMAL);
    auto recv_f_center  = face_route.irecv(locals.faces.center, MpiTag::FACE_CENTER);
    auto recv_area      = face_route.irecv(locals.faces.area, MpiTag::AREA);
    auto recv_area_alt  = face_route.irecv(locals.faces.area_alt, MpiTag::AREA_ALT);
    auto recv_face_vs   = face_route.irecv(locals.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    auto recv_vertices  = node_route.irecv(locals.verts, MpiTag::VERTICES);

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

    // Завершить отправку данных ячеек
    [&]<size_t... Is>(std::index_sequence<Is...> ) {
        ( std::get<Is>(send_data).wait(), ... );
    }(std::index_sequence_for<Args...>());

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

    // Завершить получение данных ячеек
    [&]<size_t... Is>(std::index_sequence<Is...> ) {
        ( std::get<Is>(recv_data).wait(), ... );
    }(std::index_sequence_for<Args...>());

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
}

#endif // ZEPHYR_MPI

} // namespace zephyr::mesh