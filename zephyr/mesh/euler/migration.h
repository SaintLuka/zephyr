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
    /// миграции следует вызвать build_ghosts().
    template <typename... Args>
    void migrate(Tourism& tourists, AmrCells& locals, Args&&... vars);

protected:
    /// @brief Составить матрицу пересылок перед миграцей.
    /// Проверяет новые ранги ячеек и подсчитывает число пересылок ячеек,
    /// граней и вершин. Заполняет cell_router, face_router, vert_route.
    /// Требует коллективной MPI-операции, по типу all-to-all.
    void fill_router(AmrCells& locals);

    // Новая индексация ячеек (какая будет после миграции), пересылка
    // и получение новых index и rank в ghost-слой.
    void reindexing(Tourism& tourism,
                    AmrCells& locals,
                    AmrCells& ghosts);

    // Копировать геометрию и поля данных в хранилище migrants
    // Все аргументы Vars должны иметь тип Storable<T>
    template <typename... Vars>
    void fill_migrants(AmrCells& locals,
            const std::tuple<Vars...>& loc_vars,
            const std::tuple<Vars...>& mig_vars);

    /// @brief Вспомогательный буфер ячеек
    AmrCells cell_buffer_;

    // Маршрутизаторы для отправки примитивов из cells
    Router cell_router_;
    Router face_router_;
    Router vert_router_;

    /// @brief Вспомогательный буфер узлов
    AmrNodes node_buffer_;

    // Маршрутизаторы для отправки примитивов из cells
    Router node_router_;
    Router inct_router_;
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
    cell_buffer_.resize(
            locals.n_cells(),
            locals.n_faces(),
            locals.n_verts());

    // Сортировка migrants по rank, получается нормальный AmrCells
    // стартовые индексы (?)
    auto cell_index = cell_router_.send_offset();
    auto face_index = face_router_.send_offset();
    auto vert_index = vert_router_.send_offset();

    // Кортеж указателей на буферы данных
    std::array<Buffer*, n_vars> data_src = locals.data[loc_vars];
    std::array<Buffer*, n_vars> data_dst = cell_buffer_.data[mig_vars];

    for (index_t ic = 0; ic < cell_buffer_.size(); ++ic) {
        int r = locals.rank[ic];

        index_t jc = cell_index[r];

        locals.copy_geom(ic, cell_buffer_, jc, face_index[r], vert_index[r]);

        // Копирование данных
        //data_dst[jc] = data_src[ic];
        for (int v = 0; v < n_vars; ++v) {
            data_src[v]->copy_data(ic, *data_dst[v], jc);
        }

        cell_index[r] += 1;
        face_index[r] += locals.faces.max_count(ic);
        vert_index[r] += locals.verts.max_count(ic);
    }

    /*
    if (mpi::master()) {
        std::cout << "locals.fb: " << locals.faces.offsets << "\n";
        std::cout << "locals.nb: " << locals.verts.offsets << "\n";

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
void Migration::migrate(Tourism& tourists, AmrCells& locals, Vars&&... vars) {
    using utils::Buffer;
    AmrCells& ghosts = tourists.ghosts();

    // Все дополнительные переменные имеют тип Storable<T>
    soa::assert_storable<Vars...>();

    // Число переменных для пересылки
    constexpr int n_vars = sizeof...(Vars);

    // Кортежи переменных std::tuple{Storable<T1>, Storable<T2>, ...}
    // loc_vars - переменные в locals
    // mig_vars - переменные в migrants
    auto loc_vars = std::tuple{std::forward<Vars>(vars)...};
    auto mig_vars = cell_buffer_.data.add_replace(locals.data, loc_vars);

    // Составим полные матрицы пересылок (совпадают на всех процессах).
    // То есть заполним cell_router, face_router и vert_router.
    fill_router(locals);

    // Переиндексировать ячейки, отправить новые индексы и ранги в ghost-ячейки
    reindexing(tourists, locals, ghosts);

    // Перенести ячейки из locals в migrants в нужном порядке
    fill_migrants(locals, loc_vars, mig_vars);

    // Меняем размеры locals
    locals.resize(cell_router_.recv_buffer_size(),
                  face_router_.recv_buffer_size(),
                  vert_router_.recv_buffer_size());

    // ========================== Пересылка ячеек =============================

    // Указатели на буферы данных
    std::array<Buffer*, n_vars> data_src = cell_buffer_.data[mig_vars];
    std::array<Buffer*, n_vars> data_dst = locals.data[loc_vars];

    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в faces.offsets и verts.offsets количество элементов на ячейку
    for (index_t ic = 0; ic < cell_buffer_.size(); ++ic) {
        cell_buffer_.faces.offsets[ic] = cell_buffer_.faces.offsets[ic + 1] - cell_buffer_.faces.offsets[ic];
        cell_buffer_.verts.offsets[ic] = cell_buffer_.verts.offsets[ic + 1] - cell_buffer_.verts.offsets[ic];
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
    verts_send += vert_router_.isend(cell_buffer_.verts.coords, MpiTag::VERT_COORD);
    if (cell_buffer_.verts.unique()) {
        verts_send += vert_router_.isend(cell_buffer_.verts.index, MpiTag::VERT_INDEX);
        verts_send += vert_router_.isend(cell_buffer_.verts.ghost, MpiTag::VERT_GHOST);
    }

    // ============================= IRECV ====================================

    // Получить геометрию ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += cell_router_.irecv(locals.rank, MpiTag::RANK);
    cells_recv += cell_router_.irecv(locals.next, MpiTag::NEXT);
    cells_recv += cell_router_.irecv(locals.index, MpiTag::INDEX);
    cells_recv += cell_router_.irecv(locals.flag, MpiTag::FLAG);
    cells_recv += cell_router_.irecv(locals.level, MpiTag::LEVEL);
    cells_recv += cell_router_.irecv(locals.b_idx, MpiTag::B_IDX);
    cells_recv += cell_router_.irecv(locals.z_idx, MpiTag::Z_IDX);
    cells_recv += cell_router_.irecv(locals.center, MpiTag::CENTER);
    cells_recv += cell_router_.irecv(locals.volume, MpiTag::VOLUME);
    cells_recv += cell_router_.irecv(locals.volume_alt, MpiTag::VOLUME_ALT);

    // Получить данные ячеек
    RequestsList data_recv; data_recv.reserve(n_vars);
    for (int v = 0; v < n_vars; ++v) {
        data_recv += cell_router_.irecv(*data_dst[v], static_cast<MpiTag>(12345 + v));
    };

    // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += cell_router_.irecv(locals.faces.offsets.data() + 1, MpiTag::FACE_BEG);
    faces_recv += cell_router_.irecv(locals.verts.offsets.data() + 1, MpiTag::VERT_BEG);

    // Получить данные граней
    faces_recv += face_router_.irecv(locals.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += face_router_.irecv(locals.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += face_router_.irecv(locals.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_recv += face_router_.irecv(locals.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += face_router_.irecv(locals.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_recv += face_router_.irecv(locals.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += face_router_.irecv(locals.faces.normal, MpiTag::NORMAL);
    faces_recv += face_router_.irecv(locals.faces.center, MpiTag::FACE_CENTER);
    faces_recv += face_router_.irecv(locals.faces.area, MpiTag::AREA);
    faces_recv += face_router_.irecv(locals.faces.area_alt, MpiTag::AREA_ALT);
    faces_recv += face_router_.irecv(locals.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    RequestsList verts_recv; verts_recv.reserve(3);
    verts_recv += vert_router_.irecv(locals.verts.coords, MpiTag::VERT_COORD);
    if (cell_buffer_.verts.unique()) {
        verts_recv += vert_router_.irecv(locals.verts.index, MpiTag::VERT_INDEX);
        verts_recv += vert_router_.irecv(locals.verts.ghost, MpiTag::VERT_GHOST);
    }

    // =========================== WAIT ISEND =================================

    cells_send.wait();  // Завершить отправку ячеек
    data_send.wait();   // Завершить отправку данных ячеек
    faces_send.wait();  // Завершить отправку граней
    verts_send.wait();  // Завершить отправку вершин

    // =========================== WAIT IRECV =================================

    cells_recv.wait();  // Завершить получение ячеек
    data_recv.wait();   // Завершить получение данных ячеек
    faces_recv.wait();  // Завершить получение граней
    verts_recv.wait();  // Завершить получение вершин

    // Восстановить индексацию граней
    locals.faces.offsets[0] = 0;
    locals.verts.offsets[0] = 0;
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        locals.faces.offsets[ic + 1] += locals.faces.offsets[ic];
        locals.verts.offsets[ic + 1] += locals.verts.offsets[ic];
    }

    // ========================================================================

    /*
    mpi::for_each([&]() {
        std::cout << "Rank " << mpi::rank() << "\n";
        std::cout << "  fb: " << locals.faces.offsets << "\n";
        std::cout << "  nb: " << locals.verts.offsets << "\n";
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