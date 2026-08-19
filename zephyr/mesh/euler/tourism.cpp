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
    //unique_border_indices_.shrink_to_fit();
    border_indices_.shrink_to_fit();
    border_.shrink_to_fit();
    ghosts_.shrink_to_fit();
}

void Tourism::init_types(const AmrCells& locals) {
    border_ = locals.same();
    ghosts_ = locals.same();
}

void Tourism::resize_border() {
    index_t n_border_cells = cell_router_.send_buffer_size();
    index_t n_border_faces = face_router_.send_buffer_size();
    index_t n_border_verts = vert_router_.send_buffer_size();

    border_.resize(n_border_cells, n_border_faces, n_border_verts);
}

void Tourism::extend_border() {
    index_t n_border_cells = cell_router_.send_buffer_size();
    index_t n_border_faces = face_router_.send_buffer_size();
    index_t n_border_verts = vert_router_.send_buffer_size();

    if (n_border_cells > border_.size()) {
        border_.resize(n_border_cells, n_border_faces, n_border_verts);
    }
}

void Tourism::resize_ghosts() {
    int n_ghost_cells = cell_router_.recv_buffer_size();
    int n_ghost_faces = face_router_.recv_buffer_size();
    int n_ghost_verts = vert_router_.recv_buffer_size();

    ghosts_.resize(n_ghost_cells, n_ghost_faces, n_ghost_verts);
}

void Tourism::extend_ghosts() {
    int n_ghost_cells = cell_router_.recv_buffer_size();
    int n_ghost_faces = face_router_.recv_buffer_size();
    int n_ghost_verts = vert_router_.recv_buffer_size();

    if (n_ghost_cells > ghosts_.size()) {
        ghosts_.resize(n_ghost_cells, n_ghost_faces, n_ghost_verts);
    }
}

void Tourism::fill_send_count(const AmrCells& locals) {
    const int size = mpi::size();
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    std::vector<index_t> last_append(size, -1);

    std::vector<index_t> cell_send_count(size, 0);
    std::vector<index_t> face_send_count(size, 0);
    std::vector<index_t> vert_send_count(size, 0);

    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (locals.is_undefined(ic)) continue;

        for (index_t iface: locals.faces.range(ic)) {
            if (locals.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = locals.faces.adjacent.rank[iface];
            if (neib_rank != rank && last_append[neib_rank] != ic) {
                last_append[neib_rank] = ic;
                cell_send_count[neib_rank] += 1;
                face_send_count[neib_rank] += locals.faces.max_count(ic);
                vert_send_count[neib_rank] += locals.verts.max_count(ic);
            }
        }
    }

    // Установить число на обмены
    cell_router_.set_send_count(cell_send_count);
    face_router_.set_send_count(face_send_count);
    vert_router_.set_send_count(vert_send_count);
}

void Tourism::fill_indices(const AmrCells& locals) {
    const int rank = mpi::rank();

    // Индекс последней учтенной ячейки
    index_t last_unique_append = -1;
    std::vector<index_t> last_append(mpi::size(), -1);

    // Смещения, по которым записываются индексы
    std::vector<index_t> cell_index = cell_router_.send_offset();

    border_indices_.resize(cell_router_.send_buffer_size());
    // unique_border_indices_.reserve(m_cell_route.send_buffer_size());

    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        // Для сетки с неактуальными ячейками
        if (locals.is_undefined(ic)) continue;

        for (index_t iface: locals.faces.range(ic)) {
            if (locals.faces.is_undefined(iface)) {
                continue;
            }

            int neib_rank = locals.faces.adjacent.rank[iface];
            if (neib_rank != rank) {
                if (last_append[neib_rank] != ic) {
                    last_append[neib_rank] = ic;
                    border_indices_[cell_index[neib_rank]++] = ic;
                }
                if (last_unique_append != ic) {
                    last_unique_append = ic;
                    // unique_border_indices_.push_back(ic);
                }
            }
        }
    }
}

void Tourism::prepare_geometry(const AmrCells& locals) {
    index_t face_idx = 0;
    index_t vert_idx = 0;
    for (index_t ic = 0; ic < border_indices_.size(); ++ic) {
        locals.copy_geom(border_indices_[ic], border_, ic, face_idx, vert_idx);

        face_idx += locals.faces.max_count(border_indices_[ic]);
        vert_idx += locals.verts.max_count(border_indices_[ic]);
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

// Инициализация индекса ghost = -1 для большинства граней
void set_undef_ghosts(AmrCells& locals, int rank) {
    threads::parallel_for(
        index_t{0}, locals.n_cells(),
        [&locals, rank](index_t ic) {
            for (index_t iface: locals.faces.range(ic)) {
                if (locals.faces.is_actual(iface) &&
                    locals.faces.adjacent.rank[iface] == rank) {
                    locals.faces.adjacent.ghost[iface] = -1;
                }
            }
        });
}

// Обходим ячейки в ghost и ищем связи
void Tourism::find_connections(AmrCells& locals, int rank) const {
    for (index_t ic = 0; ic < ghosts_.n_cells(); ++ic) {
        for (index_t iface: ghosts_.faces.range(ic)) {
            if (ghosts_.faces.is_undefined(iface)) {
                continue;
            }

            if (ghosts_.faces.adjacent.index[iface] >= 0 &&
                ghosts_.faces.adjacent.rank[iface] == rank) {
                // Индекс соседа
                index_t jc = ghosts_.faces.adjacent.index[iface];

                for (index_t l_face: locals.faces.range(jc)) {
                    if (locals.faces.adjacent.rank [l_face] == ghosts_.rank [ic] &&
                        locals.faces.adjacent.index[l_face] == ghosts_.index[ic]) {

                        locals.faces.adjacent.ghost[l_face] = ic;
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
    cell_router_.fill_partial();
    face_router_.fill_partial();
    vert_router_.fill_partial();

    // Расширить массив ghosts для получения геометрии
    resize_ghosts();

    // Отправить и получить геометрию
    sync_geometry();

    // Инициализация индекса ghost = -1 для большинства граней
    set_undef_ghosts(locals, mpi::rank());

    // Обходим ячейки в ghost и ищем связи
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
        // i - индекс в border_, border_indices_
        for (index_t i: cell_router_.send_indices(r)) {
            if (border_.flag[i] == 0) {
                border_.next[i] = next_index;
                next_index += 1;
            }
            else if (border_.flag[i] == 1) {
                // bitset<8> для дочерних ячеек
                auto children = dim == 2 ? border_children<2>(border_.faces, i, r) :
                                           border_children<3>(border_.faces, i, r);
                // Кодируем список дочерних ячеек
                border_.next[i] = amr::pack_children(next_index, children);
                next_index += static_cast<index_t>(children.count());
            }
            else {
                // Самый неприятный случай, border-ячейка огрубляется

                // Полный индекс родительской ячейки
                std::tuple<index_t, index_t, index_t> parent = {
                    border_.b_idx[i],
                    border_.level[i],
                    border_.z_idx[i] / indexing::CpC(dim),
                };

                auto parent_it = coarse_cells.find(parent);
                if (parent_it != coarse_cells.end()) {
                    border_.next[i] = parent_it->second;
                }
                else {
                    coarse_cells[parent] = next_index;
                    border_.next[i] = next_index;
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
    std::vector<index_t> prev_border_indices = border_indices_;
    border_indices_.resize(cell_router_.send_buffer_size());

    index_t last_border_next = 0;
    for (index_t i = 0; i < prev_border_indices.size(); ++i) {
        z_assert(i < border_.flag.size(), "out of range #1521");
        z_assert(i < border_.next.size(), "out of range #1522");
        z_assert(i < prev_border_indices.size(), "out of range #1523");

        if (border_.flag[i] == 0) {
            index_t border_next = border_.next[i];

            z_assert(border_next < border_indices_.size(), "out of range #1524");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1525");

            border_indices_[border_next] = locals_next[prev_border_indices[i]];
            last_border_next = std::max(last_border_next, border_next);
        }
        else if (border_.flag[i] < 0) {
            index_t border_next = border_.next[i];

            z_assert(border_next < border_indices_.size(), "out of range #1526");
            z_assert(prev_border_indices[i] < locals_next.size(), "out of range #1527");

            index_t parent_index = locals_next[prev_border_indices[i]];

            z_assert(parent_index < locals_next.size(), "out of range #1528");

            border_indices_[border_next] = locals_next[parent_index];
            last_border_next = std::max(last_border_next, border_next);
        }
        else {
            auto [border_next, children] = amr::unpack_children(border_.next[i]);
            index_t main_child = locals_next[prev_border_indices[i]];
            for (int c = 0; c < indexing::CpC(dim); ++c) {
                if (children[c]) {
                    if (border_next >= border_indices_.size()) {
                        std::cout << border_.next[i] << "; " << border_next << "; " << children << "; " << border_indices_.size() << "\n";
                    }
                    z_assert(border_next < border_indices_.size(), "out of range #1529");
                    z_assert(main_child + c < locals_next.size(), "out of range #1530");
                    border_indices_[border_next] = locals_next[main_child + c];
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

    bef_border.save(border_, pvd_counter);
    bef_ghosts.save(ghosts_, pvd_counter);
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
            border_.next[i] += cell_router_.send_offset(r);
        }
    }

    // Добавляем смещения, теперь индексы NEXT в ghost идут последовательно (за
    // исключением закодированных индексов для ячеек на разбиение). Для каждой
    // ghost-ячейки указана следующая позиция внутри нового ghost-слоя.
    for (int r = 0; r < mpi::size(); ++r) {
        if (r == rank) { continue; }
        for (index_t i: prev_router.recv_indices(r)) {
            ghosts_.next[i] += cell_router_.recv_offset(r);
        }
    }

#if WRITE_DBG
    aft_border.save(border_, pvd_counter);
    aft_ghosts.save(ghosts_, pvd_counter);
    ++pvd_counter;
#endif

    // ========================================================================
    //           Подготовим новые массивы, без актуальных данных
    // ========================================================================

    // Подготовить border массив
    extend_border();

    // Подготовить ghosts массив
    extend_ghosts();

    // Выставить корректные индексы в border_indices_
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

void Tourism::pack_border_indices() {
    // Оптимизируем использование памяти, используем повторно массивы.
    // Запишем в faces.offsets и verts.offsets количество элементов на ячейку
    for (index_t ic = 0; ic < border_.size(); ++ic) {
        border_.faces.offsets[ic] = border_.faces.offsets[ic + 1] - border_.faces.offsets[ic];
        border_.verts.offsets[ic] = border_.verts.offsets[ic + 1] - border_.verts.offsets[ic];
    }
    border_.faces.offsets.back() = -1;
    border_.verts.offsets.back() = -1;
}

void Tourism::unpack_border_indices() {
    int prev_n_face = border_.faces.offsets[0];
    int prev_n_vert = border_.verts.offsets[0];

    border_.faces.offsets[0] = 0;
    border_.verts.offsets[0] = 0;

    for (index_t ic = 1; ic <= border_.n_cells(); ++ic) {
        int temp_n_face = border_.faces.offsets[ic];
        int temp_n_vert = border_.verts.offsets[ic];

        border_.faces.offsets[ic] = border_.faces.offsets[ic - 1] + prev_n_face;
        border_.verts.offsets[ic] = border_.verts.offsets[ic - 1] + prev_n_vert;

        prev_n_face = temp_n_face;
        prev_n_vert = temp_n_vert;
    }
}

void Tourism::unpack_ghost_indices() {
    // Восстанавливаем индексацию, сейчас в массивах хранится число граней или вершин ячейки
    ghosts_.faces.offsets[0] = 0;
    ghosts_.verts.offsets[0] = 0;
    for (index_t ic = 0; ic < ghosts_.n_cells(); ++ic) {
        ghosts_.faces.offsets[ic + 1] += ghosts_.faces.offsets[ic];
        ghosts_.verts.offsets[ic + 1] += ghosts_.verts.offsets[ic];
    }
}

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

void Tourism::send_geometry(const AmrCells& locals) {
    prepare_geometry(locals);
    sync_geometry();
}

void Tourism::restore_indices(AmrCells& locals) const {
    for (index_t ic: border_indices_) {
        for (index_t iface: locals.faces.range(ic)) {
            index_t ghost_index = locals.faces.adjacent.ghost[iface];
            if (ghost_index >= 0) {
                locals.faces.adjacent.index[iface] = ghosts_.index[ghost_index];
            }
        }
    }
}

void Tourism::sync_geometry() {
    bool amr = border_.adaptive();
    bool axial = border_.axial();

    if (!amr) {
        // Оптимизируем пересылку индексов граней/вершин
        pack_border_indices();
    }

    // ============================= ISEND ====================================

    // Отправить данные ячеек
    RequestsList cells_send; cells_send.reserve(16);
    cells_send += cell_router_.isend(border_.rank, MpiTag::RANK);
    cells_send += cell_router_.isend(border_.next, MpiTag::NEXT);
    cells_send += cell_router_.isend(border_.index, MpiTag::INDEX);
    cells_send += cell_router_.isend(border_.flag, MpiTag::FLAG);
    cells_send += cell_router_.isend(border_.level, MpiTag::LEVEL);
    cells_send += cell_router_.isend(border_.b_idx, MpiTag::B_IDX);
    cells_send += cell_router_.isend(border_.z_idx, MpiTag::Z_IDX);
    cells_send += cell_router_.isend(border_.center, MpiTag::CENTER);
    cells_send += cell_router_.isend(border_.volume, MpiTag::VOLUME);
    if (axial) {
        cells_send += cell_router_.isend(border_.volume_alt, MpiTag::VOLUME_ALT);
    }
    if (!amr) {
        cells_send += cell_router_.isend(border_.faces.offsets, MpiTag::FACE_BEG);
        cells_send += cell_router_.isend(border_.verts.offsets, MpiTag::VERT_BEG);
    }

    // Отправить данные граней
    RequestsList faces_send; faces_send.reserve(16);
    faces_send += face_router_.isend(border_.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_send += face_router_.isend(border_.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_send += face_router_.isend(border_.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_send += face_router_.isend(border_.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_send += face_router_.isend(border_.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_send += face_router_.isend(border_.faces.boundary, MpiTag::BOUNDARY);
    faces_send += face_router_.isend(border_.faces.normal, MpiTag::NORMAL);
    faces_send += face_router_.isend(border_.faces.center, MpiTag::FACE_CENTER);
    faces_send += face_router_.isend(border_.faces.area, MpiTag::AREA);
    if (axial) {
        faces_send += face_router_.isend(border_.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_send += face_router_.isend(border_.faces.vertices, MpiTag::FACE_VERTS);

    // Отправить вершины
    RequestsList verts_send; verts_send.reserve(3);
    verts_send += vert_router_.isend(border_.verts.coords, MpiTag::VERT_COORD);
    if (border_.verts.unique()) {
        verts_send += vert_router_.isend(border_.verts.index, MpiTag::VERT_INDEX);
        verts_send += vert_router_.isend(border_.verts.ghost, MpiTag::VERT_GHOST);
    }

    // ============================= IRECV ====================================

    // Получить данные ячеек
    RequestsList cells_recv; cells_recv.reserve(16);
    cells_recv += cell_router_.irecv(ghosts_.rank, MpiTag::RANK);
    cells_recv += cell_router_.irecv(ghosts_.next, MpiTag::NEXT);
    cells_recv += cell_router_.irecv(ghosts_.index, MpiTag::INDEX);
    cells_recv += cell_router_.irecv(ghosts_.flag, MpiTag::FLAG);
    cells_recv += cell_router_.irecv(ghosts_.level, MpiTag::LEVEL);
    cells_recv += cell_router_.irecv(ghosts_.b_idx, MpiTag::B_IDX);
    cells_recv += cell_router_.irecv(ghosts_.z_idx, MpiTag::Z_IDX);
    cells_recv += cell_router_.irecv(ghosts_.center, MpiTag::CENTER);
    cells_recv += cell_router_.irecv(ghosts_.volume, MpiTag::VOLUME);
    if (axial) {
        cells_recv += cell_router_.irecv(ghosts_.volume_alt, MpiTag::VOLUME_ALT);
    }

    if (!amr) {
        // При получении используем сдвиг на единицу, чтобы записать нулевой первый элемент
        cells_recv += cell_router_.irecv(ghosts_.faces.offsets.data() + 1, MpiTag::FACE_BEG);
        cells_recv += cell_router_.irecv(ghosts_.verts.offsets.data() + 1, MpiTag::VERT_BEG);
    }

    // Получить данные граней
    RequestsList faces_recv; faces_recv.reserve(16);
    faces_recv += face_router_.irecv(ghosts_.faces.adjacent.rank, MpiTag::ADJ_RANK);
    faces_recv += face_router_.irecv(ghosts_.faces.adjacent.index, MpiTag::ADJ_INDEX);
    faces_recv += face_router_.irecv(ghosts_.faces.adjacent.ghost, MpiTag::ADJ_GHOST);
    faces_recv += face_router_.irecv(ghosts_.faces.adjacent.basic, MpiTag::ADJ_BASIC);
    faces_recv += face_router_.irecv(ghosts_.faces.adjacent.rotation, MpiTag::ADJ_ROTATION);
    faces_recv += face_router_.irecv(ghosts_.faces.boundary, MpiTag::BOUNDARY);
    faces_recv += face_router_.irecv(ghosts_.faces.normal, MpiTag::NORMAL);
    faces_recv += face_router_.irecv(ghosts_.faces.center, MpiTag::FACE_CENTER);
    faces_recv += face_router_.irecv(ghosts_.faces.area, MpiTag::AREA);
    if (axial) {
        faces_recv += face_router_.irecv(ghosts_.faces.area_alt, MpiTag::AREA_ALT);
    }
    faces_recv += face_router_.irecv(ghosts_.faces.vertices, MpiTag::FACE_VERTS);

    // Получить вершины
    RequestsList verts_recv; verts_recv.reserve(3);
    verts_recv += vert_router_.irecv(ghosts_.verts.coords, MpiTag::VERT_COORD);
    if (ghosts_.verts.unique()) {
        verts_recv += vert_router_.irecv(ghosts_.verts.index, MpiTag::VERT_INDEX);
        verts_recv += vert_router_.irecv(ghosts_.verts.ghost, MpiTag::VERT_GHOST);
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
        unpack_ghost_indices();

        // Поддерживать индексацию граней в border массиве не обязательно
        // unpack_border_indices();
    }
    else {
        if (border_.dim() == 2) {
            set_amr_indices<2>(ghosts_.faces.offsets, ghosts_.verts.offsets);
        }
        else {
            set_amr_indices<3>(ghosts_.faces.offsets, ghosts_.verts.offsets);
        }
    }
}

} // namespace zephyr::mesh

#endif