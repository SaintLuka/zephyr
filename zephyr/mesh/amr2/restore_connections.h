/// @file Файл содержит реализацию функции restore_connections, которая создает связи
/// для новых ячеек, созданных после адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/siblings.h>
#include <zephyr/mesh/amr2/statistics.h>

namespace zephyr::mesh::amr2 {

/// @brief Восстанавливает связи у одной ячейки
/// @param cell Целевая ячейка
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param rank Ранг процесса
/// @param count Статистика адаптации
/// @details Старые ячейки (существовшие до процедуры адаптации) в поле amrData.next
/// содержат данные о расположении ячеек, которые придут им на смену
/// (см. setup_positions). Если ячейка не меняется, тогда cells[ic].next = ic.
/// Если ячейка разбивается, тогда amrData.next содержит индекс первой дочерней
/// ячейки (остальные располагаются за ней). Если ячейка огрубляется, тогда
/// amrData.next содержит индекс родительской ячейки.
/// После проведения геометрической части адаптаиции (см. setup_geometry)
/// новые ячейки уже созданы и располагаются в хранилище, но ссылаются на
/// ячейки по старым индексам. К примеру, пусть ячейка ic не меняется и имеет
/// справа соседа такого же уровня с индексом in, который хочет разбиться.
/// После проведения setup_geometry ячейка ic продолжит ссылаться на in,
/// в то время как cells[in].next будет содержать индекс ребенка.
/// Таким образом, в соответствии с amrData.next ищутся актуальные ссылки на
/// реальных соседей, полученных после адаптации.
/// Грани через процессы остаются без изменений. Линковка граней между процессами
/// происходит на последнем этапе алгоритма адапатции.
template<int dim>
void restore_connections_one(index_t ic, AmrCells& locals, AmrCells& aliens, int rank, const Statistics &count) {
    // Ячейка неактуальна, пропускаем
    if (locals.is_undefined(ic)) { return; }

    auto& faces = locals.faces;

    for (index_t iface: locals.faces_range(ic)) {
        if (faces.is_undefined(iface) ||
            faces.is_boundary(iface)) {
            continue;
        }

        bool local_cell = faces.adjacent.is_local(iface);

#if SCRUTINY
        // Проверим индексы
        if (local_cell) {
            if (faces.adjacent.index[iface] < 0 ||
                faces.adjacent.index[iface] >= locals.size()) {
                throw std::runtime_error("Restore connections: out of range #1");
            }
            if (faces.adjacent.rank[iface] != rank) {
                throw std::runtime_error("Restore connections: bad adjacent rank #1");
            }
        }
        else {
            if (faces.adjacent.alien[iface] < 0) {
                throw std::runtime_error("Restore connections: out of range #2");
            }
            if (faces.adjacent.rank[iface] == rank) {
                locals.print_info(ic);
                throw std::runtime_error("Restore connections: bad adjacent rank #2");
            }
        }
#endif

        // Пропустить новые ячейки, если на них есть ссылка, значит всё
        // нормально, нам нужно переставить ссылки с ячеек, которые пропали.
        if (local_cell && faces.adjacent.index[iface] > count.n_cells) {
            continue;
        }

        // Индекс старого соседа
        auto [neibs, old_jc] = faces.adjacent.get_neib(iface, locals, aliens);

        // Актуальный сосед через грань
        if (neibs.is_actual(old_jc)) {
            continue;
        }

        // Сосед через грань ничего не делал - ничего не надо,
        // ссылки актуальны
        if (neibs.flag[old_jc] == 0) {
            continue;
        }

        // Сосед через грань огрубился
        if (neibs.flag[old_jc] < 1) {
            // Соседняя ячейка огрубляется
            faces.adjacent.index[iface] = neibs.next[old_jc];
            continue;
        }

#if SCRUTINY
        // Поиск в родительской ячейке, такого не должно быть
        if (neibs.b_idx[old_jc] == locals.b_idx[ic] &&
            neibs.level[old_jc] == locals.level[ic] - 1 &&
            neibs.z_idx[old_jc] == locals.z_idx[ic] / CpC(dim)) {
            std::cout << neibs.rank[old_jc] << " " << neibs.level[old_jc] << " " << neibs.b_idx[old_jc] << " " << neibs.z_idx[old_jc] << " ""\n";
            std::cout << locals.rank[ic] << " " << locals.level[ic] << " " << locals.b_idx[ic] << " " << locals.z_idx[ic] << " ""\n";
            throw std::runtime_error("Search neib in parent cell");
        }
#endif

        // Сосед через грань адаптируется

        // Индекс вершины ячейки, ближайшей к центру грани,
        // и есть индекс смежной дочерней ячейки (гениально)
        int nearest_j = -1;
        double dist = std::numeric_limits<double>::max();

        // Прикол, так можно?
        auto shape = neibs.mapping<dim>(old_jc).reduce();
        for (int j = 0; j < CpC(dim); ++j) {
            double loc_dist = (shape[j] - faces.center[iface]).squaredNorm();
            if (loc_dist < dist) {
                dist = loc_dist;
                nearest_j = j;
            }
        }
        faces.adjacent.index[iface] = neibs.next[old_jc] + nearest_j;

        // Но это не верно для периодического замыкания, дописать потом
        scrutiny_check(faces.boundary[iface] != Boundary::PERIODIC, "Can't handle periodic")

#if SCRUTINY
        // Проверяем, что здесь не ссылка на брата, брат должен быть установлен ранее
        if (neibs.b_idx[old_jc] == locals.b_idx[ic] &&
            neibs.level[old_jc] == locals.level[ic] - 1 &&
            neibs.z_idx[old_jc] == locals.z_idx[ic] / CpC(dim)) {

            auto [bros, bro_idx] = faces.adjacent.get_neib(iface, locals, aliens);

            if (locals.b_idx[ic] == bros.b_idx[bro_idx] && bros.z_idx[bro_idx] == locals.z_idx[ic]) {
                std::cout << neibs.rank[old_jc] << " " << neibs.level[old_jc] << " " << neibs.b_idx[old_jc] << " " << neibs.z_idx[old_jc] << " ""\n";
                std::cout << locals.rank[ic] << " " << locals.level[ic] << " " << locals.b_idx[ic] << " " << locals.z_idx[ic] << " ""\n";
                std::cout << bros.rank[bro_idx] << " " << bros.level[bro_idx] << " " << bros.b_idx[bro_idx] << " " << bros.z_idx[bro_idx] << " ""\n";
                throw std::runtime_error("Found same cell");
            }
            if (locals.b_idx[ic] == bros.b_idx[bro_idx] && neibs.z_idx[old_jc] == bros.z_idx[bro_idx] / CpC(dim)) {
                throw std::runtime_error("Found brother");
            }
        }

        if (!local_cell) {
            continue;
        }

        // Если ячейка локальная, то можем проверить полным обходом всех возможных соседей
        bool found = false;

        int adj_index = -1;

        for (int j = 0; j < CpC(dim); ++j) {
            index_t jc = locals.next[old_jc] + j;

            if (jc > count.n_cells_large) {
                continue;
            }

            if (ic == jc) { continue; }

            for (const auto &jface: locals.faces_range(jc)) {
                if (faces.is_undefined(jface))
                    continue;

                if ((faces.center[iface] - faces.center[jface]).norm() < 1.0e-5 * locals.linear_size(ic)) {
                    adj_index = jc;
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
        }
        if (!found) {
            Side<dim> side = iface - locals.face_begin[ic];
            std::cout << "Can't find neighbor through the " << side << "\n";
            locals.print_info(ic);
            throw std::runtime_error("Can't find neighbor");
        }

        if (faces.adjacent.index[iface] != adj_index) {
            std::cout << "mine.level: " << locals.level[ic] << "\n";
            std::cout << "neib.level: " << locals.level[old_jc] << "\n";

            std::cout << std::setprecision(2);
            std::cout << "nearest_j: " << nearest_j << "; " << faces.adjacent.index[iface] - locals.next[old_jc] << "\n";
            std::cout << "face1.n: " << faces.normal[iface].transpose() << "\n";
            std::cout << "face1.c: " << (faces.center[iface] - locals.center[old_jc]).transpose() << "\n";
            std::cout << "face1.c: " << (faces.center[iface]).transpose() << "\n";
            std::cout << "neib.c: " << locals.center[old_jc].transpose() << "\n";
            for (int k = 0; k < 8; ++k) {
                std::cout << "  " << k << " child: " << (locals.center[locals.next[old_jc] + k] - locals.center[old_jc]).transpose() << "\n";
            }
        }
#endif
    }
}

/// @brief Устанавливает поле AmrCell.next в неопределенное состояние
inline void set_undefined_next(index_t ic, std::vector<index_t> &next) {
    next[ic] = -1;
}

/// @brief Восстанавливает связи (без MPI)
/// @details см. restore_connections_one
template<int dim>
void restore_connections(AmrCells &locals, AmrCells& aliens, int rank, const Statistics &count) {
    threads::parallel_for(
            index_t{0}, index_t{count.n_cells_large},
            restore_connections_one<dim>,
            std::ref(locals), std::ref(aliens), rank, std::ref(count)
    );

    threads::parallel_for(
            index_t{0}, index_t{count.n_cells_large},
            set_undefined_next,
            std::ref(locals.next)
    );
}

} // namespace zephyr