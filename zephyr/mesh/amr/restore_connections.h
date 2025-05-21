/// @file Файл содержит реализацию функции restore_connections, которая создает связи
/// для новых ячеек, созданных после адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr::mesh::amr {

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
void restore_connections_one(AmrStorage::Item& cell,
        AmrStorage &locals, AmrStorage &aliens,
        int rank, const Statistics &count) {

    // Ячейка неактуальна, пропускаем
    if (cell.is_undefined()) {
        return;
    }

    for (int i = 0; i < BFaces::max_count; ++i) {
        auto &face1 = cell.faces[i];
        if (face1.is_undefined() || face1.is_boundary()) {
            continue;
        }

        bool local_cell = face1.adjacent.alien < 0;

#if SCRUTINY
        // Проверим индексы
        if (local_cell) {
            if (face1.adjacent.index < 0 || face1.adjacent.index >= locals.size()) {
                throw std::runtime_error("Restore connections: out of range #1");
            }
        }
        else {
            if (face1.adjacent.alien < 0 || face1.adjacent.alien >= aliens.size()) {
                throw std::runtime_error("Restore connections: out of range #2");
            }
        }
#endif

        // Пропустить новые ячейки, если на них есть ссылка, значит всё
        // нормально, нам нужно переставить ссылки с ячеек, которые пропали.
        if (local_cell && face1.adjacent.index > count.n_cells) {
            continue;
        }

        const auto &old_neib = local_cell ?
                locals[face1.adjacent.index] :
                aliens[face1.adjacent.alien];

        // Актуальный сосед через грань
        if (old_neib.is_actual()) {
            continue;
        }

        // Сосед через грань ничего не делал - ничего не надо,
        // ссылки актуальны
        if (old_neib.flag == 0) {
            continue;
        }

        // Сосед через грань огрубился
        if (old_neib.flag < 1) {
            // Соседняя ячейка огрубляется
            face1.adjacent.index = old_neib.next;
            continue;
        }

#if SCRUTINY
        // Поиск в родительской ячейке, такого не должно быть
        if (old_neib.b_idx == cell.b_idx &&
            old_neib.level == cell.level - 1 &&
            old_neib.z_idx == cell.z_idx / CpC(dim)) {
            std::cout << old_neib.rank << " " << old_neib.level << " " << old_neib.b_idx << " " << old_neib.z_idx << " ""\n";
            std::cout << cell.rank << " " << cell.level << " " << cell.b_idx << " " << cell.z_idx << " ""\n";
            throw std::runtime_error("Search neib in parent cell");
        }
#endif

        // Сосед через грань адаптируется

        // Индекс вершины ячейки, ближайшей к центру грани,
        // и есть индекс смежной дочерней ячейки (гениально)
        int nearest_j = -1;
        double dist = std::numeric_limits<double>::max();
        if (dim == 2) {
            Quad quad = old_neib.vertices.as2D().reduce();
            for (int j = 0; j < CpC(dim); ++j) {
                double loc_dist = (quad[j] - face1.center).squaredNorm();
                if (loc_dist < dist) {
                    dist = loc_dist;
                    nearest_j = j;
                }
            }
        }
        else {
            Cube cube = old_neib.vertices.reduce();
            for (int j = 0; j < CpC(dim); ++j) {
                double loc_dist = (cube[j] - face1.center).squaredNorm();
                if (loc_dist < dist) {
                    dist = loc_dist;
                    nearest_j = j;
                }
            }
        }
        face1.adjacent.index = old_neib.next + nearest_j;

        // Но это не верно для периодического замыкания, дописать потом
        scrutiny_check(face1.boundary != Boundary::PERIODIC, "Can't handle periodic")

#if SCRUTINY
        // Проверяем, что здесь не ссылка на брата, брат должен быть установлен ранее
        if (old_neib.b_idx == cell.b_idx &&
            old_neib.level == cell.level - 1 &&
            old_neib.z_idx == cell.z_idx / CpC(dim)) {
            auto& bro = locals[face1.adjacent.index];

            if (cell.b_idx == bro.b_idx && bro.z_idx == cell.z_idx) {
                std::cout << old_neib.rank << " " << old_neib.level << " " << old_neib.b_idx << " " << old_neib.z_idx << " ""\n";
                std::cout << cell.rank << " " << cell.level << " " << cell.b_idx << " " << cell.z_idx << " ""\n";
                std::cout << bro.rank << " " << bro.level << " " << bro.b_idx << " " << bro.z_idx << " ""\n";
                throw std::runtime_error("Found same cell");
            }
            if (cell.b_idx == bro.b_idx && old_neib.z_idx == bro.z_idx / CpC(dim)) {
                throw std::runtime_error("Found brother");
            }
        }

        if (!local_cell) {
            continue;
        }

        // Если ячейка локальная, то можем проверить полным обходом
        // всех возможных соседей
        bool found = false;

        int adj_index = -1;

        for (int j = 0; j < CpC(dim); ++j) {
            auto jc = old_neib.next + j;

            if (jc > count.n_cells_large) {
                continue;
            }

            const auto &neib = locals[jc];
            if (&neib == &cell) {
                continue;
            }

            for (const auto &face2: neib.faces) {
                if (face2.is_undefined())
                    continue;

                if ((face1.center - face2.center).norm() < 1.0e-5 * cell.linear_size()) {
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
            std::cout << "Can't find neighbor through the " << side_to_string(i, dim)
                      << " face (" << i / 6 << ")\n";
            //amr::print_cell_info(locals, aliens, ic);
            throw std::runtime_error("Can't find neighbor");
        }

        if (face1.adjacent.index != adj_index) {
            std::cout << "mine.level: " << cell.level << "\n";
            std::cout << "neib.level: " << old_neib.level << "\n";

            std::cout << std::setprecision(2);
            std::cout << "nearest_j: " << nearest_j << "; " << face1.adjacent.index - old_neib.next << "\n";
            std::cout << "face1.n: " << face1.normal.transpose() << "\n";
            std::cout << "face1.c: " << (face1.center - old_neib.center).transpose() << "\n";
            std::cout << "face1.c: " << (face1.center).transpose() << "\n";
            std::cout << "neib.c: " << (old_neib.center).transpose() << "\n";
            for (int k = 0; k < 8; ++k) {
                std::cout << "  " << k << " child: " << (locals[old_neib.next + k].center - old_neib.center).transpose() << "\n";
                std::cout << "  " << k << " vertx: " << (old_neib.vertices.reduce()[k] - old_neib.center).transpose() << "\n";
            }
        }
#endif
    }
}

/// @brief Устанавливает поле AmrCell.next в неопределенное состояние
void set_undefined_next(AmrStorage::Item& cell) {
    cell.next = -1;
}

/// @brief Восстанавливает связи (без MPI)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(
        AmrStorage &locals, AmrStorage& aliens,
        int rank, const Statistics &count) {

    threads::for_each<20>(
            locals.begin(),
            locals.iterator(count.n_cells_large),
            restore_connections_one<dim>,
            std::ref(locals),
            std::ref(aliens),
            rank, std::ref(count)
    );

    threads::for_each(
            locals.begin(),
            locals.iterator(count.n_cells_large),
            set_undefined_next);
}

} // namespace zephyr