/// @file Файл содержит реализацию одной из базовых функций адаптации coarse_cell,
/// которая объединяет дочерние ячейки в родительскую.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/cell.h>
#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>
#include <zephyr/mesh/amr/siblings.h>

namespace zephyr { namespace mesh { namespace amr {


/// @return Массив неопределенных итераторов на дочерние ячейки
template <int dim>
std::array<Storage::Item, CpC(dim)> undefined_iterators(Storage& cells);

/// @return Массив неопределенных итераторов на дочерние ячейки (2D)
template <>
std::array<Storage::Item, CpC(2)> undefined_iterators<2>(Storage& cells) {
    return {cells.end(), cells.end(), cells.end(), cells.end()};
}

/// @return Массив неопределенных итераторов на дочерние ячейки (3D)
template <>
std::array<Storage::Item, CpC(3)> undefined_iterators<3>(Storage& cells) {
    return {cells.end(), cells.end(), cells.end(), cells.end(),
            cells.end(), cells.end(), cells.end(), cells.end()};
}

/// @brief Возвращает итераторы дочерних ячеек по их индексам в хранилище
/// @param cells Хранилище ячеек
/// @param ic Индекс главной из дочерних ячеек (z_loc = 0)
/// @return Массив итераторов дочерних ячеек
template <int dim>
std::array<Storage::Item, CpC(dim)> select_children(Storage& cells, int ic) {
    auto sibs = get_siblings<dim>(cells, ic);

    // Дочерние ячейки, упорядоченные по локальному z-индексу
    auto children = undefined_iterators<dim>(cells);
    children[0] = cells[ic];
    for (auto sib: sibs) {
        scrutiny_check(sib < cells.size(), "CellsAround: Sibling index out of range")

        auto z_loc = cells[sib].z_idx() % CpC(dim);
        children[z_loc] = cells[sib];
    }

#if SCRUTINY
    if (cells[ic].z_idx() % CpC(dim) != 0) {
        throw std::runtime_error("Not main child collect siblings");
    }
    auto main_lvl = cells[ic].level();

    std::set<int> found;
    found.insert(0);
    for (auto sib: sibs) {
        if (cells[sib].level() != main_lvl) {
            throw std::runtime_error("Different levels (siblings)");
        }
        int loc_z = cells[sib].z_idx() % CpC(dim);
        found.insert(loc_z);
    }
    int counter = 0;
    for (int i = 0; i < CpC(dim); ++i) {
        counter += found.count(i);
    }
    if (counter != CpC(dim)) {
        throw std::runtime_error("Wrong siblings");
    }
#endif

    return children;
}

/// @brief Список вершин родительской ячейки
LargeList2D parent_vs(const std::array<Storage::Item, 4>& children) {
     return {
            children[0].geom().vertices[iww(0, 0)],
            children[0].geom().vertices[iww(2, 0)],
            children[1].geom().vertices[iww(2, 0)],
            children[0].geom().vertices[iww(0, 2)],
            children[0].geom().vertices[iww(2, 2)],
            children[1].geom().vertices[iww(2, 2)],
            children[2].geom().vertices[iww(0, 2)],
            children[2].geom().vertices[iww(2, 2)],
            children[3].geom().vertices[iww(2, 2)]
    };
}

/// @brief Список вершин родительской ячейки
ShortList3D parent_vs(const std::array<Storage::Item, 8>& children) {
    return {
            children[iss(0, 0, 0)].geom().vertices[isw(0, 0, 0)],
            children[iss(1, 0, 0)].geom().vertices[isw(1, 0, 0)],
            children[iss(0, 1, 0)].geom().vertices[isw(0, 1, 0)],
            children[iss(1, 1, 0)].geom().vertices[isw(1, 1, 0)],
            children[iss(0, 0, 1)].geom().vertices[isw(0, 0, 1)],
            children[iss(1, 0, 1)].geom().vertices[isw(1, 0, 1)],
            children[iss(0, 1, 1)].geom().vertices[isw(0, 1, 1)],
            children[iss(1, 1, 1)].geom().vertices[isw(1, 1, 1)]
    };
}

/// @brief Создает родительскую ячейку с учетом окружения, то есть полученная
/// ячейка будет иметь актуальные ссылки на (старых) соседей, также полученная
/// ячейка может иметь более одной грани по каждой из сторон
/// @param locals Локальное хранилище ячеек
/// @param rank Ранг текущего процесса
/// @param children Массив дочерних ячеек
/// @return Полностью готовая родительская ячейка
template<int dim>
Cell get_parent(Storage &locals, Storage &aliens, int rank,
                     std::array<Storage::Item, CpC(dim)> children) {
    Cell parent(parent_vs(children));

    parent.b_idx = children[0].b_idx();
    parent.flag = 0;
    parent.level = children[0].level() - 1;
    parent.z_idx = children[0].z_idx() / CpC(dim);

    auto children_by_side = get_children_by_side<dim>();

    /// Линкуем грани
    for (int side = 0; side < FpC(dim); ++side) {
        auto some_child = children[children_by_side[side][0]];
        auto &some_face = some_child.geom().faces[side];

#if SCRUTINY
        if (some_face.boundary == FaceFlag::UNDEFINED) {
            throw std::runtime_error("Undefined boundary (coarse cell");
        }
        for (int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto &face = child.geom().faces[side];

            if (face.boundary != some_face.boundary) {
                throw std::runtime_error("Different boundary conditions");
            }
        }
#endif
        parent.faces[side].boundary = some_face.boundary;

        // Внешняя граница, не требуется линковать
        if (some_face.boundary != FaceFlag::ORDINARY &&
            some_face.boundary != FaceFlag::PERIODIC) {
            parent.faces[side].adjacent.rank = rank;
            parent.faces[side].adjacent.ghost = -1;
            continue;
        }

        auto some_neib_rank = some_face.adjacent.rank;
        auto some_neib_index = some_face.adjacent.index;
        auto some_neib_ghost = some_face.adjacent.ghost;
#if SCRUTINY
        if (some_neib_rank == rank && some_neib_index >= locals.size()) {
            std::cout << "Child has no local neighbor through the " <<
                      side_to_string(side % 6) << " side #1\n";
            some_child.print_info();
            throw std::runtime_error("Child has no local neighbor (coarse_cell) #1");
        }
        if (some_neib_rank != rank && some_neib_ghost >= aliens.size()) {
            std::cout << "Child has no remote neighbor through the " <<
                      side_to_string(side % 6) << " side #1\n";
            some_child.print_info();
            throw std::runtime_error("Child has no remote neighbor (coarse_cell) #1");
        }
#endif

        auto some_neib = some_neib_rank == rank ? locals[some_neib_index] : aliens[some_neib_ghost];
        auto some_neib_wanted_lvl = some_neib.level() + some_neib.flag();

#if SCRUTINY
        for (int i = 1; i < VpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto adj = child.geom().faces[side].adjacent;

            if (adj.rank == rank && adj.index >= locals.size()) {
                std::cout << "Child has no local neighbor through the " <<
                          side_to_string(side % 6) << " side #2\n";
                child.print_info();
                throw std::runtime_error("Child has no local neighbor (coarse_cell) #2");
            }
            if (adj.rank != rank && adj.ghost >= aliens.size()) {
                std::cout << "Child has no remote neighbor through the " <<
                          side_to_string(side % 6) << " side #2\n";
                child.print_info();
                throw std::runtime_error("Child has no remote neighbor (coarse_cell) #2");
            }

            auto neib = adj.rank == rank ? locals[adj.index] : aliens[adj.ghost];

            auto neib_wanted_lvl = neib.level() + neib.flag();
            if (neib_wanted_lvl != some_neib_wanted_lvl) {
                throw std::runtime_error("Different wanted level (coarsing cell neighbor)");
            }
        }
#endif

        // Простой случай: грань родительской ячейки должна быть простой
        if (some_neib_wanted_lvl <= parent.level) {
            parent.faces[side].adjacent.rank = some_neib_rank;
            parent.faces[side].adjacent.index = some_neib_index;
            parent.faces[side].adjacent.ghost = some_neib_ghost;
            continue;
        }

        // Если мы здесь, то сторона side от родителя должна адаптироваться
        // Далее потребуется связать грани

        split_face<dim>(parent, side);

        // Центры подграней родительской ячейки
        std::array<Vector3d, FpF(dim)> pfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            Face& face = parent.faces[side + 6*i];
            pfaces[i] = face.center<dim>(parent.vertices);
        }

        // Центры граней дочерних ячеек
        std::array<Vector3d, FpF(dim)> cfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            if (child.geom().faces[side + 6].is_undefined()) {
                // Простая грань
                throw std::runtime_error("olololshja");
                //? cfaces[i] = child.geom().faces[side].center<dim>(child.geom().vertices);
            } else {
                // Сложная грань, считаем центр по основным вершинам
                cfaces[i] = Vector3d(0.0, 0.0, 0.0);
                for (int j = 0; j < FpF(dim); ++j) {
                    //cfaces[i] += face_center<dim>(child.geom().faces[side + 6 * j], child[vertices]);
                    cfaces[i] += (Vector3d&) child.geom().vertices[child.geom().faces[side + 6 *j].vertices[j]];
                }
                cfaces[i] /= FpF(dim);

            }
        }

        /// Находим соответстиве между pfaces и cfaces
        double eps = 1.0e-6 * parent.size;
        for (int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto& child_face = child.geom().faces[side];
            for (int j = 0; j < FpF(dim); ++j) {
                if ((pfaces[i] - cfaces[j]).norm() < eps) {
                    parent.faces[side + 6*j].adjacent.rank = child_face.adjacent.rank;
                    parent.faces[side + 6*j].adjacent.index = child_face.adjacent.index;
                    parent.faces[side + 6*j].adjacent.ghost = child_face.adjacent.ghost;
                    break;
                }
            }
        }

#if SCRUTINY
        std::set<int> found;
        for (int i = 0; i < FpF(dim); ++i) {
            for (int j = 0; j < FpF(dim); ++j) {
                if ((pfaces[i] - cfaces[j]).norm() < eps) {
                    found.insert(j);
                    break;
                }
            }
        }
        /*
        if (found.size() != FpF(dim)) {
            std::cout << "side: " << side_to_string(side) << "\n";
            std::cout << std::scientific << std::setprecision(6) << "\n";
            std::cout << parent2.vertices.list[parent2.faces[side].vertices[0]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[1]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[2]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[3]] << "\nchildrend:\n";

            for (int ck = 0; ck < 8; ++ck) {
                auto child = children[ck];
                std::cout << child.geom().vertices[child.geom().faces[side].vertices[0]] << ", "
                          << child.geom().vertices[child.geom().faces[side].vertices[1]] << ", "
                          << child.geom().vertices[child.geom().faces[side].vertices[2]] << ", "
                          << child.geom().vertices[child.geom().faces[side].vertices[3]] << "\n";
            }

            for (int i = 0; i < FpF(dim); ++i) {
                std::cout << "pf: " << pfaces[i] << "\n";
                for (int j = 0; j < FpF(dim); ++j) {
                    std::cout << "  cf: " << cfaces[j] << "\n";
                    if (distance(pfaces[i], cfaces[j]) < eps) {
                        found.insert(j);
                        break;
                    }
                }
            }
            throw std::runtime_error("Can't link faces");
        }
         */
#endif
    }

    return parent;
}

/// @brief Добавляет в массив сиблингов главного ребенка, сортирует сиблингов
/// по номеру в родительской ячейке
/// @param cells Хранилище ячеек
/// @param ic Номер главного ребенка в хранилище
/// @return Массив индексов дочерних ячеек
template <int dim>
std::array<int, CpC(dim)> sorted_siblings(Storage& cells,  int ic) {
}

/// @brief Основная функция огрубления ячеек, состоит из сбора сиблингов в
/// один массив, вызова функции get_parent и переноса полученных данных в
/// хранилище на место родительской ячейки.
/// @param locals Хранилище ячеек
/// @param ic Индекс родительской ячейки в хранилище
/// @param op Оператор огрубления данных
template<int dim>
void coarse_cell(Storage &locals, Storage &aliens, int rank, int ic, const Distributor& op) {
    locals[ic].set_undefined();

    // главный ребенок всем заведует
    if (locals[ic].z_idx() % CpC(dim) != 0) {
        return;
    }

    auto children = select_children<dim>(locals, ic);

    auto parent = get_parent<dim>(locals, aliens, rank, children);

    int ip = locals[ic].geom().next;

    locals[ic].copy_to(locals[ip]);

    locals[ip].geom() = parent;

    op.merge<dim>(children, locals[ip]);
}

} // namespace amr
} // namespace mesh
} // namespace zephyr