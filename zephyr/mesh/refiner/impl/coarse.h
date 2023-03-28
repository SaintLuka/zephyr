/// @file Файл содержит реализацию одной из базовых функций адаптации coarse_cell,
/// которая объединяет дочерние ячейки в родительскую.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/math/geom/cell.h>
#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/faces.h>
#include <zephyr/mesh/refiner/impl/siblings.h>

namespace zephyr { namespace mesh { namespace impl {

using namespace zephyr::data;
using namespace zephyr::math;
using amrData;

/// @return Массив неопределенных итераторов на дочерние ячейки
template <unsigned int dim>
std::array<Storage::iterator, CpC(dim)> undefined_iterators(Storage& cells);

/// @return Массив неопределенных итераторов на дочерние ячейки (2D)
template <>
std::array<Storage::iterator, CpC(2)> undefined_iterators<2>(Storage& cells) {
    return {cells.end(), cells.end(), cells.end(), cells.end()};
}

/// @return Массив неопределенных итераторов на дочерние ячейки (3D)
template <>
std::array<Storage::iterator, CpC(3)> undefined_iterators<3>(Storage& cells) {
    return {cells.end(), cells.end(), cells.end(), cells.end(),
            cells.end(), cells.end(), cells.end(), cells.end()};
}

/// @brief Возвращает итераторы дочерних ячеек по их индексам в хранилище
/// @param cells Хранилище ячеек
/// @param ic Индекс главной из дочерних ячеек (z_loc = 0)
/// @return Массив итераторов дочерних ячеек
template <unsigned int dim>
std::array<Storage::iterator, CpC(dim)> select_children(Storage& cells, size_t ic) {
    auto sibs = get_siblings<dim>(cells, ic);

    // Дочерние ячейки, упорядоченные по локальному z-индексу
    auto children = undefined_iterators<dim>(cells);
    children[0] = cells[ic];
    for (auto sib: sibs) {
        scrutiny_check(sib < cells.size(), "CellsAround: Sibling index out of range")

        auto z_loc = cells[sib][amrData].z % CpC(dim);
        children[z_loc] = cells[sib];
    }

#if SCRUTINY
    if (cells[ic][amrData].z % CpC(dim) != 0) {
        throw std::runtime_error("Not main child collect siblings");
    }
    auto main_lvl = cells[ic][amrData].level;

    std::set<int> found;
    found.insert(0);
    for (auto sib: sibs) {
        if (cells[sib][amrData].level != main_lvl) {
            throw std::runtime_error("Different levels (siblings)");
        }
        int loc_z = cells[sib][amrData].z % CpC(dim);
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

/// @brief Создает родительскую ячейку с учетом окружения, то есть полученная
/// ячейка будет иметь актуальные ссылки на (старых) соседей, также полученная
/// ячейка может иметь более одной грани по каждой из сторон
/// @param locals Локальное хранилище ячеек
/// @param rank Ранг текущего процесса
/// @param children Массив дочерних ячеек
/// @return Полностью готовая родительская ячейка
template<unsigned int dim>
geom::Cell get_parent(Storage &locals, Storage &aliens, unsigned int rank,
                     std::array<Storage::iterator, CpC(dim)> children) {
    geom::Cell parent(children);

    parent.amrData.base_id = children[0][amrData].base_id;
    parent.amrData.flag = 0;
    parent.amrData.level = children[0][amrData].level - 1;
    parent.amrData.z = children[0][amrData].z / CpC(dim);

    auto children_by_side = get_children_by_side<dim>();

    /// Линкуем грани
    for (unsigned int side = 0; side < FpC(dim); ++side) {
        auto some_child = children[children_by_side[side][0]];
        auto &some_face = some_child[faces].list[side];

#if SCRUTINY
        if (some_face.boundary == FaceFlag::UNDEFINED) {
            throw std::runtime_error("Undefined boundary (coarse cell");
        }
        for (unsigned int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto &face = child[faces].list[side];

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
            parent.faces[side].adjacent.ghost = std::numeric_limits<unsigned int>::max();
            continue;
        }

        auto some_neib_rank = some_face.adjacent.rank;
        auto some_neib_index = some_face.adjacent.index;
        auto some_neib_ghost = some_face.adjacent.ghost;
#if SCRUTINY
        if (some_neib_rank == rank && some_neib_index >= locals.size()) {
            std::cout << "Child has no local neighbor through the " <<
                      side_to_string(side(side % 6)) << " side #1\n";
            print_cell_info(some_child);
            throw std::runtime_error("Child has no local neighbor (coarse_cell) #1");
        }
        if (some_neib_rank != rank && some_neib_ghost >= aliens.size()) {
            std::cout << "Child has no remote neighbor through the " <<
                      side_to_string(side(side % 6)) << " side #1\n";
            print_cell_info(some_child);
            throw std::runtime_error("Child has no remote neighbor (coarse_cell) #1");
        }
#endif

        auto some_neib = some_neib_rank == rank ? locals[some_neib_index] : aliens[some_neib_ghost];
        auto some_neib_wanted_lvl = some_neib[amrData].level + some_neib[amrData].flag;

#if SCRUTINY
        for (unsigned int i = 1; i < VpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto adj = child[faces].list[side].adjacent;

            if (adj.rank == rank && adj.index >= locals.size()) {
                std::cout << "Child has no local neighbor through the " <<
                          side_to_string(side(side % 6)) << " side #2\n";
                print_cell_info(child);
                throw std::runtime_error("Child has no local neighbor (coarse_cell) #2");
            }
            if (adj.rank != rank && adj.ghost >= aliens.size()) {
                std::cout << "Child has no remote neighbor through the " <<
                          side_to_string(side(side % 6)) << " side #2\n";
                print_cell_info(child);
                throw std::runtime_error("Child has no remote neighbor (coarse_cell) #2");
            }

            auto neib = adj.rank == rank ? locals[adj.index] : aliens[adj.ghost];

            auto neib_wanted_lvl = neib[amrData].level + neib[amrData].flag;
            if (neib_wanted_lvl != some_neib_wanted_lvl) {
                throw std::runtime_error("Different wanted level (coarsing cell neighbor)");
            }
        }
#endif

        // Простой случай: грань родительской ячейки должна быть простой
        if (some_neib_wanted_lvl <= parent.amrData.level) {
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
        for (unsigned int i = 0; i < FpF(dim); ++i) {
            _face_& face = parent.faces[side + 6*i];
            pfaces[i] = face_center<dim>(face, parent.vertices);
        }

        // Центры граней дочерних ячеек
        std::array<Vector3d, FpF(dim)> cfaces;
        for (unsigned int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            if (child[faces].list[side + 6].is_undefined()) {
                // Простая грань
                cfaces[i] = face_center<dim>(child[faces].list[side], child[vertices]);
            } else {
                // Сложная грань, считаем центр по основным вершинам
                cfaces[i] = Vector3d(0.0, 0.0, 0.0);
                for (unsigned int j = 0; j < FpF(dim); ++j) {
                    //cfaces[i] += face_center<dim>(child[faces].list[side + 6 * j], child[vertices]);
                    cfaces[i] += (Vector3d&) child[vertices].list[child[faces].list[side + 6 *j].vertices[j]];
                }
                cfaces[i] /= FpF(dim);

            }
        }

        /// Находим соответстиве между pfaces и cfaces
        double eps = 1.0e-6 * parent.size;
        for (unsigned int i = 0; i < FpF(dim); ++i) {
            auto child = children[children_by_side[side][i]];
            auto& child_face = child[faces].list[side];
            for (unsigned int j = 0; j < FpF(dim); ++j) {
                if (distance(pfaces[i], cfaces[j]) < eps) {
                    parent.faces[side + 6*j].adjacent.rank = child_face.adjacent.rank;
                    parent.faces[side + 6*j].adjacent.index = child_face.adjacent.index;
                    parent.faces[side + 6*j].adjacent.ghost = child_face.adjacent.ghost;
                    break;
                }
            }
        }

#if SCRUTINY
        std::set<unsigned int> found;
        for (unsigned int i = 0; i < FpF(dim); ++i) {
            for (unsigned int j = 0; j < FpF(dim); ++j) {
                if (distance(pfaces[i], cfaces[j]) < eps) {
                    found.insert(j);
                    break;
                }
            }
        }
        if (found.size() != FpF(dim)) {
            std::cout << "side: " << side_to_string(side) << "\n";
            std::cout << std::scientific << std::setprecision(6) << "\n";
            std::cout << parent2.vertices.list[parent2.faces[side].vertices[0]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[1]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[2]] << ", "
                      << parent2.vertices.list[parent2.faces[side].vertices[3]] << "\nchildrend:\n";

            for (unsigned int ck = 0; ck < 8; ++ck) {
                auto child = children[ck];
                std::cout << child[vertices].list[child[faces].list[side].vertices[0]] << ", "
                          << child[vertices].list[child[faces].list[side].vertices[1]] << ", "
                          << child[vertices].list[child[faces].list[side].vertices[2]] << ", "
                          << child[vertices].list[child[faces].list[side].vertices[3]] << "\n";
            }

            for (unsigned int i = 0; i < FpF(dim); ++i) {
                std::cout << "pf: " << pfaces[i] << "\n";
                for (unsigned int j = 0; j < FpF(dim); ++j) {
                    std::cout << "  cf: " << cfaces[j] << "\n";
                    if (distance(pfaces[i], cfaces[j]) < eps) {
                        found.insert(j);
                        break;
                    }
                }
            }
            throw std::runtime_error("Can't link faces");
        }
#endif
    }

    return parent;
}

/// @brief Добавляет в массив сиблингов главного ребенка, сортирует сиблингов
/// по номеру в родительской ячейке
/// @param cells Хранилище ячеек
/// @param ic Номер главного ребенка в хранилище
/// @return Массив индексов дочерних ячеек
template <unsigned int dim>
std::array<size_t, CpC(dim)> sorted_siblings(Storage& cells,  size_t ic) {
}

/// @brief Основная функция огрубления ячеек, состоит из сбора сиблингов в
/// один массив, вызова функции get_parent и переноса полученных данных в
/// хранилище на место родительской ячейки.
/// @param locals Хранилище ячеек
/// @param ic Индекс родительской ячейки в хранилище
/// @param op Оператор огрубления данных
template<unsigned int dim>
void coarse_cell(Storage &locals, Storage &aliens, unsigned int rank, size_t ic, const DataDistributor& op) {
    locals[ic][element].set_undefined();

    // главный ребенок всем заведует
    if (locals[ic][amrData].z % CpC(dim) != 0) {
        return;
    }

    auto children = select_children<dim>(locals, ic);

    auto parent = get_parent<dim>(locals, aliens, rank, children);

    size_t ip = locals[ic][amrData].next;

    copy_data(locals[ic], locals[ip]);

    locals[ip][coords] = parent.coords;
    locals[ip][vertices] = parent.vertices;
    locals[ip][faces] = parent.faces;
    locals[ip][size] = parent.size;
    locals[ip][amrData] = parent.amrData;
    locals[ip][element].kind = kind::EULER;
    locals[ip][element].dimension = dim;
    locals[ip][element].rank = rank;
    locals[ip][element].index = std::numeric_limits<size_t>::max();

    op.merge<dim>(children, locals[ip]);
}

} // namespace impl
} // namespace mesh
} // namespace zephyr