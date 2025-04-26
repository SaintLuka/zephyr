/// @file Файл содержит реализацию одной из базовых функций адаптации coarse_cell,
/// которая объединяет дочерние ячейки в родительскую.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/faces.h>
#include <zephyr/mesh/amr2/siblings.h>

namespace zephyr::mesh::amr2 {

template<int dim>
using Children2 = std::array<index_t, CpC(dim)>;

/// @brief Возвращает итераторы дочерних ячеек по их индексам в хранилище
/// @param cells Хранилище ячеек
/// @param ic Индекс главной из дочерних ячеек (z_loc = 0)
/// @return Массив итераторов дочерних ячеек
template <int dim>
Children2<dim> select_children(SoaCell& cells, int ic) {
    auto sibs = get_siblings<dim>(cells, ic);

    // Дочерние ячейки, упорядоченные по локальному z-индексу
    Children2<dim> children;
    children[0] = ic;
    for (auto sib: sibs) {
        scrutiny_check(sib < cells.n_cells(), "Vicinity: Sibling index out of range")

        auto z_loc = cells.z_idx[sib] % CpC(dim);
        children[z_loc] = sib;
    }

#if SCRUTINY
    if (cells.z_idx[ic] % CpC(dim) != 0) {
        throw std::runtime_error("Not main child collect siblings");
    }
    auto main_lvl = cells.level[ic];

    std::set<int> found;
    found.insert(0);
    for (auto sib: sibs) {
        if (cells.level[sib] != main_lvl) {
            throw std::runtime_error("Different levels (siblings)");
        }
        int loc_z = cells.z_idx[sib] % CpC(dim);
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

/// @brief Вершины родительской ячейки (2D)
template <int dim>
std::enable_if_t<dim == 2, SqQuad>
parent_vs(SoaCell& cells, Children2<dim>& children) {
    return {
        cells.verts[cells.node_begin[children[0]] + SqQuad::iss<-1, -1>()],
        cells.verts[cells.node_begin[children[0]] + SqQuad::iss<+1, -1>()],
        cells.verts[cells.node_begin[children[1]] + SqQuad::iss<+1, -1>()],
        cells.verts[cells.node_begin[children[0]] + SqQuad::iss<-1, +1>()],
        cells.verts[cells.node_begin[children[0]] + SqQuad::iss<+1, +1>()],
        cells.verts[cells.node_begin[children[1]] + SqQuad::iss<+1, +1>()],
        cells.verts[cells.node_begin[children[2]] + SqQuad::iss<-1, +1>()],
        cells.verts[cells.node_begin[children[2]] + SqQuad::iss<+1, +1>()],
        cells.verts[cells.node_begin[children[3]] + SqQuad::iss<+1, +1>()]
    };
}

/// @brief Вершины родительской ячейки (3D)
template <int dim>
std::enable_if_t<dim == 3, Cube>
parent_vs(SoaCell& cells, Children2<dim>& children) {
#define subtr3D(i, j, k) (cells.verts[cells.node_begin[children[Cube::iss<i, j, k>()]] + SqCube::iss<i, j, k>()])
    return {
        subtr3D(-1, -1, -1),
        subtr3D(+1, -1, -1),
        subtr3D(-1, +1, -1),
        subtr3D(+1, +1, -1),
        subtr3D(-1, -1, +1),
        subtr3D(+1, -1, +1),
        subtr3D(-1, +1, +1),
        subtr3D(+1, +1, +1)
    };
}

/// @brief Создает родительскую ячейку с учетом окружения, то есть полученная
/// ячейка будет иметь актуальные ссылки на (старых) соседей, также полученная
/// ячейка может иметь более одной грани по каждой из сторон
/// @param locals Локальное хранилище ячеек
/// @param rank Ранг текущего процесса
/// @param children Массив дочерних ячеек
/// @return Полностью готовая родительская ячейка
template<int dim, bool axial=false>
void make_parent(SoaCell& cells, int rank, Children2<dim>& children, index_t ip) {
    const auto children_by_side = get_children_by_side<dim>();

    if constexpr (dim == 2 && axial) {
        cells.add_cell(ip, parent_vs<dim>(cells, children), axial);
    } else {
        cells.add_cell(ip, parent_vs<dim>(cells, children));
    }

    auto& adj = cells.faces.adjacent;

    cells.rank[ip]        = cells.rank[children[0]];
    cells.owner_index[ip] = cells.next[children[0]];
    cells.next[ip]        = cells.next[children[0]];

    cells.b_idx[ip] = cells.b_idx[children[0]];
    cells.flag[ip]  = 0;
    cells.level[ip] = cells.level[children[0]] - 1;
    cells.z_idx[ip] = cells.z_idx[children[0]] / CpC(dim);

    /// Линкуем грани
    for (int side = 0; side < FpC(dim); ++side) {
        index_t pface = cells.face_begin[ip] + side;

        index_t some_ch = children[children_by_side[side][0]];
        index_t some_face = cells.face_begin[some_ch] + side;


#if SCRUTINY
        if (cells.faces.boundary[some_face] == Boundary::UNDEFINED) {
            throw std::runtime_error("Undefined boundary (coarse cell");
        }
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t iface = cells.face_begin[ich] + side;

            if (cells.faces.boundary[iface] != cells.faces.boundary[some_face]) {
                throw std::runtime_error("Different boundary conditions");
            }
        }
#endif
        cells.faces.boundary[pface + side] = cells.faces.boundary[some_face];

        // Внешняя граница, не требуется линковать
        if (cells.faces.is_boundary(some_face)) {
            adj.rank[pface + side] = rank;
            continue;
        }

        auto some_neib_rank  = adj.rank[some_face];
        auto some_neib_owner = adj.owner_index[some_face];
        auto some_neib_local = adj.local_index[some_face];
#if SCRUTINY
        if (some_neib_rank == rank && some_neib_owner >= cells.n_cells()) {
            std::cout << "Child has no local neighbor through the " <<
                      side_to_string(side % 6) << " side #1\n";
            cells.print_info(some_ch);
            throw std::runtime_error("Child has no local neighbor (coarse_cell) #1");
        }
        if (some_neib_rank != rank && some_neib_local >= cells.n_cells()) {
            std::cout << "Child has no remote neighbor through the " <<
                      side_to_string(side % 6) << " side #1\n";
            cells.print_info(some_ch);
            throw std::runtime_error("Child has no remote neighbor (coarse_cell) #1");
        }
#endif

        index_t some_neib = adj.local_index[some_face];
        auto some_neib_wanted_lvl = cells.level[some_neib] + cells.flag[some_neib];

#if SCRUTINY
        for (int i = 1; i < VpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t iface = cells.face_begin[ich] + side;

            if (adj.rank[iface] == rank && adj.local_index[iface] >= cells.n_cells()) {
                std::cout << "Child has no local neighbor through the " <<
                          side_to_string(side % 6) << " side #2\n";
                cells.print_info(ich);
                throw std::runtime_error("Child has no local neighbor (coarse_cell) #2");
            }
            if (adj.rank[iface] != rank && adj.local_index[iface] >= cells.n_cells()) {
                std::cout << "Child has no remote neighbor through the " <<
                          side_to_string(side % 6) << " side #2\n";
                cells.print_info(ich);
                throw std::runtime_error("Child has no remote neighbor (coarse_cell) #2");
            }

            index_t neib = adj.local_index[iface];

            auto neib_wanted_lvl = cells.level[neib] + cells.flag[neib];
            if (neib_wanted_lvl != some_neib_wanted_lvl) {
                throw std::runtime_error("Different wanted level (coarsing cell neighbor)");
            }
        }
#endif

        // Простой случай: грань родительской ячейки должна быть простой
        if (some_neib_wanted_lvl <= cells.level[ip]) {
            adj.rank[pface + side] = rank;
            adj.owner_index[pface + side] = some_neib_owner;
            adj.local_index[pface + side] = some_neib_local;

            // Обнуляем неактивные подграни
            for (int i = 1; i < FpF(dim); ++i) {
                cells.faces.set_undefined(pface + side + 6 * i);
            }
            continue;
        }

        // Если мы здесь, то сторона side от родителя должна адаптироваться
        // Далее потребуется связать грани

        split_face<dim>(pface, cells.faces, cells.get_vertices<dim>(ip), side, axial);

        // Центры подграней родительской ячейки
        std::array<Vector3d, FpF(dim)> pfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            pfaces[i] = cells.faces.center[pface + side + 6 * i];
        }

        // Центры граней дочерних ячеек
        std::array<Vector3d, FpF(dim)> cfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t iface = cells.face_begin[ich] + side;
            if (cells.faces.is_undefined(iface + side + 6)) {
                // Простая грань
                cfaces[i] = cells.faces.center[iface];
            } else {
                // Сложная грань, считаем центр тяжести
                cfaces[i] = Vector3d(0.0, 0.0, 0.0);
                double full_area = 0.0;
                for (int j = 0; j < FpF(dim); ++j) {
                    double area = cells.faces.get_area(iface + side + 6 * j, axial);
                    full_area += area;
                    cfaces[i] += area * cells.faces.center[side + 6 * j];
                }
                cfaces[i] /= full_area;

            }
        }

        /// Находим соответстиве между pfaces и cfaces
        double eps = 1.0e-3 * cells.linear_size(ip);
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t ch_face = cells.face_begin[ich] + side;
            for (int j = 0; j < FpF(dim); ++j) {
                if ((pfaces[i] - cfaces[j]).norm() < eps) {
                    adj.rank[pface + side + 6*j]        = adj.rank[ch_face];
                    adj.local_index[pface + side + 6*j] = adj.owner_index[ch_face];
                    adj.owner_index[pface + side + 6*j] = adj.local_index[ch_face];
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

        if (found.size() != FpF(dim)) {
            std::cout << "Side: " << side_to_string(side) << "\n";
            std::cout << std::scientific << std::setprecision(6) << "\n";
            std::cout << "Parent:\n";
            /*
            std::cout << parent.vertices[parent.faces[side].vertices[0]].transpose() << ", "
                      << parent.vertices[parent.faces[side].vertices[1]].transpose() << ", "
                      << parent.vertices[parent.faces[side].vertices[2]].transpose() << ", "
                      << parent.vertices[parent.faces[side].vertices[3]].transpose() << "\n";

            std::cout << "Children:\n";
            for (int ck = 0; ck < CpC(dim); ++ck) {
                auto& child = children[ck];
                std::cout << child.vertices[child.faces[side].vertices[0]].transpose() << ", "
                          << child.vertices[child.faces[side].vertices[1]].transpose() << ", "
                          << child.vertices[child.faces[side].vertices[2]].transpose() << ", "
                          << child.vertices[child.faces[side].vertices[3]].transpose() << "\n";
            }
            */

            for (int i = 0; i < FpF(dim); ++i) {
                std::cout << "pface: " << pfaces[i].transpose() << "\n";
                for (int j = 0; j < FpF(dim); ++j) {
                    std::cout << "  cface: " << cfaces[j].transpose() << "\n";
                    std::cout << "  norm:  " << (pfaces[i] - cfaces[j]).norm() << "\n";
                    if ((pfaces[i] - cfaces[j]).norm() < eps) {
                        found.insert(j);
                        break;
                    }
                }
            }

            cells.visualize(ip, "parent.py");
            for (auto& child: children) {
                cells.visualize(child, "child" + std::to_string(cells.z_idx[child]) + ".py");
            }

            throw std::runtime_error("Can't link faces");
        }
#endif
    }
}

/// @brief Основная функция огрубления ячеек, состоит из сбора сиблингов в
/// один массив, вызова функции get_parent и переноса полученных данных в
/// хранилище на место родительской ячейки.
/// @param locals Хранилище ячеек
/// @param op Оператор огрубления данных
template<int dim>
void coarse_cell(index_t ich, SoaCell& cells, int rank, const Distributor& op) {
    // Функцию выполняет главный ребенок, остальные
    // выставляются на undefined и отдыхают
    if (cells.z_idx[ich] % CpC(dim) != 0) {
        cells.set_undefined(ich);
        return;
    }

    auto children = select_children<dim>(cells, ich);

    index_t ip = cells.next[ich];
    if (dim == 2 && cells.axial) {
        make_parent<dim, true>(cells, rank, children, ip);
    }
    else {
        make_parent<dim>(cells, rank, children, ip);
    }
    //op.merge(children, parent);

    cells.set_undefined(ich);
}

} // namespace zephyr::mesh::amr2