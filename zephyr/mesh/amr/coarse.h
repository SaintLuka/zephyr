// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>
#include <zephyr/mesh/amr/siblings.h>

namespace zephyr::mesh::amr {

/// @brief Возвращает итераторы дочерних ячеек по их индексам в хранилище
/// @param locals Хранилище ячеек
/// @param ic Индекс главной из дочерних ячеек (z_loc = 0)
/// @return Массив итераторов дочерних ячеек
template <int dim>
Children select_children(AmrCells& locals, index_t ic) {
    auto sibs = get_siblings<dim>(locals, ic);

    // Дочерние ячейки, упорядоченные по локальному z-индексу
    Children children(&locals);
    children.index[0] = ic;
    for (auto sib: sibs) {
        scrutiny_check(sib < locals.size(), "Vicinity: Sibling index out of range")

        auto z_loc = locals.z_idx[sib] % CpC(dim);
        children.index[z_loc] = sib;
    }

#if SCRUTINY
    if (locals.z_idx[ic] % CpC(dim) != 0) {
        throw std::runtime_error("Not main child collect siblings");
    }
    auto main_lvl = locals.level[ic];

    std::set<int> found;
    found.insert(0);
    for (auto sib: sibs) {
        if (locals.level[sib] != main_lvl) {
            throw std::runtime_error("Different levels (siblings)");
        }
        int loc_z = locals.z_idx[sib] % CpC(dim);
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

/// @brief Вершины родительской ячейки (2D или 3D)
template <int dim>
SqMap<dim> parent_vs(const AmrCells& cells, Children& children) {
#define subs_vertex_3D(i, j, k) (cells.verts[cells.node_begin[children.index[Cube::iss<i, j, k>()]] + SqCube::iss<i, j, k>()])
    if constexpr (dim == 2) {
        return {
            cells.verts[cells.node_begin[children.index[0]] + SqQuad::iss<-1, -1>()],
            cells.verts[cells.node_begin[children.index[0]] + SqQuad::iss<+1, -1>()],
            cells.verts[cells.node_begin[children.index[1]] + SqQuad::iss<+1, -1>()],
            cells.verts[cells.node_begin[children.index[0]] + SqQuad::iss<-1, +1>()],
            cells.verts[cells.node_begin[children.index[0]] + SqQuad::iss<+1, +1>()],
            cells.verts[cells.node_begin[children.index[1]] + SqQuad::iss<+1, +1>()],
            cells.verts[cells.node_begin[children.index[2]] + SqQuad::iss<-1, +1>()],
            cells.verts[cells.node_begin[children.index[2]] + SqQuad::iss<+1, +1>()],
            cells.verts[cells.node_begin[children.index[3]] + SqQuad::iss<+1, +1>()]
        };
    }
    else {
        return {
            subs_vertex_3D(-1, -1, -1),
            subs_vertex_3D(+1, -1, -1),
            subs_vertex_3D(-1, +1, -1),
            subs_vertex_3D(+1, +1, -1),
            subs_vertex_3D(-1, -1, +1),
            subs_vertex_3D(+1, -1, +1),
            subs_vertex_3D(-1, +1, +1),
            subs_vertex_3D(+1, +1, +1)
        };
    }
}

/// @brief Создает родительскую ячейку. Ячейка ссылается на соседей после
/// перемещения, то есть учитывает проставленные индексы next.
/// @param locals Локальное хранилище ячеек
/// @param aliens Хранилище ячеек с другого процесса
/// @param children Массив дочерних ячеек
/// @param ip Индекс в хранилище locals, по которому следует разместить родительскую ячейку
/// @param rank Ранг текущего процесса
template<int dim>
void make_parent(AmrCells& locals, AmrCells& aliens, Children& children, index_t ip, int rank) {
    if constexpr (dim == 2) {
        locals.set_cell(ip, parent_vs<dim>(locals, children), locals.axial());
    }
    else {
        locals.set_cell(ip, parent_vs<dim>(locals, children));
    }

    // Выставляем основные свойства
    locals.rank [ip] = locals.rank[children.index[0]];
    locals.next [ip] = ip;
    locals.index[ip] = ip;

    locals.flag [ip] = 0;
    locals.b_idx[ip] = locals.b_idx[children.index[0]];
    locals.level[ip] = locals.level[children.index[0]] - 1;
    locals.z_idx[ip] = locals.z_idx[children.index[0]] / CpC(dim);

    // Связываем грани
    auto& adj = locals.faces.adjacent;
    for (Side<dim> side: Side<dim>::items()) {
        // Некоторая дочерняя у грани и её грань
        index_t some_ch = children.index[side.child()];
        index_t some_ch_face = locals.face_begin[some_ch] + side;
#if SCRUTINY
        if (locals.faces.is_undefined(some_ch_face)) {
            throw std::runtime_error("Undefined boundary (coarse cell");
        }
        for (auto subface: side.subfaces()) {
            index_t ich = children.index[subface.child()];
            index_t iface = locals.face_begin[ich] + side;

            if (locals.faces.boundary[iface] != locals.faces.boundary[some_ch_face]) {
                throw std::runtime_error("Different boundary conditions");
            }
        }
#endif
        index_t face_beg = locals.face_begin[ip];
        locals.faces.boundary[face_beg + side] = locals.faces.boundary[some_ch_face];

        // Внешняя граница, не требуется линковать
        if (locals.faces.is_boundary(some_ch_face)) {
            adj.rank [face_beg + side] = rank;
            adj.index[face_beg + side] = ip;
            adj.alien[face_beg + side] = -1;
            adj.basic[face_beg + side] = ip;
            continue;
        }

        auto some_neib_rank  = adj.rank [some_ch_face];
        auto some_neib_alien = adj.alien[some_ch_face];
        auto some_neib_index = adj.index[some_ch_face];
#if SCRUTINY
        if (some_neib_rank == rank && (some_neib_alien >= 0 || some_neib_index >= locals.size())) {
            std::cout << "Child has no local neighbor through the " <<
                      side_to_string(side, dim) << " side #1\n";
            locals.print_info(some_ch);
            throw std::runtime_error("Child has no local neighbor (coarse_cell) #1");
        }
        if (some_neib_rank != rank && (some_neib_alien < 0 || some_neib_alien >= aliens.size())) {
            std::cout << "Child has no remote neighbor through the " <<
                      side_to_string(side, dim) << " side #1\n";
            locals.print_info(some_ch);
            throw std::runtime_error("Child has no remote neighbor (coarse_cell) #1");
        }
#endif
        auto [some_neibs, some_neib] = adj.get_neib(some_ch_face, locals, aliens);
        auto some_neib_wanted_lvl = some_neibs.level[some_neib] + some_neibs.flag[some_neib];
#if SCRUTINY
        auto children_by_side = side.children();
        for (int i = 1; i < VpF(dim); ++i) {
            index_t ich = children.index[children_by_side[i]];
            index_t iface = locals.face_begin[ich] + side;

            if (adj.rank[iface] == rank && (adj.alien[iface] >= 0 ||
                    adj.index[iface] < 0 || adj.index[iface] >= locals.size())) {
                std::cout << "Child has no local neighbor through the " <<
                          side_to_string(side, dim) << " side #2\n";
                locals.print_info(ich);
                throw std::runtime_error("Child has no local neighbor (coarse_cell) #2");
            }
            if (adj.rank[iface] != rank && (adj.alien[iface] < 0 || adj.alien[iface] >= aliens.size())) {
                std::cout << "Child has no remote neighbor through the " <<
                          side_to_string(side, dim) << " side #2\n";
                locals.print_info(ich);
                throw std::runtime_error("Child has no remote neighbor (coarse_cell) #2");
            }

            auto [neibs2, neib2] = adj.get_neib(iface, locals, aliens);

            auto neib_wanted_lvl = neibs2.level[neib2] + neibs2.flag[neib2];
            if (neib_wanted_lvl != some_neib_wanted_lvl) {
                throw std::runtime_error("Different wanted level (coarsing cell neighbor)");
            }
        }
#endif
        // Простой случай: грань родительской ячейки должна быть простой
        if (some_neib_wanted_lvl <= locals.level[ip]) {
            // TODO: MPI VERSION
            adj.rank [face_beg + side] = some_neib_rank;
            adj.alien[face_beg + side] = some_neib_alien;
            adj.index[face_beg + side] = some_neibs.next[some_neib];
            adj.basic[face_beg + side] = ip;

            // Обнуляем неактивные подграни
            for (int i = 1; i < FpF(dim); ++i) {
                locals.faces.set_undefined(face_beg + side[i]);
            }
            continue;
        }

        // Если мы здесь, то сторона side от родителя должна адаптироваться
        // Далее потребуется связать грани

        split_face<dim>(locals, ip, side);

        // lvl_c > lvl_n & flag_n == 1
        // Сосед выше и ничего не делает
        // Там есть сосед, который делает coarse
        for (auto subface: side.subfaces()) {
            // TODO: MPI VERSION
            index_t ich = children.index[subface.child()];
            index_t ch_face = locals.face_begin[ich] + side;
            auto[neibs, jc] = adj.get_neib(ch_face, locals, aliens);
            scrutiny_check(0 <= jc && jc < neibs.size(), "Out-of-bounds #3");

            adj.rank [face_beg + subface] = adj.rank [ch_face];
            adj.index[face_beg + subface] = -1;
            adj.alien[face_beg + subface] = adj.alien[ch_face];
            adj.basic[face_beg + subface] = ip;

            if (locals.level[ich] > neibs.level[jc]) {
                scrutiny_check(neibs.flag[jc] > 0, "Wrong assumption #4");
                adj.index[face_beg + subface] = neibs.next[jc] + subface.neib_child();
            }
            else {
                scrutiny_check(neibs.flag[jc] <= 0, "Wrong assumption #5");
                adj.index[face_beg + subface] = neibs.next[jc];
            }

            scrutiny_check(adj.index[face_beg + subface] >= 0, "Not all cases");
        }

#if SCRUTINY
        // Центры подграней родительской ячейки
        std::array<Vector3d, FpF(dim)> pfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            pfaces[i] = locals.faces.center[face_beg + side[i]];
        }

        // Центры граней дочерних ячеек
        std::array<Vector3d, FpF(dim)> cfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children.index[children_by_side[i]];
            index_t iface = locals.face_begin[ich];
            if (locals.simple_face(ich, side)) {
                // Простая грань
                cfaces[i] = locals.faces.center[iface + side];
            } else {
                // Сложная грань, считаем центр тяжести
                cfaces[i] = Vector3d::Zero();
                double full_area = 0.0;
                for (int j = 0; j < FpF(dim); ++j) {
                    double area = locals.faces.get_area(iface + side[j], locals.axial());
                    full_area += area;
                    cfaces[i] += area * locals.faces.center[iface + side[j]];
                }
                cfaces[i] /= full_area;
            }
        }

        /// Находим соответстиве между pfaces и cfaces
        double eps = 1.0e-3 * locals.linear_size(ip);
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children.index[children_by_side[i]];
            index_t ch_face = locals.face_begin[ich] + side;
            for (int j = 0; j < FpF(dim); ++j) {
                if ((pfaces[i] - cfaces[j]).norm() < eps) {
                    break;
                }
            }
        }

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
            std::cout << "Side: " << side_to_string(side, dim) << "\n";
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
                //std::cout << "pface: " << pfaces[i].transpose() << "\n";
                for (int j = 0; j < FpF(dim); ++j) {
                    //std::cout << "  cface: " << cfaces[j].transpose() << "\n";
                    //std::cout << "  norm:  " << (pfaces[i] - cfaces[j]).norm() << "\n";
                    if ((pfaces[i] - cfaces[j]).norm() < eps) {
                        found.insert(j);
                        break;
                    }
                }
            }

            locals.visualize(ip, "parent.py");
            for (auto child: children.index) {
                if (child < 0) continue;
                locals.visualize(child, "child" + std::to_string(locals.z_idx[child]) + ".py");
            }

            locals.print_info(ip);

            throw std::runtime_error("Can't link faces");
        }
#endif
    }
}

/// @brief Основная функция огрубления ячеек, состоит из сбора сиблингов в
/// один массив, вызова функции make_parent и переноса полученных данных в
/// хранилище на место родительской ячейки.
/// @param locals Хранилище ячеек
/// @param aliens Хранилище ячеек с других процессов
/// @param ich Индекс дочерней ячейки, для которой выполняется огрубление
/// @param op Оператор огрубления данных
/// @param rank Ранг текущего процесса
template<int dim>
void coarse_cell(AmrCells& locals, AmrCells& aliens, index_t ich, const Distributor& op, int rank) {
    // Функцию выполняет главный ребенок, остальные выставляются на undefined и отдыхают
    if (locals.z_idx[ich] % CpC(dim) != 0) {
        locals.set_undefined(ich);
        return;
    }

    auto children = select_children<dim>(locals, ich);

    index_t ip = locals.next[ich];
    make_parent<dim>(locals, aliens, children, ip, rank);

    EuCell parent(&locals, ip);
    op.merge(children, parent);

    locals.set_undefined(ich);
}

} // namespace zephyr::mesh::amr