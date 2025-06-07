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


/// @brief Возвращает итераторы дочерних ячеек по их индексам в хранилище
/// @param locals Хранилище ячеек
/// @param ic Индекс главной из дочерних ячеек (z_loc = 0)
/// @return Массив итераторов дочерних ячеек
template <int dim>
SoaChildren select_children(AmrCells& locals, int ic) {
    auto sibs = get_siblings<dim>(locals, ic);

    // Дочерние ячейки, упорядоченные по локальному z-индексу
    SoaChildren children(&locals);
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

/// @brief Вершины родительской ячейки (2D)
template <int dim>
std::enable_if_t<dim == 2, SqQuad>
parent_vs(AmrCells& cells, SoaChildren& children) {
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

/// @brief Вершины родительской ячейки (3D)
template <int dim>
std::enable_if_t<dim == 3, Cube>
parent_vs(AmrCells& cells, SoaChildren& children) {
#define subtr3D(i, j, k) (cells.verts[cells.node_begin[children.index[Cube::iss<i, j, k>()]] + SqCube::iss<i, j, k>()])
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
void make_parent(AmrCells& locals, AmrCells& aliens, int rank, SoaChildren& children, index_t ip) {
    const auto children_by_side = get_children_by_side<dim>();

    if constexpr (dim == 2 && axial) {
        locals.set_cell(ip, parent_vs<dim>(locals, children), axial);
    } else {
        locals.set_cell(ip, parent_vs<dim>(locals, children));
    }

    auto& adj = locals.faces.adjacent;

    locals.rank[ip]  = locals.rank[children.index[0]];
    locals.index[ip] = locals.next[children.index[0]];
    locals.next[ip]  = locals.next[children.index[0]];

    locals.b_idx[ip] = locals.b_idx[children.index[0]];
    locals.flag[ip]  = 0;
    locals.level[ip] = locals.level[children.index[0]] - 1;
    locals.z_idx[ip] = locals.z_idx[children.index[0]] / CpC(dim);

    /// Линкуем грани
    for (Side<dim> side: Side<dim>::items()) {
        index_t pface = locals.face_begin[ip];

        index_t some_ch = children.index[children_by_side[side][0]];
        index_t some_face = locals.face_begin[some_ch] + side;

#if SCRUTINY
        if (locals.faces.boundary[some_face] == Boundary::UNDEFINED) {
            throw std::runtime_error("Undefined boundary (coarse cell");
        }
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t iface = locals.face_begin[ich] + side;

            if (locals.faces.boundary[iface] != locals.faces.boundary[some_face]) {
                throw std::runtime_error("Different boundary conditions");
            }
        }
#endif
        locals.faces.boundary[pface + side] = locals.faces.boundary[some_face];

        // Внешняя граница, не требуется линковать
        if (locals.faces.is_boundary(some_face)) {
            adj.rank[pface + side] = rank;
            continue;
        }

        auto some_neib_rank  = adj.rank[some_face];
        auto some_neib_owner = adj.alien[some_face];
        auto some_neib_local = adj.index[some_face];
#if SCRUTINY
        if (some_neib_rank == rank && some_neib_owner >= locals.size()) {
            std::cout << "Child has no local neighbor through the " <<
                      side_to_string(side, dim) << " side #1\n";
            locals.print_info(some_ch);
            throw std::runtime_error("Child has no local neighbor (coarse_cell) #1");
        }
        if (some_neib_rank != rank && some_neib_local >= locals.size()) {
            std::cout << "Child has no remote neighbor through the " <<
                      side_to_string(side, dim) << " side #1\n";
            locals.print_info(some_ch);
            throw std::runtime_error("Child has no remote neighbor (coarse_cell) #1");
        }
#endif

        auto [neibs, some_neib] = adj.get_neib(some_face, locals, aliens);
        auto some_neib_wanted_lvl = neibs.level[some_neib] + neibs.flag[some_neib];

#if SCRUTINY
        for (int i = 1; i < VpF(dim); ++i) {
            index_t ich = children[children_by_side[side][i]];
            index_t iface = locals.face_begin[ich] + side;

            if (adj.rank[iface] == rank && adj.index[iface] >= locals.size()) {
                std::cout << "Child has no local neighbor through the " <<
                          side_to_string(side, dim) << " side #2\n";
                locals.print_info(ich);
                throw std::runtime_error("Child has no local neighbor (coarse_cell) #2");
            }
            if (adj.rank[iface] != rank && adj.index[iface] >= locals.size()) {
                std::cout << "Child has no remote neighbor through the " <<
                          side_to_string(side, dim) << " side #2\n";
                locals.print_info(ich);
                throw std::runtime_error("Child has no remote neighbor (coarse_cell) #2");
            }

            index_t neib = adj.index[iface];

            auto neib_wanted_lvl = locals.level[neib] + locals.flag[neib];
            if (neib_wanted_lvl != some_neib_wanted_lvl) {
                throw std::runtime_error("Different wanted level (coarsing cell neighbor)");
            }
        }
#endif

        // Простой случай: грань родительской ячейки должна быть простой
        if (some_neib_wanted_lvl <= locals.level[ip]) {
            adj.rank[pface + side] = rank;
            adj.alien[pface + side] = some_neib_owner;
            adj.index[pface + side] = some_neib_local;
            adj.basic[pface + side] = ip;

            // Обнуляем неактивные подграни
            for (int i = 1; i < FpF(dim); ++i) {
                locals.faces.set_undefined(pface + side[i]);
            }
            continue;
        }

        // Если мы здесь, то сторона side от родителя должна адаптироваться
        // Далее потребуется связать грани

        split_face<dim>(pface, locals.faces, locals.mapping<dim>(ip), side, axial);

        // Центры подграней родительской ячейки
        std::array<Vector3d, FpF(dim)> pfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            pfaces[i] = locals.faces.center[pface + side[i]];
        }

        // Центры граней дочерних ячеек
        std::array<Vector3d, FpF(dim)> cfaces;
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children.index[children_by_side[side][i]];
            index_t iface = locals.face_begin[ich];
            if (locals.simple_face(ich, side)) {
                // Простая грань
                cfaces[i] = locals.faces.center[iface + side];
            } else {
                // Сложная грань, считаем центр тяжести
                cfaces[i] = Vector3d::Zero();
                double full_area = 0.0;
                for (int j = 0; j < FpF(dim); ++j) {
                    double area = locals.faces.get_area(iface + side[j], axial);
                    full_area += area;
                    cfaces[i] += area * locals.faces.center[iface + side[j]];
                }
                cfaces[i] /= full_area;

            }
        }

        /// Находим соответстиве между pfaces и cfaces
        double eps = 1.0e-3 * locals.linear_size(ip);
        for (int i = 0; i < FpF(dim); ++i) {
            index_t ich = children.index[children_by_side[side][i]];
            index_t ch_face = locals.face_begin[ich] + side;
            for (int j = 0; j < FpF(dim); ++j) {
                if ((pfaces[i] - cfaces[j]).norm() < eps) {
                    adj.rank [pface + side[j]] = adj.rank[ch_face];
                    adj.index[pface + side[j]] = adj.index[ch_face];
                    adj.alien[pface + side[j]] = adj.alien[ch_face];
                    adj.basic[pface + side[j]] = ip;
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
            for (auto& child: children) {
                locals.visualize(child, "child" + std::to_string(locals.z_idx[child]) + ".py");
            }

            locals.print_info(ip);

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
void coarse_cell(AmrCells& locals, AmrCells& aliens, index_t ich, int rank, const Distributor& op) {
    // Функцию выполняет главный ребенок, остальные
    // выставляются на undefined и отдыхают
    if (locals.z_idx[ich] % CpC(dim) != 0) {
        locals.set_undefined(ich);
        return;
    }

    auto children = select_children<dim>(locals, ich);

    index_t ip = locals.next[ich];
    if (dim == 2 && locals.axial()) {
        make_parent<dim, true>(locals, aliens, rank, children, ip);
    }
    else {
        make_parent<dim>(locals, aliens, rank, children, ip);
    }
    QCell parent(&locals, ip);
    op.merge_soa(children, parent);

    locals.set_undefined(ich);
}

} // namespace zephyr::mesh::amr2