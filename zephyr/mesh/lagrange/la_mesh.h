#pragma once

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/generator.h>

#include <zephyr/mesh/primitives/mov_node.h>
#include <zephyr/mesh/primitives/mov_cell.h>

#include <zephyr/mesh/storage.h>


namespace zephyr::mesh {

using namespace zephyr::geom;

/// @brief Хранилище с подвижными ячейками
using CellStorage = Storage<MovCell>;

/// @brief Хранилище с узлами подвижной сетки
using NodeStorage = Storage<MovNode>;

class LaMesh {
public:

    template <class U, class V>
    LaMesh(const Grid& grid, const U &u, const V &v) : m_cells(u), m_nodes(v) {
        initialize(grid);
    }

    template<class U, class V>
    LaMesh(Generator *gen, const U &u, const V &v) : m_cells(0, u), m_nodes(0, v) {
        initialize(gen->make());
    }

    template <class U>
    LaMesh(const Grid& grid, const U &u) : m_cells(0, u), m_nodes(0, true, 0) {
        initialize(grid);
    }

    template<class U>
    LaMesh(Generator *gen, const U &u) : m_cells(0, u), m_nodes(0, true, 0) {
        initialize(gen->make());
    }


    CellStorage& locals() {
        return m_cells;
    }

    NodeStorage& nodes() {
        return m_nodes;
    }

    geom::Box bbox();

private:
    void initialize(const Grid& gen);

    CellStorage m_cells;
    NodeStorage m_nodes;
};


} // namespace zephyr::mesh