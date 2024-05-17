#pragma once

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/generator.h>

#include <zephyr/geom/primitives/mov_node.h>
#include <zephyr/geom/primitives/mov_cell.h>

#include <zephyr/mesh/storage.h>


namespace zephyr::mesh {

using namespace zephyr::geom;

/// @brief Хранилище с подвижными ячейками
using CellStorage = Storage<MovCell>;

/// @brief Хранилище с узлами подвижной сетки
using NodeStorage = Storage<MovNode>;

class LaMesh {
public:

    template<class U, class V>
    LaMesh(const U &u, const V &v, Generator *gen)
            : m_cells(u, 0), m_nodes(v, 0) {
        initialize(gen->make());
    }

    template <class U, class V>
    LaMesh(const U &u, const V &v, const Grid& grid)
            : m_cells(u), m_nodes(v) {
        initialize(grid);
    }

    template<class U>
    LaMesh(const U &u, Generator *gen)
            : m_cells(u), m_nodes() {
        initialize(gen->make());
    }

    template <class U>
    LaMesh(const U &u, const Grid& grid)
            : m_cells(u), m_nodes() {
        initialize(grid);
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