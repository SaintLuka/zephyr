#include <iostream>

#include <zephyr/geom/geom.h>
#include <zephyr/math/funcs.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/amr_cells.h>

using namespace zephyr::geom;

namespace zephyr::mesh {

AmrCells::AmrCells(const generator::Strip& gen) {
    throw std::runtime_error("AmrCells::AmrCells strip gen");
    /*
    bool x_period = periodic_along_x();

    double m_ymin = y_min();
    double m_ymax = y_max();

    index_t m_ny = 1;

    double hx = (m_xmax - m_xmin) / m_nx;
    double hy = (m_ymax - m_ymin) / m_ny;

    auto get_vertex = [=](index_t i, index_t j) -> Vector3d {
        return {
            m_xmin + ((m_xmax - m_xmin) * i) / m_nx,
            m_ymin + ((m_ymax - m_ymin) * j) / m_ny,
            0.0
    };
    };

    auto neib_index = [=](index_t i, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  i : (i - 1 + m_nx) % m_nx;
        }
        else if (side == Side2D::RIGHT) {
            return i == m_nx - 1 && !x_period ? i : (i + 1 + m_nx) % m_nx;
        }
        else {
            throw std::runtime_error("Strange side #153");
        }
    };

    cells.set_dimension(2);
    cells.set_adaptive(true);
    cells.set_linear(true);
    cells.set_axial(false);

    cells.resize_amr(m_nx);

    int n_faces = 8;
    int n_nodes = 9;

    for (index_t ic = 0; ic < m_nx; ++ic) {
        cells.next[ic] = ic;
        cells.rank[ic] = 0;
        cells.index[ic] = ic;

        cells.flag[ic] = 0;
        cells.level[ic] = 0;
        cells.b_idx[ic] = ic;
        cells.z_idx[ic] = 0;

        SqQuad quad(
                get_vertex(ic, 0),
                get_vertex(ic + 1, 0),
                get_vertex(ic, 1),
                get_vertex(ic + 1, 1));

        cells.center[ic] = quad.vs<0, 0>();
        cells.volume[ic] = hx * hy;
        cells.volume_alt[ic] = NAN;
        cells.face_begin[ic] = ic * n_faces;
        cells.node_begin[ic] = ic * n_nodes;
        cells.face_begin[ic + 1] = cells.face_begin[ic] + n_faces;
        cells.node_begin[ic + 1] = cells.face_begin[ic] + n_nodes;

        // INIT FACES
        for (index_t iface: cells.faces_range(ic)) {
            cells.faces.set_undefined(iface);
            cells.faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        cells.faces.boundary[iface + Side2D::L] = ic > 0 ? Boundary::INNER : m_bounds.left;
        cells.faces.boundary[iface + Side2D::R] = ic < m_nx - 1 ? Boundary::INNER : m_bounds.right;

        for (auto side: {Side2D::LEFT, Side2D::RIGHT}) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(ic, side);
            cells.faces.adjacent.alien[iface + side] = -1;
            cells.faces.adjacent.basic[iface + side] = ic;
            cells.faces.vertices[iface + side].fill(-1);
        }

        cells.faces.normal[iface + Side2D::L] = -Vector3d::UnitX();
        cells.faces.normal[iface + Side2D::R] =  Vector3d::UnitX();
        cells.faces.normal[iface + Side2D::B] = -Vector3d::UnitY();
        cells.faces.normal[iface + Side2D::T] =  Vector3d::UnitY();

        cells.faces.center[iface + Side2D::L] = quad.vs<-1, 0>();
        cells.faces.center[iface + Side2D::R] = quad.vs<+1, 0>();
        cells.faces.center[iface + Side2D::B] = quad.vs<0, -1>();
        cells.faces.center[iface + Side2D::T] = quad.vs<0, +1>();

        cells.faces.area[iface + Side2D::L] = hy;
        cells.faces.area[iface + Side2D::R] = hy;
        cells.faces.area[iface + Side2D::B] = hx;
        cells.faces.area[iface + Side2D::T] = hx;

        cells.faces.vertices[iface + Side2D::L] = Side2D::L.sf();
        cells.faces.vertices[iface + Side2D::R] = Side2D::R.sf();
        cells.faces.vertices[iface + Side2D::B] = Side2D::B.sf();
        cells.faces.vertices[iface + Side2D::T] = Side2D::T.sf();

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            cells.verts[ic * n_nodes + jn] = quad[jn];
        }
    }
    */
}

AmrCells::AmrCells(const generator::Rectangle& rect) {
    // ------- Считываем параметры Rectangle -------
    const bool x_period = rect.periodic_along_x();
    const bool y_period = rect.periodic_along_y();
    
    const int nx = rect.nx();
    const int ny = rect.ny();
    
    const double x_min = rect.x_min();
    const double x_max = rect.x_max();
    const double y_min = rect.y_min();
    const double y_max = rect.y_max();

    const double hx = (x_max - x_min) / nx;
    const double hy = (y_max - y_min) / ny;

    const auto [left, right, bottom, top] = rect.bounds();

    // ------- Индексация -------

    auto get_index = [ny]
    (index_t i, index_t j) -> index_t {
        return i * ny + j;
    };

    auto get_index_pair = [ny]
    (index_t n) -> std::array<index_t, 2> {
        return {n / ny, n % ny};
    };

    auto get_vertex = [x_min, x_max, y_min, y_max, nx, ny]
    (index_t i, index_t j) -> Vector3d {
        return {
            x_min + ((x_max - x_min) * i) / nx,
            y_min + ((y_max - y_min) * j) / ny,
            0.0
        };
    };

    auto neib_index = [x_period, y_period, nx, ny, get_index]
    (index_t i, index_t j, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j) : get_index((i - 1 + nx) % nx, j);
        }
        else if (side == Side2D::RIGHT) {
            return i == nx - 1 && !x_period ? get_index(i, j) : get_index((i + 1 + nx) % nx, j);
        }
        else if (side == Side2D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j) : get_index(i, (j - 1 + ny) % ny);
        }
        else if (side == Side2D::TOP) {
            return j == ny - 1 && !y_period ? get_index(i, j): get_index(i, (j + 1 + ny) % ny);
        }
        else {
            throw std::runtime_error("Strange side #142");
        }
    };

    // ------- Поехали заполнять -------
    m_dim = 2;
    m_adaptive = true;
    m_linear = true;
    m_axial = rect.axial();

    resize_amr(nx * ny);

    for (index_t ic = 0; ic < m_size; ++ic) {
        constexpr int n_nodes = 9;
        constexpr int n_faces = 8;

        auto[i, j] = get_index_pair(ic);

        next[ic] = ic;
        rank[ic] = 0;
        index[ic] = ic;

        flag[ic] = 0;
        level[ic] = 0;
        b_idx[ic] = ic;
        z_idx[ic] = 0;

        SqQuad quad(get_vertex(i, j),
                    get_vertex(i + 1, j),
                    get_vertex(i, j + 1),
                    get_vertex(i + 1, j + 1));

        center[ic] = quad.vs<0, 0>();
        volume[ic] = hx * hy;
        volume_alt[ic] = NAN;
        face_begin[ic] = ic * n_faces;
        node_begin[ic] = ic * n_nodes;
        face_begin[ic + 1] = (ic + 1) * n_faces;
        node_begin[ic + 1] = (ic + 1) * n_nodes;

        // INIT FACES
        for (auto iface: faces_range(ic)) {
            faces.set_undefined(iface);
            faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        faces.boundary[iface + Side2D::L] = i > 0 ? Boundary::INNER : left;
        faces.boundary[iface + Side2D::R] = i < nx - 1 ? Boundary::INNER : right;
        faces.boundary[iface + Side2D::B] = j > 0 ? Boundary::INNER : bottom;
        faces.boundary[iface + Side2D::T] = j < ny - 1 ? Boundary::INNER : top;

        for (auto side: Side2D::items()) {
            faces.adjacent.rank[iface + side] = 0;
            faces.adjacent.index[iface + side] = neib_index(i, j, side);
            faces.adjacent.alien[iface + side] = -1;
            faces.adjacent.basic[iface + side] = ic;
            faces.vertices[iface + side].fill(-1);
        }

        faces.normal[iface + Side2D::L] = -Vector3d::UnitX();
        faces.normal[iface + Side2D::R] =  Vector3d::UnitX();
        faces.normal[iface + Side2D::B] = -Vector3d::UnitY();
        faces.normal[iface + Side2D::T] =  Vector3d::UnitY();

        faces.center[iface + Side2D::L] = quad.vs<-1, 0>();
        faces.center[iface + Side2D::R] = quad.vs<+1, 0>();
        faces.center[iface + Side2D::B] = quad.vs<0, -1>();
        faces.center[iface + Side2D::T] = quad.vs<0, +1>();

        faces.area[iface + Side2D::L] = hy;
        faces.area[iface + Side2D::R] = hy;
        faces.area[iface + Side2D::B] = hx;
        faces.area[iface + Side2D::T] = hx;

        faces.vertices[iface + Side2D::L] = indexing::amr::sf(Side2D::L);
        faces.vertices[iface + Side2D::R] = indexing::amr::sf(Side2D::R);
        faces.vertices[iface + Side2D::B] = indexing::amr::sf(Side2D::B);
        faces.vertices[iface + Side2D::T] = indexing::amr::sf(Side2D::T);

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            verts[ic * n_nodes + jn] = quad[jn];
        }

        if (m_axial) {
            // "Альтернативный" объем ячейки и площади граней
            volume_alt[ic] = hx * hy * quad.vs<0, 0>().y();
            faces.area_alt[iface + Side2D::L] = hy * quad.vs<-1, 0>().y();
            faces.area_alt[iface + Side2D::R] = hy * quad.vs<+1, 0>().y();
            faces.area_alt[iface + Side2D::B] = hx * quad.vs< 0,-1>().y();
            faces.area_alt[iface + Side2D::T] = hx * quad.vs< 0,+1>().y();

            // Смещения барицентров, есть необходимость?
            center[ic].y() += hy*hy / (12.0 * quad.vs<0, 0>().y());
            faces.center[iface + Side2D::L].y() += hy*hy / (12.0 * quad.vs<-1, 0>().y());
            faces.center[iface + Side2D::R].y() += hy*hy / (12.0 * quad.vs<+1, 0>().y());
        }
    }
}

AmrCells::AmrCells(const generator::Cuboid& c) {
    bool x_period = c.periodic_along_x();
    bool y_period = c.periodic_along_y();
    bool z_period = c.periodic_along_z();

    double hx = (c.x_max() - c.x_min()) / c.nx();
    double hy = (c.y_max() - c.y_min()) / c.ny();
    double hz = (c.z_max() - c.z_min()) / c.nz();

    const auto [left, right, bottom, top, back, front] = c.bounds();

    auto get_index = [c]
    (index_t i, index_t j, index_t k) -> index_t {
        return c.nz() * (c.ny() * i + j) + k;
    };

    auto get_index_pair = [c]
    (index_t n) -> std::array<index_t, 3> {
        return {(n / c.nz()) / c.ny(), (n / c.nz()) % c.ny(), n % c.nz()};
    };

    auto get_vertex = [c]
    (index_t i, index_t j, index_t k) -> Vector3d {
        return {
                c.x_min() + ((c.x_max() - c.x_min()) * i) / c.nx(),
                c.y_min() + ((c.y_max() - c.y_min()) * j) / c.ny(),
                c.z_min() + ((c.z_max() - c.z_min()) * k) / c.nz()
        };
    };

    auto neib_index = [c, x_period, y_period, z_period, get_index]
    (index_t i, index_t j, index_t k, Side3D side) -> index_t {
        if (side == Side3D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j, k) : get_index((i - 1 + c.nx()) % c.nx(), j, k);
        }
        else if (side == Side3D::RIGHT) {
            return i == c.nx() - 1 && !x_period ? get_index(i, j, k) : get_index((i + 1) % c.nx(), j, k);
        }
        else if (side == Side3D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j, k) : get_index(i, (j - 1 + c.ny()) % c.ny(), k);
        }
        else if (side == Side3D::TOP) {
            return j == c.ny() - 1 && !y_period ? get_index(i, j, k): get_index(i, (j + 1) % c.ny(), k);
        }
        else if (side == Side3D::BACK) {
            return k == 0 && !z_period ? get_index(i, j, k) : get_index(i, j, (k - 1 + c.nz()) % c.nz());
        }
        else if (side == Side3D::FRONT) {
            return k == c.nz() - 1 && !z_period ? get_index(i, j, k): get_index(i, j, (k + 1) % c.nz());
        }
        else {
            throw std::runtime_error("Strange side #265");
        }
    };

    set_dimension(3);
    set_adaptive(true);
    set_linear(true);
    set_axial(false);

    resize_amr(c.nx() * c.ny() * c.nz());

    for (index_t ic = 0; ic < m_size; ++ic) {
        constexpr int n_faces = 24;
        constexpr int n_nodes = 27;

        auto[i, j, k] = get_index_pair(ic);

        next[ic] = ic;
        rank[ic] = 0;
        index[ic] = ic;

        flag[ic] = 0;
        level[ic] = 0;
        b_idx[ic] = ic;
        z_idx[ic] = 0;

        SqCube cube(get_vertex(i,   j,   k),
                    get_vertex(i+1, j,   k),
                    get_vertex(i,   j+1, k),
                    get_vertex(i+1, j+1, k),
                    get_vertex(i,   j,   k+1),
                    get_vertex(i+1, j,   k+1),
                    get_vertex(i,   j+1, k+1),
                    get_vertex(i+1, j+1, k+1));

        center[ic] = cube.vs<0, 0, 0>();
        volume[ic] = hx * hy * hz;
        volume_alt[ic] = NAN;
        face_begin[ic] = ic * n_faces;
        node_begin[ic] = ic * n_nodes;
        face_begin[ic + 1] = (ic + 1) * n_faces;
        node_begin[ic + 1] = (ic + 1) * n_nodes;

        // INIT FACES
        for (auto iface: faces_range(ic)) {
            faces.set_undefined(iface);
            faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        faces.boundary[iface + Side3D::L] = i > 0 ? Boundary::INNER : left;
        faces.boundary[iface + Side3D::R] = i < c.nx() - 1 ? Boundary::INNER : right;
        faces.boundary[iface + Side3D::B] = j > 0 ? Boundary::INNER : bottom;
        faces.boundary[iface + Side3D::T] = j < c.ny() - 1 ? Boundary::INNER : top;
        faces.boundary[iface + Side3D::Z] = k > 0 ? Boundary::INNER : back;
        faces.boundary[iface + Side3D::F] = k < c.nz() - 1 ? Boundary::INNER : front;

        for (auto side: Side3D::items()) {
            faces.adjacent.rank[iface + side] = 0;
            faces.adjacent.index[iface + side] = neib_index(i, j, k, side);
            faces.adjacent.alien[iface + side] = -1;
            faces.adjacent.basic[iface + side] = ic;
            faces.vertices[iface + side].fill(-1);
        }

        faces.normal[iface + Side3D::L] = -Vector3d::UnitX();
        faces.normal[iface + Side3D::R] =  Vector3d::UnitX();
        faces.normal[iface + Side3D::B] = -Vector3d::UnitY();
        faces.normal[iface + Side3D::T] =  Vector3d::UnitY();
        faces.normal[iface + Side3D::Z] = -Vector3d::UnitZ();
        faces.normal[iface + Side3D::F] =  Vector3d::UnitZ();

        faces.center[iface + Side3D::L] = cube.vs<-1, 0, 0>();
        faces.center[iface + Side3D::R] = cube.vs<+1, 0, 0>();
        faces.center[iface + Side3D::B] = cube.vs< 0,-1, 0>();
        faces.center[iface + Side3D::T] = cube.vs< 0,+1, 0>();
        faces.center[iface + Side3D::Z] = cube.vs< 0, 0,-1>();
        faces.center[iface + Side3D::F] = cube.vs< 0, 0,+1>();

        faces.area[iface + Side3D::L] = hy * hz;
        faces.area[iface + Side3D::R] = hy * hz;
        faces.area[iface + Side3D::B] = hx * hz;
        faces.area[iface + Side3D::T] = hx * hz;
        faces.area[iface + Side3D::Z] = hx * hy;
        faces.area[iface + Side3D::F] = hx * hy;

        faces.vertices[iface + Side3D::L] = indexing::amr::sf(Side3D::L);
        faces.vertices[iface + Side3D::R] = indexing::amr::sf(Side3D::R);
        faces.vertices[iface + Side3D::B] = indexing::amr::sf(Side3D::B);
        faces.vertices[iface + Side3D::T] = indexing::amr::sf(Side3D::T);
        faces.vertices[iface + Side3D::Z] = indexing::amr::sf(Side3D::Z);
        faces.vertices[iface + Side3D::F] = indexing::amr::sf(Side3D::F);

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            verts[ic * n_nodes + jn] = cube[jn];
        }
    }
}

AmrCells::AmrCells(const Grid& grid) {
    m_dim = grid.dimension();
    m_adaptive = grid.adaptive();
    m_linear = true;
    m_axial = false;

    if (!grid.has_faces()) {
        throw std::runtime_error("Grid was built with wrong options, need faces per cell");
    }

    if (!m_adaptive) {
        resize(grid.n_cells(), grid.total_faces_per_cell(), grid.total_nodes_per_cell());
    }
    else {
        resize_amr(grid.n_cells());
    }

    const auto& nodes = grid.nodes();
    const auto& cells = grid.cells();

    face_begin[0] = 0;
    node_begin[0] = 0;
    for (index_t ic = 0; ic < m_size; ++ic) {
        rank[ic] = 0;
        index[ic] = ic;

        flag[ic] = 0;
        b_idx[ic] = ic;
        z_idx[ic] = 0;
        level[ic] = 0;

        volume[ic] = cells[ic].volume();
        center[ic] = cells[ic].centroid();

        if (m_axial) {
            // m_locals.volume_alt[ic] = geom.cell_volumes_alt[ic];
        }

        // Число узлов и граней ячейки
        int n_nodes = cells[ic].n_nodes();
        int n_faces = cells[ic].n_faces();
        int max_faces = n_faces;
        if (m_adaptive) {
            max_faces = m_dim < 3 ? Side2D::n_subfaces() : Side3D::n_subfaces();
        }

        // Выставить индексы грани
        face_begin[ic + 1] = face_begin[ic] + max_faces;
        if (cells[ic].type() != CellType::POLYHEDRON) {
            faces.insert(face_begin[ic], cells[ic].type(), max_faces);
        }
        else {
            auto iface = face_begin[ic];
            for (int i = 0; i < n_faces; ++i) {
                const auto& face = cells[ic].get_face(i);
                faces.vertices[iface + i].fill(-1);
                for (int j = 0; j < face.n_nodes(); ++j) {
                    faces.vertices[iface + i][j] = face.node_idx(j);
                }
                faces.set_undefined(iface + i);
            }
        }
        const auto& node_ids = cells[ic].nodes();

        node_begin[ic + 1] = node_begin[ic] + n_nodes;
        for (int i = 0; i < n_nodes; ++i) {
            verts[node_begin[ic] + i] = nodes[node_ids[i]].pos;
        }

        // Геометрия
        for (int i = 0; i < n_faces; ++i) {
            const auto& face = cells[ic].get_face(i);

            int iface = face_begin[ic] + i;
            faces.boundary[iface] = face.bc();
            faces.area[iface]     = face.area();
            faces.center[iface]   = face.center();
            faces.normal[iface]   = face.normal();

            if (m_axial) {
                //faces.area_alt[iface] = grid_geom.face_areas_alt[jface];
            }
        }

        // Смежность
        for (int i = 0; i < n_faces; ++i) {
            const auto& face = cells[ic].get_face(i);

            auto iface = face_begin[ic] + i;

            faces.adjacent.rank[iface]  = 0;
            faces.adjacent.index[iface] = face.neib();
            faces.adjacent.alien[iface] = -1;
            faces.adjacent.basic[iface] = ic;

            if (faces.is_boundary(iface)) {
                faces.adjacent.index[iface] = ic;
            }
        }
    }
}

} // namespace zephyr::mesh