#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/mesh/euler/soa_mesh.h>

namespace zephyr::mesh {



CellIt AmrCells::begin() { return {this, 0}; }

CellIt AmrCells::end() { return {this, size()}; }

QCell AmrCells::operator[](index_t cell_idx) {
    return {this, cell_idx};
}


void AmrCells::initialize(const Strip& gen) {
    bool x_period = gen.periodic_along_x();

    double xmin = gen.x_min();
    double ymin = gen.y_min();
    double xmax = gen.x_max();
    double ymax = gen.y_max();

    index_t nx = gen.nx();
    index_t ny = 1;

    double hx = (xmax - xmin) / nx;
    double hy = (ymax - ymin) / ny;

    auto get_vertex = [=](index_t i, index_t j) -> Vector3d {
        return {
                xmin + ((xmax - xmin) * i) / nx,
                ymin + ((ymax - ymin) * j) / ny,
                0.0
        };
    };

    auto neib_index = [=](index_t i, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  i : (i - 1 + nx) % nx;
        }
        else if (side == Side2D::RIGHT) {
            return i == nx - 1 && !x_period ? i : (i + 1 + nx) % nx;
        }
        else {
            throw std::runtime_error("Strange side #153");
        }
    };

    m_size = nx * ny;

    dim = 2;
    adaptive = true;

    linear = true;
    axial = false;

    int n_faces = 8;
    int n_nodes = 9;

    resize(m_size);

    for (index_t ic = 0; ic < m_size; ++ic) {
        next[ic] = ic;
        rank[ic] = 0;
        index[ic] = ic;

        flag[ic] = 0;
        level[ic] = 0;
        b_idx[ic] = ic;
        z_idx[ic] = 0;

        SqQuad quad(
                get_vertex(ic, 0),
                get_vertex(ic + 1, 0),
                get_vertex(ic, 1),
                get_vertex(ic + 1, 1));

        center[ic] = quad.vs<0, 0>();
        volume[ic] = hx * hy;
        volume_alt[ic] = NAN;
        face_begin[ic] = ic * n_faces;
        node_begin[ic] = ic * n_nodes;
        face_begin[ic + 1] = face_begin[ic] + n_faces;
        node_begin[ic + 1] = face_begin[ic] + n_nodes;

        // INIT FACES
        for (index_t iface: faces_range(ic)) {
            faces.set_undefined(iface);
            faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        faces.boundary[iface + Side2D::L] = ic > 0 ? Boundary::ORDINARY : gen.bounds().left;
        faces.boundary[iface + Side2D::R] = ic < nx - 1 ? Boundary::ORDINARY : gen.bounds().right;

        for (auto side: {Side2D::LEFT, Side2D::RIGHT}) {
            faces.adjacent.rank[iface + side] = 0;
            faces.adjacent.index[iface + side] = neib_index(ic, side);
            faces.adjacent.alien[iface + side] = -1;
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

        faces.vertices[iface + Side2D::L] = Side2D::L.sf();
        faces.vertices[iface + Side2D::R] = Side2D::R.sf();
        faces.vertices[iface + Side2D::B] = Side2D::B.sf();
        faces.vertices[iface + Side2D::T] = Side2D::T.sf();

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            verts[ic * n_nodes + jn] = quad[jn];
        }
    }
}

void AmrCells::initialize(const Rectangle& gen) {
    bool x_period = gen.periodic_along_x();
    bool y_period = gen.periodic_along_y();

    double xmin = gen.x_min();
    double ymin = gen.y_min();
    double xmax = gen.x_max();
    double ymax = gen.y_max();

    index_t nx = gen.nx();
    index_t ny = gen.ny();

    double hx = (xmax - xmin) / nx;
    double hy = (ymax - ymin) / ny;

    auto get_index = [=](index_t i, index_t j) -> index_t {
        return i * ny + j;
    };

    auto get_index_pair = [=](index_t n) -> std::array<index_t, 2> {
        return {n / ny, n % ny};
    };

    auto get_vertex = [=](index_t i, index_t j) -> Vector3d {
        return {
                xmin + ((xmax - xmin) * i) / nx,
                ymin + ((ymax - ymin) * j) / ny,
                0.0
        };
    };

    auto neib_index = [=](index_t i, index_t j, Side2D side) -> index_t {
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

    m_size = nx * ny;

    dim = 2;
    adaptive = true;
    linear = true;
    axial = false;

    resize(m_size);

    int n_faces = 8;
    int n_nodes = 9;

    for (index_t ic = 0; ic < m_size; ++ic) {
        auto[i, j] = get_index_pair(ic);

        next[ic] = ic;
        rank[ic] = 0;
        index[ic] = ic;

        flag[ic] = 0;
        level[ic] = 0;
        b_idx[ic] = ic;
        z_idx[ic] = 0;

        SqQuad quad(
                get_vertex(i, j),
                get_vertex(i + 1, j),
                get_vertex(i, j + 1),
                get_vertex(i + 1, j + 1));

        center[ic] = quad.vs<0, 0>();
        volume[ic] = hx * hy;
        volume_alt[ic] = NAN;
        face_begin[ic] = ic * n_faces;
        node_begin[ic] = ic * n_nodes;
        face_begin[ic + 1] = face_begin[ic] + n_faces;
        node_begin[ic + 1] = face_begin[ic] + n_nodes;

        // INIT FACES
        for (index_t iface: faces_range(ic)) {
            faces.set_undefined(iface);
            faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        faces.boundary[iface + Side2D::L] = i > 0 ? Boundary::ORDINARY : gen.bounds().left;
        faces.boundary[iface + Side2D::R] = i < nx - 1 ? Boundary::ORDINARY : gen.bounds().right;
        faces.boundary[iface + Side2D::B] = j > 0 ? Boundary::ORDINARY : gen.bounds().bottom;
        faces.boundary[iface + Side2D::T] = j < ny - 1 ? Boundary::ORDINARY : gen.bounds().top;

        for (auto side: Side2D::items()) {
            faces.adjacent.rank[iface + side] = 0;
            faces.adjacent.index[iface + side] = neib_index(i, j, side);
            faces.adjacent.alien[iface + side] = -1;
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

        faces.vertices[iface + Side2D::L] = Side2D::L.sf();
        faces.vertices[iface + Side2D::R] = Side2D::R.sf();
        faces.vertices[iface + Side2D::B] = Side2D::B.sf();
        faces.vertices[iface + Side2D::T] = Side2D::T.sf();

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            verts[ic * n_nodes + jn] = quad[jn];
        }
    }
}

void AmrCells::initialize(const Cuboid& gen) {
    throw std::runtime_error("AMRCells Initialize cuboid");
}

} // namespace zephyr::mesh