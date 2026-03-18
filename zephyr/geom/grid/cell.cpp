#include <limits>
#include <stdexcept>
#include <utility>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/indexing.h>

#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>
#include <zephyr/geom/primitives/polyhedron.h>

namespace zephyr::geom {

// ------------------------------------------------------- NODE -------------------------------------------------------

Node::Node(const Vector3d& v) : pos(v) { }

Node::Ptr Node::create(const Vector3d& v) {
    return std::make_shared<Node>(v);
}

Node::Ptr Node::create(double x, double y) {
    return std::make_shared<Node>(Vector3d{x, y, 0.0});
}

bool Node::is_finite() const {
    return std::isfinite(pos[0]) && std::isfinite(pos[1]) && std::isfinite(pos[2]);
}

bool Node::operator==(const Node& other) const {
    return m_id == other.m_id && pos == other.pos && bc == other.bc;
}

bool Node::operator!=(const Node& other) const {
    return !(*this == other);
}

// ------------------------------------------------------- FACE -------------------------------------------------------

void Face::set_nodes(std::span<const int> node_ids) {
    m_nodes.resize(node_ids.size());
    for (int i = 0; i < node_ids.size(); ++i) {
        m_nodes[i] = static_cast<node_id_t>(node_ids[i]);
    }
}

void Face::set_bc(Boundary bc) {
    m_bc = m_neib != invalid_id ? Boundary::INNER : bc;
}

void Face::set_neib(id_t neib_id, int face_id) {
    if (neib_id != invalid_id) {
        m_neib = neib_id;
        m_twin = face_id;
        m_bc = Boundary::INNER;
    }
}

void Face::calc_geom_line(const std::vector<Vector3d>& vertices, const Vector3d& view) {
    if (m_nodes.size() != 2) {
        throw std::runtime_error("Face::calc_geom_line: invalid number of nodes");
    }
    Line line{vertices[m_nodes[0]], vertices[m_nodes[1]]};

    m_center = line.center();
    m_area_n = line.area_n(view);
}

void Face::calc_geom_poly(const std::vector<Vector3d>& vertices, const Vector3d& view) {

}

void Face::calc_geom_amr2d(const std::vector<Vector3d>& vertices, const Vector3d& view) {
    if (m_nodes.size() != 3) {
        throw std::runtime_error("Face::calc_geom_amr2d: invalid number of nodes");
    }
    SqLine line{vertices[m_nodes[0]], vertices[m_nodes[1]], vertices[m_nodes[2]]};
    m_center = line.center();
    m_area_n = line.length() * line.normal(view);
}

void Face::calc_geom_amr3d(const std::vector<Vector3d>& vertices, const Vector3d& view) {

}

// ------------------------------------------------------- CELL -------------------------------------------------------

Cell::Cell(CellType type, std::vector<id_t>&& node_ids)
    : m_type(type), m_nodes(std::move(node_ids)) {
    int default_n_nodes = indexing::n_nodes(m_type);
    if (default_n_nodes >= 0 && default_n_nodes != m_nodes.size()) {
        throw std::runtime_error("Cell::Cell: bad nodes count.");
    }
}

Cell::Cell(CellType type, std::initializer_list<id_t> node_ids)
    : m_type(type) {
    std::ranges::copy(node_ids, std::back_inserter(m_nodes));
    int default_n_nodes = indexing::n_nodes(m_type);
    if (default_n_nodes >= 0 && default_n_nodes != m_nodes.size()) {
        throw std::runtime_error("Cell::Cell: bad nodes count.");
    }
}

void Cell::set_face_bc(const std::vector<Boundary>& face_bc) {
    if (face_bc.empty()) return;

    if (m_faces.empty()) {
        if (m_type != CellType::POLYHEDRON) {
            init_faces();
        }
        else {
            m_faces.resize(face_bc.size());
        }
    }
    if (m_faces.size() != face_bc.size()) {
        throw std::runtime_error("Cell::set_face_bc: bad boundary conditions count.");
    }
    for (int i = 0; i < face_bc.size(); ++i) {
        m_faces[i].set_bc(face_bc[i]);
    }
}

void Cell::init_faces() {
    int n_faces = indexing::n_faces(m_type);
    if (n_faces < 0) {
        // dynamic face type
        if (m_type == CellType::POLYGON) {
            n_faces = n_nodes();
        }
        else {
            throw std::runtime_error("Can't init faces for type POLYHEDRON");
        }
    }

    m_faces.resize(n_faces);

    using namespace indexing;
    switch (m_type) {
        // 2D
        case CellType::TRIANGLE: {
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(tri::face_nodes[i], 2));
            }
        } break;

        case CellType::QUAD: {
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(quad::face_nodes[i], 2));
            }
        } break;

        case CellType::POLYGON:{
            if (n_faces != m_nodes.size()) {
                throw std::runtime_error("Cell::Cell: bad faces count (polygon.n_nodes != polygon.n_faces).");
            }
            for (int i = 0; i < n_faces; ++i) {
                std::array node_ids{i, (i + 1) % n_faces};
                m_faces[i].set_nodes(node_ids);
            }
        } break;

        case CellType::AMR2D: {
            for (int i = 0; i < n_faces; ++i) {
                // криволинейные грани на 3 узлах
                m_faces[i].set_nodes(std::span(amr::sqface_nodes2d[i], 3));
            }
        } break;

        // 3D classic
        case CellType::TETRA: {
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(tetra::face_nodes[i], 3));
            }
        } break;

        case CellType::PYRAMID: {
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(pyramid::face_nodes[i], pyramid::face_size[i]));
            }
        } break;

        case CellType::WEDGE: {
            // 3 quads + 2 triangles
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(wedge::face_nodes[i], wedge::face_size[i]));
            }
        } break;

        case CellType::HEXAHEDRON: {
            for (int i = 0; i < n_faces; ++i) {
                m_faces[i].set_nodes(std::span(hex::face_nodes[i], 4));
            }
        } break;

        case CellType::AMR3D: {
            for (int i = 0; i < n_faces; ++i) {
                // криволинейные грани на 9 узлах
                m_faces[i].set_nodes(std::span(amr::sqface_nodes3d[i], 9));
            }
        } break;

        default: throw std::runtime_error("Cell::init_faces: unsupported cell type.");
    }
}

void Cell::set_faces(const std::vector<std::vector<int>>& face_ids) {
    if (m_type != CellType::POLYHEDRON) {
        throw std::runtime_error("Cell::set_faces: set_faces only for polyhedron type.");
    }
    m_faces.resize(face_ids.size());
    for (int i = 0; i < face_ids.size(); ++i) {
        m_faces[i].set_nodes(face_ids[i]);
    }
}

void Cell::set_neib(int iface, id_t neib_id, int face_id) {
    m_faces[iface].set_neib(neib_id, face_id);
}

void Cell::replace_nodes(std::vector<id_t>&& new_nodes) {
    if (new_nodes.size() != m_nodes.size()) {
        throw std::runtime_error("Cell:replace_nodes must have the same number of nodes");
    }
    m_nodes = std::move(new_nodes);
}

void Cell::mirror() {
    // Необходимо развернуть ячейки
    if (indexing::get_dimension(m_type) == 2) {
        if (m_type == CellType::TRIANGLE || m_type == CellType::POLYGON) {
            // Узлы в обратном порядке, первый узел на месте
            std::reverse(m_nodes.begin() + 1, m_nodes.end());

            // Индексы на гранях остаются те же самые,
            // но граничные условия надо развернуть
            for (int i = 0; i < m_faces.size() / 2; ++i) {
                int j = m_faces.size() - i - 1;
                auto bc_i = m_faces[i].bc();
                auto bc_j = m_faces[j].bc();
                m_faces[i].set_bc(bc_j);
                m_faces[j].set_bc(bc_i);
            }
        }
        else if (m_type == CellType::QUAD) {
            std::swap(m_nodes[0], m_nodes[1]);
            std::swap(m_nodes[2], m_nodes[3]);
            auto bc_L = m_faces[Side2D::L].bc();
            auto bc_R = m_faces[Side2D::R].bc();
            m_faces[Side2D::L].set_bc(bc_R);
            m_faces[Side2D::R].set_bc(bc_L);
        }
        else if (m_type == CellType::AMR2D) {
            throw std::runtime_error("Grid::mirror: no amr implementation");
        }
        else {
            throw std::runtime_error("Grid::mirror: bad cell type.");
        }
    }
    else {
        throw std::runtime_error("Cell::mirror: 3D mirror not implemented");
    }
}

Vector3d Cell::face_center(const std::vector<Node>& grid_nodes, int iface) const {
    if (m_faces.empty()) { return nanvec(); }
    if (!m_faces[iface].has_nodes()) { return nanvec(); }

    Vector3d fc = Vector3d::Zero();
    for (auto nid: m_faces[iface].nodes()) {
        fc += grid_nodes[m_nodes[nid]].pos;
    }
    fc /= m_faces[iface].n_nodes();
    return fc;
}

Vector3d Cell::center(const std::vector<Node>& grid_nodes) const {
    Vector3d c = Vector3d::Zero();
    for (id_t j: m_nodes) {
        c += grid_nodes[j].pos;
    }
    c /= m_nodes.size();
    return c;
}

void Cell::calc_geom(const std::vector<Node>& grid_nodes) {
    if (indexing::get_dimension(m_type) == 2) {
        if (m_type != CellType::AMR2D) {
            std::vector<Vector3d> vertices;
            vertices.reserve(m_nodes.size());
            for (const auto& nid: m_nodes) {
                vertices.push_back(grid_nodes[nid].pos);
            }
            Polygon poly(std::move(vertices));
            m_volume = poly.area();
            m_center = poly.centroid(m_volume);
            for (auto& face: m_faces) {
                face.calc_geom_line(poly.vertices(), m_center);
            }
        }
        else {
            std::vector<Vector3d> vertices;
            vertices.reserve(m_nodes.size());
            for (const auto& nid: m_nodes) {
                vertices.push_back(grid_nodes[nid].pos);
            }
            SqQuad quad(std::span<const Vector3d, 9>{vertices});
            m_volume = quad.area();
            m_center = quad.centroid(m_volume);
            for (auto& face: m_faces) {
                face.calc_geom_amr2d(vertices, m_center);
            }
        }
    }
    else {
        if (m_type != CellType::AMR3D) {
            std::vector<Vector3d> vertices;
            vertices.reserve(m_nodes.size());
            for (const auto& nid: m_nodes) {
                vertices.push_back(grid_nodes[nid].pos);
            }
            std::vector<std::vector<int>> face_ids(m_faces.size());
            for (int i = 0; i < m_faces.size(); ++i) {
                const auto& face = m_faces[i];
                face_ids[i].resize(face.n_nodes());
                for (int k = 0; k < face.n_nodes(); ++k) {
                    face_ids[i][k] = face.node_idx(k);
                }
            }
            Polyhedron poly(vertices, face_ids);
            m_volume = poly.volume();
            m_center = poly.centroid(m_volume);
            for (int i = 0; i < m_faces.size(); ++i) {
                m_faces[i].set_center(poly.face_center(i));
                m_faces[i].set_area_n(poly.face_area_n(i));
            }
        }
        else {
            throw std::runtime_error("not implemented #1654");
        }
    }
}

} // namespace zephyr::geom