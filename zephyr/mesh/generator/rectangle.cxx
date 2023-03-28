#include <iostream>
#include <algorithm>

#include <zephyr/geom/cell.h>

#include <zephyr/geom/box.h>
#include <zephyr/mesh/generator/rectangle.h>


namespace zephyr { namespace mesh { namespace generator {

#ifdef ZEPHYR_ENABLE_MPI
using zephyr::network::mpi::Network;
#endif

#ifdef ZEPHYR_ENABLE_YAML
Rectangle::Rectangle(YAML::Node config)
    : Base("rectangle"),
      m_xmin(0.0), m_xmax(1.0), m_ymin(0.0), m_ymax(1.0), m_nx(0), m_ny(0),
      m_left_flag(FaceFlag::WALL), m_right_flag(FaceFlag::WALL),
      m_bottom_flag(FaceFlag::WALL), m_top_flag(FaceFlag::WALL),
      m_voronoi(false) {

    if (!config["geometry"]) {
        throw std::runtime_error("Mesh config doesn't contain 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();

    if (!config["boundary"]) {
        throw std::runtime_error("Mesh config doesn't contain 'boundary'");
    }
    m_left_flag = boundary_from_string(config["boundary"]["left"].as<std::string>());
    m_right_flag = boundary_from_string(config["boundary"]["right"].as<std::string>());
    m_bottom_flag = boundary_from_string(config["boundary"]["bottom"].as<std::string>());
    m_top_flag = boundary_from_string(config["boundary"]["top"].as<std::string>());

    if (config["voronoi"]) {
        m_voronoi = config["voronoi"].as<bool>();
    }

    if (config["cells"]) {
        set_size(config["cells"].as<int>());
    }
    else {
        int nx(0), ny(0);
        if (config["cells_per_x"]) {
            nx = config["cells_per_x"].as<int>();
        }
        if (config["cells_per_y"]) {
            ny = config["cells_per_y"].as<int>();
        }
        if (nx > 0 && ny > 0) {
            set_sizes(nx, ny);
        } else if (nx > 0) {
            set_nx(nx);
        } else if (ny > 0) {
            set_ny(ny);
        } else {
            throw std::runtime_error("Rectangle error: Strange mesh sizes");
        }
    }
}
#endif

Rectangle::Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi) :
        Base("rectangle"),
        m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax), m_nx(0), m_ny(0),
        m_left_flag(FaceFlag::UNDEFINED), m_right_flag(FaceFlag::UNDEFINED),
        m_bottom_flag(FaceFlag::UNDEFINED), m_top_flag(FaceFlag::UNDEFINED),
        m_voronoi(voronoi) {
    check_params();
}

Box Rectangle::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void Rectangle::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Rectangle: Nx < 1");
    }
    m_nx = nx;
    m_ny = int(round(nx * (m_ymax - m_ymin) / (m_xmax - m_xmin)));
    compute_size();
}

void Rectangle::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Rectangle: Ny < 1");
    }
    m_ny = ny;
    m_nx = int(round(ny * (m_xmax - m_xmin) / (m_ymax - m_ymin)));
    compute_size();
}

void Rectangle::set_sizes(int nx, int ny) {
    if (nx < 1 || ny < 1) {
        throw std::runtime_error("Rectangle: Nx < 1 or Ny < 1");
    }
    m_nx = nx;
    m_ny = ny;
    compute_size();

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    double dmax = std::max(dx, dy);
    double dmin = std::min(dx, dy);
    if (dmax / dmin > 2.0) {
        std::cerr << "Rectangle Warning: Large aspect ratio (> 2)\n";
    }
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Rectangle Warning: Huge aspect ratio (> 1000)\n";
    }
}

void Rectangle::set_size(int N) {
    double d = std::sqrt((m_xmax - m_xmin) * (m_ymax - m_ymin) / N);
    m_nx = int(round((m_xmax - m_xmin) / d));
    m_ny = int(round((m_ymax - m_ymin) / d));
    compute_size();
}

void Rectangle::set_boundary_flags(FaceFlag left, FaceFlag right, FaceFlag bottom, FaceFlag top) {
    m_left_flag = left;
    m_right_flag = right;
    m_bottom_flag = bottom;
    m_top_flag = top;
    if (periodic_along_x()) {
        m_left_flag = m_right_flag = FaceFlag::PERIODIC;
    }
    if (periodic_along_y()) {
        m_bottom_flag = m_top_flag = FaceFlag::PERIODIC;
    }
}

double Rectangle::x_min() const {
    return m_xmin;
}

double Rectangle::x_max() const {
    return m_xmax;
}

double Rectangle::y_min() const {
    return m_ymin;
}

double Rectangle::y_max() const {
    return m_ymax;
}

int Rectangle::nx() const {
    return m_nx;
}

int Rectangle::ny() const {
    return m_ny;
}

bool Rectangle::periodic_along_x() const {
    return m_left_flag == FaceFlag::PERIODIC || m_right_flag == FaceFlag::PERIODIC;
}

bool Rectangle::periodic_along_y() const {
    return m_bottom_flag == FaceFlag::PERIODIC || m_top_flag == FaceFlag::PERIODIC;
}

void Rectangle::check_params() const {
    if (m_xmin >= m_xmax) {
        std::cerr << "Rectangle Error: x_min >= x_max\n";
        throw std::runtime_error("Rectangle Error: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        std::cerr << "Rectangle Error: y_min >= y_max\n";
        throw std::runtime_error("Rectangle Error: y_min >= y_max");
    }
}

void Rectangle::compute_size() {
    if (m_voronoi) {
        int Nb = m_nx * m_ny;
        double DX = m_xmax - m_xmin;
        double DY = m_ymax - m_ymin;
        int Ny = int(std::floor(std::sqrt(std::sqrt(3.0) * DY * Nb / 2.0 / DX))) + 1;
        double h = DY / 2.0 / Ny;
        double D = h / std::sqrt(0.75);

        int Nx = int(std::floor(DX / 1.5 / D));

        m_size = (Nx + 1) * Ny + (Nx + 1) / 2;
    } else {
        m_size = m_nx * m_ny;
    }
}

void Rectangle::initialize(Storage &cells) {
    cells.resize(m_size);

    if (m_voronoi) {
        init_voronoi(cells, Part(m_size));
    } else {
        init_classic(cells, Part(m_size));
    }
}

void Rectangle::initialize(Storage &cells, Part part) {
    if (m_voronoi) {
        init_voronoi(cells, part);
    } else {
        init_classic(cells, part);
    }
}

void Rectangle::init_classic(Storage &cells, Part part) const {
    using zephyr::geom::Side;
    using zephyr::geom::Cell;
    using zephyr::geom::ShortList2D;

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    auto get_n = [this](int i, int j) -> int {
        return ((i + m_nx) % m_nx) * m_ny + (j + m_ny) % m_ny;
    };

    int counter = 0;
    for (int n = part.from; n < part.to; ++n, ++counter) {
        int i = n / m_ny;
        int j = n % m_ny;

        // Границы ячейки
        double x1 = m_xmin + i * dx;
        double x2 = m_xmin + (i + 1.0) * dx;

        double y1 = m_ymin + j * dy;
        double y2 = m_ymin + (j + 1.0) * dy;

        ShortList2D verts = {
                Vector3d(x1, y1, 0.0), Vector3d(x2, y1, 0.0),
                Vector3d(x1, y2, 0.0), Vector3d(x2, y2, 0.0)
        };

        Cell cell(verts);

        // Данные AMR
        cell.amr.base_id = n;
        cell.amr.level   = 0;
        cell.amr.flag    = 0;
        cell.amr.z       = 0;

        // Флаги граничных условий
        auto ordinary = FaceFlag::ORDINARY;

        cell.faces[Side::L].boundary = i > 0 ? ordinary : m_left_flag;
        cell.faces[Side::R].boundary = i < m_nx - 1 ? ordinary : m_right_flag;
        cell.faces[Side::B].boundary = j > 0 ? ordinary : m_bottom_flag;
        cell.faces[Side::T].boundary = j < m_ny - 1 ? ordinary : m_top_flag;

        cell.faces[Side::L].adjacent.rank = 0;
        cell.faces[Side::R].adjacent.rank = 0;
        cell.faces[Side::B].adjacent.rank = 0;
        cell.faces[Side::T].adjacent.rank = 0;

        cell.faces[Side::L].adjacent.ghost = 0;
        cell.faces[Side::R].adjacent.ghost = 0;
        cell.faces[Side::B].adjacent.ghost = 0;
        cell.faces[Side::T].adjacent.ghost = 0;

        cell.faces[Side::L].adjacent.index = get_n(i - 1, j);
        cell.faces[Side::R].adjacent.index = get_n(i + 1, j);
        cell.faces[Side::B].adjacent.index = get_n(i, j - 1);
        cell.faces[Side::T].adjacent.index = get_n(i, j + 1);

        cells[counter].geom() = cell;
    }
}

void Rectangle::init_voronoi(Storage &cells, Part part) const {
    int Nb = m_nx * m_ny;
    double DX = m_xmax - m_xmin;
    double DY = m_ymax - m_ymin;
    int Ny = int(std::floor(std::sqrt(std::sqrt(3.0) * DY * Nb / 2.0 / DX))) + 1;
    double h = DY / 2.0 / Ny;
    double D = h / std::sqrt(0.75);

    int Nx = int(std::floor(DX / 1.5 / D));

    double x_shift = (DX - Nx * 1.5 * D) / 2.0;

    int n = 0;
    int counter = 0;

    /*
    std::vector<Vector3d> centers(cells.size());

    for (int i = 0; i < (Nx + 1) / 2; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if (part.from <= n && n < part.to) {
                centers[counter].x = 3 * D * i + x_shift;
                centers[counter].y = (2.0 * j + 1.0) * h;
                centers[counter].z = 0.0;
                ++counter;
            }

            ++n;
            if (n > m_size) {
                break;
            }
        }
        if (n > m_size) {
            break;
        }

        for (int j = 0; j <= Ny; ++j) {
            if (part.from <= n && n < part.to) {
                centers[counter].x = 3 * D * (i + 0.5) + x_shift;
                centers[counter].y = 2.0 * j * h;
                centers[counter].z = 0.0;
                ++counter;
            }

            ++n;
            if (n > m_size) {
                break;
            }
        }
        if (n > m_size) {
            break;
        }
    }
    for (int j = 0; j < Ny; ++j) {
        if (part.from <= n && n < part.to) {
            centers[counter].x = 3 * D * ((Nx + 1) / 2) + x_shift;
            centers[counter].y = (2.0 * j + 1.0) * h;
            centers[counter].z = 0.0;
            ++counter;
        }

        ++n;
        if (n > m_size) {
            break;
        }
    }

    double d = 0.5 * D;

    using zephyr::math::geom::Cell;
    using zephyr::math::geom::VerticesList;

    for (int ic = 0; ic < cells.size(); ++ic) {
        using Vector3d;

        Vector3d c = centers[ic];

        VerticesList vlist;
        if (c.x - D < m_xmin) {
            /// LEFT
            vlist = {
                    Vector3d(c.x + D, c.y, 0.0),
                    Vector3d(c.x + d, c.y + h, 0.0),
                    Vector3d(m_xmin, c.y + h, 0.0),
                    Vector3d(m_xmin, c.y - h, 0.0),
                    Vector3d(c.x + d, c.y - h, 0.0)
            };
        } else if (c.x + D > m_xmax) {
            /// RIGHT
            if (c.y - 0.5 * h < m_ymin) {
                /// RIGHT TOP
                vlist = {
                        Vector3d(m_xmax, m_ymax, 0.0),
                        Vector3d(c.x - D, m_ymax, 0.0),
                        Vector3d(c.x - d, m_ymax - h, 0.0),
                        Vector3d(m_xmax, m_ymax - h, 0.0)
                };
            } else if (c.y + 0.5 * h > m_ymax) {
                /// RIGHT BOTTOM
                vlist = {
                        Vector3d(m_xmax, m_ymin, 0.0),
                        Vector3d(m_xmax, m_ymin + h, 0.0),
                        Vector3d(c.x - d, m_ymin + h, 0.0),
                        Vector3d(c.x - D, m_ymin, 0.0)
                };
            } else {
                /// ORDINARY RIGHT
                vlist = {
                        Vector3d(m_xmax, c.y + h, 0.0),
                        Vector3d(c.x - d, c.y + h, 0.0),
                        Vector3d(c.x - D, c.y, 0.0),
                        Vector3d(c.x - d, c.y - h, 0.0),
                        Vector3d(m_xmax, c.y - h, 0.0)
                };
            }
        } else if (c.y - 0.5 * h < m_ymin) {
            /// BOTTOM
            vlist = {
                    Vector3d(c.x + D, c.y, 0.0),
                    Vector3d(c.x + d, c.y + h, 0.0),
                    Vector3d(c.x - d, c.y + h, 0.0),
                    Vector3d(c.x - D, c.y, 0.0)
            };
        } else if (c.y + 0.5 * h > m_ymax) {
            /// TOP
            vlist = {
                    Vector3d(c.x + D, c.y, 0.0),
                    Vector3d(c.x - D, c.y, 0.0),
                    Vector3d(c.x - d, c.y - h, 0.0),
                    Vector3d(c.x + d, c.y - h, 0.0)
            };
        } else {
            /// INNER CELL
            vlist = {
                    Vector3d(c.x + D, c.y, 0.0),
                    Vector3d(c.x + d, c.y + h, 0.0),
                    Vector3d(c.x - d, c.y + h, 0.0),
                    Vector3d(c.x - D, c.y, 0.0),
                    Vector3d(c.x - d, c.y - h, 0.0),
                    Vector3d(c.x + d, c.y - h, 0.0)
            };
        }


        Cell g_cell(vlist);

        // Граничные условия
        for (auto& face: g_cell.faces) {
            if (face.is_undefined()) {
                continue;
            }

            Vector3d v1 = (Vector3d &) g_cell.vertices.list[face.vertices[0]];
            Vector3d v2 = (Vector3d &) g_cell.vertices.list[face.vertices[1]];

            Vector3d fc = (v1 + v2) / 2.0;

            FaceFlag flag = FaceFlag::ORDINARY;
            if (fc.x + 1.0e-8 * DX > m_xmax) {
                flag = m_right_flag;
            } else if (fc.x - 1.0e-8 * DX < m_xmin) {
                flag = m_left_flag;
            } else if (fc.y + 1.0e-8 * DY > m_ymax) {
                flag = m_top_flag;
            } else if (fc.y - 1.0e-8 * DY < m_ymin) {
                flag = m_bottom_flag;
            }

            face.boundary = flag;
        }

        auto cell = cells[ic];

        cell[coords]   = g_cell.coords;
        cell[size]     = g_cell.size;
        cell[faces]    = g_cell.faces;
        cell[vertices] = g_cell.vertices;

        cell[element].dimension = 2;
        cell[element].kind = kind::EULER;
        cell[neibsSearchRadius].value = 3.0 * D;
    }
     */
}

} // generator
} // mesh
} // zephyr