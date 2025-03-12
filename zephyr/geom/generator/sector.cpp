#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/utils/json.h>
#include <zephyr/geom/generator/sector.h>

namespace zephyr::geom::generator {

Sector::Sector(const Json& config)
        : Generator("sector"), m_r2(-1.0), m_r1(-1.0), m_angle(-1.0), m_hole(false) {

    if (!config["geometry"]) {
        throw std::runtime_error("Sector config doesn't contain key 'geometry'");
    }

    m_r1 = config["geometry"]["r_min"].as<double>();
    m_r2 = config["geometry"]["r_max"].as<double>();

    m_angle = config["geometry"]["angle"].as<double>();

    m_hole = false;
    if (config["geometry"]["hole"]) {
        m_hole = config["geometry"]["hole"].as<bool>();
    }

    init_params();

    if (!config["bounds"]) {
        throw std::runtime_error("Sector config doesn't contain 'bounds'");
    }

    m_bounds.left = Boundary::UNDEFINED;
    if (config["bounds"]["left"]) {
        m_bounds.left = boundary_from_string(config["bounds"]["left"].as<std::string>());
    }

    m_bounds.right = Boundary::UNDEFINED;
    if (config["bounds"]["right"]) {
        m_bounds.right = boundary_from_string(config["bounds"]["right"].as<std::string>());
    }

    m_bounds.inner = Boundary::UNDEFINED;
    if (config["boundary"]["bottom"]) {
        m_bounds.inner = boundary_from_string(config["bounds"]["bottom"].as<std::string>());
    }
    if (config["boundary"]["inner"]) {
        m_bounds.inner = boundary_from_string(config["bounds"]["inner"].as<std::string>());
    }

    m_bounds.outer = Boundary::UNDEFINED;
    if (config["boundary"]["bottom"]) {
        m_bounds.outer = boundary_from_string(config["bounds"]["top"].as<std::string>());
    }
    if (config["boundary"]["inner"]) {
        m_bounds.outer = boundary_from_string(config["bounds"]["outer"].as<std::string>());
    }

    set_n_phi(config["n_phi"].as<int>());
}

Sector::Sector(double r_max, double r_min, double angle, bool hole) :
        Generator("sector"),
        m_r2(r_max),
        m_r1(r_min),
        m_angle(angle),
        m_hole(hole),
        m_bounds() {
    check_params();
    init_params();
}

int Sector::size() const {
    return m_size;
}

Box Sector::bbox() const {
    Vector3d vmin(0.0, 0.0, 0.0);
    Vector3d vmax(0.0, 0.0, 0.0);

    if (m_angle < M_PI) {
        vmin.x() = -m_r2 * std::sin(0.5 * m_angle);
        vmax.x() = +m_r2 * std::sin(0.5 * m_angle);

        vmin.y() = m_hole ? m_r1 * std::cos(0.5 * m_angle) : 0.0;
        vmax.y() = m_r2;
    }
    else {
        vmin.x() = -m_r2;
        vmax.x() = +m_r2;

        vmax.y() = m_r2;
        vmin.y() = -std::cos(0.5 * m_angle);
    }

    return {vmin, vmax};
}

void Sector::init_params() {
    m_angle = std::max(0.0, std::min(m_angle, 2.0 * M_PI));
    m_parts = std::max(1, static_cast<int>(std::round(2.0 * m_angle / M_PI)));
    m_part_angle = m_angle / m_parts;
    m_part_cos = std::cos(m_part_angle);
    m_part_sin = std::sin(m_part_angle);

    fit_angle();

    if (!m_hole && m_angle < 0.1 * M_PI) {
        std::cerr << "Sector constructor warning: small angle, it's recommended to set option 'hole = true'.\n";
    }

    set_n_phi(m_parts * 4);
}

void Sector::fit_angle() {
    if (std::fabs(m_part_angle - M_PI_2) < 1.0e-4) {
        m_part_angle = M_PI_2;
        m_part_cos = 0.0;
        m_part_sin = 1.0;
    }
    m_angle = m_parts * m_part_angle;
}

void Sector::set_n_phi(int N) {
    if (m_hole) {
        m_nx = 0;
        m_Nx = N > 0 ? N : 1;
    }
    else {
        m_nx = static_cast<int>(std::round(double(N) / 2.0 / m_parts));
        m_Nx = 2 * m_parts * m_nx;
        if (N % (2 * m_parts) != 0) {
            std::cerr << "Sector::set_n_phi warning.\n"
                      << "\tN_phi should be a multiple of " << 2 * m_parts << "\n"
                      << "\tThe number has been rounded to " << m_Nx << "\n";
        }
    }

    double Ny = std::log(m_r2 / m_r1) / std::log(1.0 + m_angle / m_Nx);
    m_Ny = static_cast<int>(std::ceil(Ny));

    if (m_hole) {
        // Сектор с отверстием
        m_nx = 0;
        m_ny = 0;

        m_zeta = std::log(m_r2 / m_r1) / std::log(1.0 + m_angle / m_Nx) / m_Ny;
    }
    else {
        // Полный сектор
        m_zeta = 1.0;

        // Подгон внутреннего радиуса
        m_r1 = m_r2 / std::pow(1.0 + m_angle / m_Nx, m_Ny);

        // Пока я не придумал как выбрать оптимальное значение
        m_ny = m_nx;
    }

    m_size = m_parts * m_nx * (m_nx + 2 * m_ny) + m_Nx * m_Ny;
}

void Sector::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}

void Sector::set_axial(bool axial) {
    m_axial = axial;
}

void Sector::check_params() const {
    if (m_r2 <= 0.0 || m_r1 <= 0.0) {
        std::string message = "Sector generator error #1. r_max < 0.0 or r_min <= 0.0";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    if (m_r1 >= m_r2) {
        std::string message = "Sector generator error #2. r_max < r_min";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
}

inline Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

Vector3d Sector::rotate(const Vector3d &v, int r) const {
    double sin = m_part_sin;
    double cos = m_part_cos;

    double SIN, COS;
    switch (r) {
        case 1:
            SIN = sin;
            COS = cos;
            break;

        case 2:
            SIN = 2.0 * sin * cos;
            COS = cos * cos - sin * sin;
            break;

        case 3:
            SIN = sin * (3.0 * cos * cos - sin * sin);
            COS = cos * (cos * cos - 3.0 * sin * sin);
            break;

        default:
            SIN = 0.0;
            COS = 1.0;
            break;
    }

    return {
            COS * v.x() - SIN * v.y(),
            SIN * v.x() + COS * v.y(),
            0.0
    };
}

Grid Sector::make() {
    // Вершины центральной части
    std::vector<std::vector<GNode::Ptr>> CV;
    std::vector<std::vector<GNode::Ptr>> CV2;
    std::vector<std::vector<GNode::Ptr>> CV3;
    std::vector<std::vector<GNode::Ptr>> CV4;
    if (!m_hole) {
        CV.resize(m_nx + 1, std::vector<GNode::Ptr>(m_nx + 1, nullptr));
        if (m_parts > 1) {
            CV2.resize(m_nx + 1, std::vector<GNode::Ptr>(m_nx + 1, nullptr));
            if (m_parts > 2) {
                CV3.resize(m_nx + 1, std::vector<GNode::Ptr>(m_nx + 1, nullptr));
                if (m_parts > 3) {
                    CV4.resize(m_nx + 1, std::vector<GNode::Ptr>(m_nx + 1, nullptr));
                }
            }
        }
    }

    // Вершины секторальной части
    std::vector<std::vector<GNode::Ptr>> SV;
    SV.resize(m_Nx + 1, std::vector<GNode::Ptr>(m_ny + m_Ny + 1, nullptr));

    if (!m_hole) {
        double m_r0 = m_r1 / std::pow(1.0 + M_PI_4 / m_nx, m_ny);

        // Вершины базового квадрата
        Vector3d v0 = Vector3d::Zero();
        Vector3d v1 = to_cartesian(m_r0, 0.0);
        Vector3d v2 = to_cartesian(m_r0, M_PI_4);
        Vector3d v3 = to_cartesian(m_r0, M_PI_2);

        // Вектора граней
        Vector3d a = v1 - v0;
        Vector3d b = v3 - v0;
        Vector3d c = v2 - v3;
        Vector3d d = v2 - v1;

        Vector3d skew = ((c + d) - (a + b)) / 2;

        for (size_t i = 0; i <= m_nx; ++i) {
            double x = double(i) / m_nx;
            for (size_t j = 0; j <= m_nx; ++j) {
                double y = double(j) / m_nx;
                CV[i][j] = GNode::create(v0 + a * x + b * y + skew * x * y);
            }
        }

        // Создаем вершины секторальной части
        double d_alpha = M_PI_4 / m_nx;
        for (size_t i = 0; i <= 2 * m_nx; ++i) {
            if (i <= m_nx) {
                SV[i][0] = CV[m_nx][i];
            } else {
                SV[i][0] = CV[2 * m_nx - i][m_nx];
            }

            double phi = (M_PI_4 * i) / m_nx;
            for (size_t j = 1; j < m_ny; ++j) {
                double r = m_r0 * pow(1.0 + d_alpha, j);
                SV[i][j] = GNode::create(to_cartesian(r, phi));
            }
            SV[i][m_ny] = GNode::create(to_cartesian(m_r1, phi));
        }


        // ========================================================================
        //     Сглаживание сетки
        // ========================================================================

        int smooth_iter = 2000; // число итераций сглаживания

        std::vector<GNode::Ptr> T(m_nx + 1, nullptr);
        std::vector<GNode::Ptr> R(m_nx + 1, nullptr);
        for (size_t i = 0; i <= m_nx; ++i) {
            R[i] = SV[i][1];
            T[i] = SV[2 * m_nx - i][1];
        }

        for (int k = 0; k < smooth_iter; ++k) {
            Vector3d b1 = {0.0, 0.0, 0.0};
            Vector3d b2 = {0.0, 0.0, 0.0};
            Vector3d b3 = {0.0, 0.0, 0.0};
            Vector3d b4 = {0.0, 0.0, 0.0};

            // Сглаживание секторальной части
            for (size_t j = 1; j < m_ny; ++j) {
                b1 = SV[0][j - 1]->v;
                b2 = SV[0][j + 1]->v;
                b3 = SV[1][j]->v;

                SV[0][j]->v.x() = (b1.x() + b2.x() + 2.0 * b3.x()) / 4.0;
                SV[0][j]->v.y() = 0.0;

                for (size_t i = 1; i < 2 * m_nx; ++i) {
                    b1 = SV[i + 1][j]->v;
                    b2 = SV[i - 1][j]->v;
                    b3 = SV[i][j - 1]->v;
                    b4 = SV[i][j + 1]->v;

                    SV[i][j]->v = (b1 + b2 + b3 + b4) / 4.0;
                }

                b1 = SV[2 * m_nx][j - 1]->v;
                b2 = SV[2 * m_nx][j + 1]->v;
                b3 = SV[2 * m_nx - 1][j]->v;

                SV[2 * m_nx][j]->v.x() = 0.0;
                SV[2 * m_nx][j]->v.y() = (b1.y() + b2.y() + 2.0 * b3.y()) / 4.0;
            }

            // Сглаживание центральной части
            for (size_t j = 1; j <= m_nx; ++j) {
                b1 = CV[0][j - 1]->v;
                b2 = j < m_nx ? CV[0][j + 1]->v : T[0]->v;
                b3 = CV[1][j]->v;

                CV[0][j]->v.x() = 0.0;
                CV[0][j]->v.y() = (b1.y() + b2.y() + 2.0 * b3.y()) / 4.0;
            }

            for (size_t i = 1; i <= m_nx; ++i) {
                for (size_t j = 1; j <= m_nx; ++j) {
                    b1 = i < m_nx ? CV[i + 1][j]->v : R[j]->v;
                    b2 = CV[i - 1][j]->v;
                    b3 = CV[i][j - 1]->v;
                    b4 = j < m_nx ? CV[i][j + 1]->v : T[i]->v;

                    if (R[j] == T[i]) {
                        // особая точка
                        CV[i][j]->v = (b1 + b2 + b3) / 3.0;
                    }
                    else {
                        CV[i][j]->v = (b1 + b2 + b3 + b4) / 4.0;
                    }
                }
            }

            for (size_t i = 1; i <= m_nx; ++i) {
                b1 = CV[i - 1][0]->v;
                b2 = i < m_nx ? CV[i + 1][0]->v : R[0]->v;
                b3 = CV[i][1]->v;

                CV[i][0]->v.x() = (b1.x() + b2.x() + 2.0 * b3.x()) / 4.0;
                CV[i][0]->v.y() = 0.0;
            }
        }


        // ========================================================================
        //     Преобразуем внутреннюю часть
        // ========================================================================

        if (std::abs(m_part_angle - 0.5 * M_PI) > 1.0e-4) {

            double xi = 2.0 * m_part_angle / M_PI;

            auto transform = [this, xi](GNode::Ptr &v) {
                if (!v) { return; }
                double r = std::sqrt(v->v.x() * v->v.x() + v->v.y() * v->v.y());
                double phi = std::atan2(v->v.y(), v->v.x());

                double r2 = m_r1 * std::pow(r / m_r1, xi);
                double phi2 = phi * xi;

                v->v.x() = r2 * std::cos(phi2);
                v->v.y() = r2 * std::sin(phi2);
            };

            for (auto &series: CV) {
                for (auto &v: series) {
                    if (v.use_count() > 1) {
                        continue;
                    }
                    transform(v);
                }
            }

            for (auto &series: SV) {
                for (auto &v: series) {
                    transform(v);
                }
            }
        }


        // ========================================================================
        //     Симметризация
        // ========================================================================

        auto symmetrization = [this](GNode::Ref p1, GNode::Ref p2) {
            Vector3d v1 = p1->v;
            Vector3d v2 = p2->v;

            p1->v.x() = 0.5 * (v1.x() + m_part_cos * v2.x() + m_part_sin * v2.y());
            p1->v.y() = 0.5 * (v1.y() + m_part_sin * v2.x() - m_part_cos * v2.y());

            p2->v.x() = 0.5 * (v2.x() + m_part_cos * v1.x() + m_part_sin * v1.y());
            p2->v.y() = 0.5 * (v2.y() + m_part_sin * v1.x() - m_part_cos * v1.y());
        };

        for (size_t i = 0; i <= m_nx; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                // Внутри квадратных блоков
                symmetrization(CV[i][j], CV[j][i]);
            }

            for (size_t j = 1; j <= m_ny; ++j) {
                // Внутри сектора
                symmetrization(SV[i][j], SV[2 * m_nx - i][j]);
            }
        }
    }


    // ========================================================================
    //     Копируем внутреннюю часть
    // ========================================================================

    if (!m_hole) {
        // Секторальная часть
        for (size_t r = 1; r < m_parts; ++r) {
            for (size_t j = 0; j <= m_ny; ++j) {
                for (size_t i = 1; i <= 2 * m_nx; ++i) {
                    SV[2 * m_nx * r + i][j] = GNode::create(rotate(SV[i][j]->v, r));
                }
            }
        }

        // Центральная часть
        for (size_t i = 0; i <= m_nx; ++i) {
            for (size_t j = 0; j <= m_nx; ++j) {
                Vector3d v = CV[i][j]->v;
                if (m_parts > 1) {
                    CV2[m_nx - j][i] = GNode::create(rotate(v, 1));
                    if (m_parts > 2) {
                        CV3[m_nx - i][m_nx - j] = GNode::create(rotate(v, 2));
                        if (m_parts > 3) {
                            CV4[j][m_nx - i] = GNode::create(rotate(v, 3));
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i <= m_nx; ++i) {
            if (m_parts > 1) {
                CV2[m_nx][i] = CV[0][i];
                CV2[0][i] = SV[4*m_nx - i][0];
                CV2[i][m_nx] = SV[3*m_nx - i][0];
            }
            if (m_parts > 2) {
                CV3[i][m_nx] = CV2[i][0];
                CV3[0][i] = SV[5 * m_nx - i][0];
                CV3[i][0] = SV[5 * m_nx + i][0];
            }
            if (m_parts > 3) {
                CV4[0][i] = CV3[m_nx][i];
                CV4[i][0] = SV[6 * m_nx + i][0];
                CV4[m_nx][i] = SV[7 * m_nx + i][0];

                // замыкание круга
                if (std::abs(m_angle - 2 * M_PI) < 1.0e-4) {
                    CV4[i][m_nx] = CV[i][0];
                }
            }
        }
    }

    // Круг не замыкается
    if (!m_hole && std::abs(m_angle - 2 * M_PI) > 1.0e-4) {
        for (size_t i = 0; i <= m_nx; ++i) {
            CV[i][0]->add_boundary(m_bounds.right);
        }
        if (m_parts == 1) {
            for (size_t i = 0; i <= m_nx; ++i) {
                CV[0][i]->add_boundary(m_bounds.left);
            }
        } else if (m_parts == 2) {
            for (size_t i = 0; i <= m_nx; ++i) {
                CV2[i][0]->add_boundary(m_bounds.left);
            }
        } else if (m_parts == 3) {
            for (size_t i = 0; i <= m_nx; ++i) {
                CV3[m_nx][i]->add_boundary(m_bounds.left);
            }
        } else {
            for (size_t i = 0; i <= m_nx; ++i) {
                CV4[i][m_nx]->add_boundary(m_bounds.left);
            }
        }
    }




    // ========================================================================
    //     Вершины с внешней части сектора
    // ========================================================================

    double D_alpha = m_angle / m_Nx;
    for (size_t i = 0; i <= m_Nx; ++i) {
        SV[i].resize(m_ny + m_Ny + 1);

        double phi = m_angle * double(i) / m_Nx;
        SV[i][m_ny] = GNode::create(to_cartesian(m_r1, phi));
        for (size_t j = 1; j < m_Ny; ++j) {
            double r = m_r1 * std::pow(1.0 + D_alpha, m_zeta * j);

            SV[i][m_ny + j] = GNode::create(to_cartesian(r, phi));
        }
        SV[i][m_ny + m_Ny] = GNode::create(to_cartesian(m_r2, phi));
    }
    if (m_hole) {
        for (size_t i = 0; i <= m_Nx; ++i) {
            SV[i][0]->add_boundary(m_bounds.inner);
        }
    }
    for (size_t i = 0; i <= m_Nx; ++i) {
        SV[i][m_ny + m_Ny]->add_boundary(m_bounds.outer);
    }

    if (std::abs(m_angle - 2 * M_PI) < 1.0e-4) {
        for (size_t j = 0; j <= m_ny + m_Ny; ++j) {
            SV[m_Nx][j] = SV[0][j];
        }
    }
    else {
        for (size_t j = 0; j <= m_ny + m_Ny; ++j) {
            SV[0][j]->add_boundary(m_bounds.right);
            SV[m_Nx][j]->add_boundary(m_bounds.left);
        }
    }


    // ========================================================================
    //     Генерируем итоговый массив вершин
    // ========================================================================

    auto get_center_vface = [&](size_t x, size_t y, size_t r)
            -> std::array<GNode::Ptr, 2> {

        assert(0 <= x && x <= m_nx);
        assert(0 <= y && y <= m_nx);

        switch (r) {
            case 1:
                return {CV2[x][y], CV2[x][y + 1]};
            case 2:
                return {CV3[x][y], CV3[x][y + 1]};
            case 3:
                return {CV4[x][y], CV4[x][y + 1]};
            default:
                return {CV[x][y], CV[x][y + 1]};
        }
    };

    auto get_sector_vface = [this, &SV](size_t x, size_t y)
            -> std::array<GNode::Ptr, 2> {

        return {SV[x][y], SV[x][y + 1]};
    };


    int max_nodes = (m_Nx + 1)*(m_ny + m_Ny + 1) + m_parts * (m_nx + 1) * (m_nx + 1);

    Grid grid;
    grid.set_axial(m_axial);
    grid.reserve_cells(m_size);
    grid.reserve_nodes(max_nodes);
    int n_cells = 0;
    int n_nodes = 0;

    // Создаем вершины центральной части
    for (size_t r = 0; r < m_parts; ++r) {
        for (size_t x = 0; x < m_nx; ++x) {
            for (size_t y = 0; y < m_nx; ++y) {
                auto l_face = get_center_vface(x, y, r);
                auto r_face = get_center_vface(x + 1, y, r);

                std::array<GNode::Ptr, 4> vlist = {
                        l_face[0],
                        r_face[0],
                        r_face[1],
                        l_face[1],
                };

                auto cell = GCell::quad(vlist);
                cell.index = n_cells;
                grid += cell;
                ++n_cells;

                for (auto&v : vlist) {
                    if (v->index < 0) {
                        v->index = n_nodes;
                        grid += v;
                        ++n_nodes;
                    }
                }
            }
        }
    }

    // Создаем вершины секторальной части
    for (size_t x = 0; x < m_Nx; ++x) {
        for (size_t y = 0; y < m_ny + m_Ny; ++y) {
            auto l_face = get_sector_vface(x, y);
            auto r_face = get_sector_vface(x + 1, y);

            std::array<GNode::Ptr, 4> vlist = {
                    r_face[0],
                    l_face[0],
                    l_face[1],
                    r_face[1]
            };

            auto cell = GCell::quad(vlist);
            cell.index = n_cells;
            grid += cell;
            ++n_cells;

            for (auto&v : vlist) {
                if (v->index < 0) {
                    v->index = n_nodes;
                    grid += v;
                    ++n_nodes;
                }
            }
        }
    }

    grid.setup_adjacency();
    return grid;
}

} // namespace zephyr::geom::generator