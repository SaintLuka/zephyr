#include <random>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>

#include <zephyr/mesh/decomp/vdiagram.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/math/random.h>
#include <zephyr/utils/numpy.h>

namespace zephyr::mesh::decomp {

VDiagram::VDiagram(const Box& domain, int size)
    : m_domain(domain), m_actual(false) {
    auto gen = domain.random2D(13);

    m_weights.resize(size, 0.0);
    m_coords.resize(size);
    for (int i = 0; i < size; ++i) {
        m_coords[i] = gen.get();
    }
    build();
}

VDiagram::VDiagram(const Box& domain, const std::vector<Vector3d>& gs)
    : m_domain(domain), m_coords(gs), m_actual(false) {
    m_weights.resize(gs.size(), 0.0);
    build();
}

void VDiagram::normalize() {
    double avg = 0.0;
    for (auto &w: m_weights) {
        avg += w;
    }
    avg /= m_weights.size();

    for (auto &w: m_weights) {
        w -= avg;
    }
}

double VDiagram::edge_function(const Vector3d& p, int iGen) const {
    double min_dist = std::numeric_limits<double>::max();

    int N = size();
    for (int j = 0; j < N; ++j) {
        if (j != iGen) {
            double dist = wdistance(p, j);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
    }
    return (min_dist - wdistance(p, iGen));
}

void VDiagram::build() {
    if (m_actual) {
        return;
    }

    // Число ячеек
    int n_gen = size();

    // ========================================================================
    //  Находим границы ячеек
    // ========================================================================

    // Точек на границе ячейки
    int M = 200;

    // Точность позиционирования границы
    double eps = 1.0e-3 * m_domain.diameter();

    m_lines.resize(n_gen);
    std::vector<double> phi(M, 0.0);
    for (int j = 0; j < M; ++j) {
        phi[j] = M_PI * (2.0 * j - M + 1.0) / (M - 1.0);
    }

    std::vector<Vector3d> inside(M, Vector3d::Zero());
    std::vector<Vector3d> outside(M, Vector3d::Zero());
    for (int iGen = 0; iGen < n_gen; ++iGen) {
        double x_c = m_coords[iGen].x();
        double y_c = m_coords[iGen].y();

        // Полярные углы углов области
        double phi1 = std::atan2(m_domain.vmin.y() - y_c, m_domain.vmin.x() - x_c);
        double phi2 = std::atan2(m_domain.vmin.y() - y_c, m_domain.vmax.x() - x_c);
        double phi3 = std::atan2(m_domain.vmax.y() - y_c, m_domain.vmax.x() - x_c);
        double phi4 = std::atan2(m_domain.vmax.y() - y_c, m_domain.vmin.x() - x_c);

        for(int j = 0; j < M; ++j) {
            inside[j].x() = x_c;
            inside[j].y() = y_c;

            if (phi1 <= phi[j] && phi[j] < phi2) {
                outside[j].y() = m_domain.vmin.y();
                outside[j].x() = x_c + (m_domain.vmin.y() - y_c) / std::tan(phi[j]);
            }
            else if (phi2 <= phi[j] && phi[j] < phi3) {
                outside[j].x() = m_domain.vmax.x();
                outside[j].y() = y_c + (m_domain.vmax.x() - x_c) * std::tan(phi[j]);
            }
            else if (phi3 <= phi[j] && phi[j] < phi4) {
                outside[j].y() = m_domain.vmax.y();
                outside[j].x() = x_c + (m_domain.vmax.y() - y_c) / std::tan(phi[j]);
            }
            else {
                outside[j].x() = m_domain.vmin.x();
                outside[j].y() = y_c + (m_domain.vmin.x() - x_c) * std::tan(phi[j]);
            }

            m_domain.shove_in( inside[iGen]);
            m_domain.shove_in(outside[iGen]);
        }

        m_lines[iGen].resize(M);
        for (int j = 0; j < M; ++j) {
            double diff = m_domain.diameter();
            while (diff >= eps) {
                Vector3d middle = 0.5 * (inside[j] + outside[j]);
                double value = edge_function(middle, iGen);
                if (value >= 0.0) {
                    inside[j] = middle;
                } else {
                    outside[j] = middle;
                }
                diff = (inside[j] - outside[j]).norm();
            }
            m_lines[iGen][j] = 0.5 * (inside[j] + outside[j]);
        }
    }

    // ========================================================================
    //  Устанавливаем центры ячеек
    // ========================================================================

    m_centers.resize(n_gen);
    for (int i = 0; i < n_gen; ++i) {
        Vector3d c = Vector3d::Zero();
        for (auto &p: m_lines[i]) {
            c += p;
        }
        m_centers[i] = c / m_lines[i].size();
    }

    // ========================================================================
    //  Устанавливаем радиус поиска
    // ========================================================================

    m_search_radii.resize(n_gen);

    enum class RadiusType { AVERAGE, MINIMUM };

    RadiusType type = RadiusType::AVERAGE;

    if (type == RadiusType::AVERAGE) {
        // Среднее расстояние
        for (int iGen = 0; iGen < n_gen; ++iGen) {
            m_search_radii[iGen] = 0.0;
            for (auto &p: m_lines[iGen]) {
                m_search_radii[iGen] += (p - m_coords[iGen]).norm();
            }
            m_search_radii[iGen] /= m_lines[iGen].size();
        }
    }
    else if (type == RadiusType::MINIMUM) {
        // Минимальное расстояние
        for (int iGen = 0; iGen < n_gen; ++iGen) {
            m_search_radii[iGen] = std::numeric_limits<double>::max();
            for (auto &p: m_lines[iGen]) {
                m_search_radii[iGen] = std::min(
                        m_search_radii[iGen],
                        (p - m_coords[iGen]).norm()
                );
            }
        }
    }

    // ========================================================================
    //  Поиск смежных ячеек
    // ========================================================================

    m_adjacency.clear();
    m_adjacency.resize(n_gen, {});
    for (int iGen = 0; iGen < n_gen - 1; ++iGen) {
        for (auto& v: m_lines[iGen]) {
            double min = std::numeric_limits<double>::max();
            int jGen = -1;
            for (int kGen = 0; kGen < n_gen; ++kGen) {
                double wdist = wdistance(v, kGen);
                if (wdist < min) {
                    min = wdist;
                    jGen = kGen;
                }
            }
            if (jGen >= 0 && jGen != iGen) {
                m_adjacency[iGen].insert(jGen);
                m_adjacency[jGen].insert(iGen);
            }
        }
    }

    m_actual = true;
}

void VDiagram::changed() {
    m_actual = false;

    m_lines.clear();
    m_adjacency.clear();
    m_search_radii.clear();

    // Обнулить цвета
    std::fill(m_colors.begin(), m_colors.end(), -1);
}

int VDiagram::size() const {
    return m_coords.size();
}

const Vector3d& VDiagram::coords(int idx) const {
    return m_coords[idx];
}

std::vector<Vector3d>& VDiagram::centers() {
    return m_centers;
}

std::vector<double> VDiagram::coords_x() const {
    return np::get_x(m_coords);
}

std::vector<double> VDiagram::coords_y() const {
    return np::get_y(m_coords);
}

std::vector<double> VDiagram::coords_z() const {
    return np::get_z(m_coords);
}

std::vector<double> VDiagram::centers_x() const {
    return np::get_x(m_centers);
}

std::vector<double> VDiagram::centers_y() const {
    return np::get_y(m_centers);
}

std::vector<double> VDiagram::centers_z() const {
    return np::get_z(m_centers);
}

std::vector<double> VDiagram::weights() const {
    return m_weights;
}

std::vector<int> VDiagram::degrees() {
    build();

    std::vector<int> deg(size());
    for (int iGen = 0; iGen < size(); ++iGen) {
        deg[iGen] = m_adjacency[iGen].size();
    }
    return deg;
}

int VDiagram::chromatic_number() const {
    std::set<int> cols;
    for (auto& col: m_colors) {
        cols.insert(col);
    }
    return cols.size();
}

const std::vector<int>& VDiagram::colors() {
    m_colors.resize(size(), -1);
    return m_colors;
}

std::vector<std::vector<double>> VDiagram::lines_x() {
    build();

    std::vector<std::vector<double>> lx(m_lines.size());
    for (int i = 0; i < m_lines.size(); ++i) {
        lx[i].resize(m_lines[i].size());
        for (int j = 0; j < m_lines[i].size(); ++j) {
            lx[i][j] = m_lines[i][j].x();
        }
    }
    return lx;
}

std::vector<std::vector<double>> VDiagram::lines_y() {
    build();

    std::vector<std::vector<double>> ly(m_lines.size());
    for (int i = 0; i < m_lines.size(); ++i) {
        ly[i].resize(m_lines[i].size());
        for (int j = 0; j < m_lines[i].size(); ++j) {
            ly[i][j] = m_lines[i][j].y();
        }
    }
    return ly;
}

std::vector<std::vector<double>> VDiagram::connections_x() {
    build();

    std::vector<std::vector<double>> lx;
    for (int iGen = 0; iGen < size(); ++iGen) {
        for (auto jGen: m_adjacency[iGen]) {
            if (jGen > iGen) {
                std::vector<double> segment = {
                        m_coords[iGen].x(), m_coords[jGen].x()
                };
                lx.push_back(segment);
            }
        }
    }
    return lx;
}

std::vector<std::vector<double>> VDiagram::connections_y() {
    build();

    std::vector<std::vector<double>> ly;
    for (int iGen = 0; iGen < size(); ++iGen) {
        for (auto jGen: m_adjacency[iGen]) {
            if (jGen > iGen) {
                std::vector<double> segment = {
                        m_coords[iGen].y(), m_coords[jGen].y()
                };
                ly.push_back(segment);
            }
        }
    }
    return ly;
}

const std::vector<double>& VDiagram::search_radii() {
    build();
    return m_search_radii;
}

double VDiagram::search_radius(int iGen) const {
    return m_search_radii[iGen];
}

std::vector<std::vector<double>> VDiagram::search_area_x() {
    build();

    const int M = 100;
    auto radii = search_radii();
    std::vector<std::vector<double>> arr_x(
            size(), std::vector<double>(M, 0.0)
    );
    for (int iGen = 0; iGen < size(); ++iGen) {
        for (int i = 0; i < M; ++i) {
            arr_x[iGen][i] = m_coords[iGen].x() +
                             m_search_radii[iGen] * cos(2 * M_PI * double(i) / M);
        }
    }
    return arr_x;
}

std::vector<std::vector<double>> VDiagram::search_area_y() {
    build();

    const int M = 100;
    auto radii = search_radii();
    std::vector<std::vector<double>> arr_y(
            size(), std::vector<double>(M, 0.0)
    );
    for (int iGen = 0; iGen < size(); ++iGen) {
        for (int i = 0; i < M; ++i) {
            arr_y[iGen][i] = m_coords[iGen].y() +
                             m_search_radii[iGen] * sin(2 * M_PI * double(i) / M);
        }
    }
    return arr_y;
}

void VDiagram::add_generator(double x, double y, double w) {
    x = std::max(m_domain.vmin.x(), std::min(x, m_domain.vmax.x()));
    y = std::max(m_domain.vmin.y(), std::min(y, m_domain.vmax.y()));
    double z = 0.0;

    m_coords.emplace_back(Vector3d{x, y, z});
    m_weights.push_back(w);
    normalize();
    changed();
}

void VDiagram::set_coords(int iGen, double x, double y) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::set_coords");
    }
    m_coords[iGen] = {x, y, 0.0};
    // Если здесь вызывать changed(), то ломается балансировка
}

void VDiagram::set_weight(int iGen, double w) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::set_weight");
    }
    m_weights[iGen] = w;
    // Если здесь вызывать changed(), то ломается балансировка
}

double VDiagram::get_weight(int iGen) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::get_weight");
    }
    return m_weights[iGen];
}

double VDiagram::get_coord_x(int iGen) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::get_coord_x");
    }
    return m_coords[iGen].x();
}

double VDiagram::get_coord_y(int iGen) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::get_coord_y");
    }
    return m_coords[iGen].y();
}

double VDiagram::get_coord_z(int iGen) {
    if (iGen < 0 || iGen >= m_weights.size()) {
        throw std::out_of_range("VDiagram::get_coord_y");
    }
    return m_coords[iGen].z();
}

void VDiagram::set_coords(const std::vector<Vector3d>& coords) {
    m_coords = coords;
    m_weights.resize(m_coords.size());
    changed();
}

void VDiagram::set_weights(const std::vector<double>& ws) {
    m_weights = ws;
    m_coords.resize(ws.size());
    normalize();
    changed();
}

void VDiagram::paint() {
    if (m_coords.empty())
        return;

    m_colors.resize(m_coords.size(), -1);
    if (m_colors[0] >= 0) {
        return;
    }

    // Раскраска графа
    build();

    int N = size();

    // Составляем порядок обхода вершин
    std::vector<int> order;

    std::map<int, std::set<int>> adj;
    for (int iGen = 0; iGen < N; ++iGen) {
        adj[iGen] = m_adjacency[iGen];
    }

    while (!adj.empty()) {
        // Находим вершину с минимальной степенью
        int min_gen = adj.begin()->first;
        int min_deg = adj.begin()->second.size();
        for (auto &gen: adj) {
            int deg = gen.second.size();
            if (deg < min_deg) {
                min_deg = deg;
                min_gen = gen.first;
            }
        }

        // Добавляем в список
        order.push_back(min_gen);

        // Исключаем вершину
        adj.erase(min_gen);
        for (auto& gen: adj) {
            gen.second.erase(min_gen);
        }
    }

    std::vector<bool> full_set(size(), true);

    std::reverse(order.begin(), order.end());
    for (auto iGen: order) {
        for (auto& jGen: m_adjacency[iGen]) {
            if (m_colors[jGen] >= 0) {
                full_set[m_colors[jGen]] = false;
            }
        }

        int col = std::distance(full_set.begin(), std::find(full_set.begin(), full_set.end(), true));
        m_colors[iGen] = col;

        std::fill(full_set.begin(), full_set.end(), true);
    }
}

int VDiagram::rank(const Vector3d& v) const {
    int res = -1;
    double min_dist = std::numeric_limits<double>::max();

    for (int iGen = 0; iGen < size(); ++iGen) {
        double dist = wdistance(v, iGen);
        if (dist < min_dist) {
            min_dist = dist;
            res = iGen;
        }
    }
    return res;
}

double VDiagram::wdistance(const Vector3d& p, const Vector3d& g, double w) {
    return (p - g).norm() - w;
}

double VDiagram::wdistance(const Vector3d& p, int iGen) const {
    return wdistance(p, m_coords[iGen], m_weights[iGen]);
}

double VDiagram::distance_gen(int i, int j) {
    return (m_coords[i] - m_coords[j]).norm();
}

double imb_func(double Imax) {
    // Значение дисбаланса, при котором возвращается 1/2
    // Приемлимое значение дисбаланса
    const double I0 = 1.0e-2;
    return Imax / (Imax + I0);
}

void VDiagram::balancing() {
    build();

    std::vector<double> ws(size());
    for (int iGen = 0; iGen < size(); ++iGen) {
        geom::Polygon poly(m_lines[iGen], true);
        ws[iGen] = poly.area();
    }
    balancing(ws);
}

void VDiagram::balancing(const std::vector<double> &loads) {
    build();

    const double xi_x = mobility;
    const double xi_w = growth_rate;
    const double sigma = centroidal;

    bool weighted = xi_w > 0.0;

    for (int iGen = 0; iGen < size(); ++iGen) {
        Vector3d dr = {0.0, 0.0, 0.0};
        double dw = 0.0;

        double min_imb = 1.0e20;
        double max_imb = 0.0;

        for (auto jGen: m_adjacency[iGen]) {
            if (iGen == jGen) {
                continue;
            }
            double imb = (loads[jGen] - loads[iGen]) / (loads[jGen] + loads[iGen]);

            // Проверка на NaN, возникает при loads[iGen] = loads[jGen] = 0.0
            if (imb != imb) {
                imb = 0.0;
            }

            Vector3d dir = (m_coords[jGen] - m_coords[iGen]).normalized();

            dr += imb * dir;
            dw += imb;

            min_imb = std::min(min_imb, std::abs(imb));
            max_imb = std::max(max_imb, std::abs(imb));
        }

        // Обновляем координаты

        // Минимальное расстояние до границы ячейки
        double DR = m_search_radii[iGen];

        double theta_x = xi_x * imb_func(max_imb);

        Vector3d new_coords = m_coords[iGen] + theta_x * DR * dr;

        // Сдвиг к центру масс
        new_coords = sigma * m_centers[iGen] + (1.0 - sigma) * new_coords;

        m_coords[iGen] = new_coords;

        if (!weighted) {
            continue;
        }

        // Предельный угол для асимптот гипербол
        const double angle_threshold = 90.0 * M_PI / 180.0;
        const double fix = std::cos(0.5 * angle_threshold);

        double min_w = -std::numeric_limits<double>::infinity();
        double max_w = +std::numeric_limits<double>::infinity();
        for (auto jGen: m_adjacency[iGen]) {
            min_w = std::max(min_w, m_weights[jGen] - fix * distance_gen(iGen, jGen));
            max_w = std::min(max_w, m_weights[jGen] + fix * distance_gen(iGen, jGen));
        }

        if (min_w > max_w) {
            double avg = 0.5 * (min_w + max_w);
            min_w = max_w = avg;
        }

        double DW = dw > 0.0 ?
                    std::max(0.0, max_w - m_weights[iGen]) :
                    std::max(0.0, m_weights[iGen] - min_w);

        double theta_w = xi_w * imb_func(max_imb);

        double new_weight = m_weights[iGen] + theta_w * DW * dw;

        new_weight = std::min(max_w, std::max(min_w, new_weight));

        set_weight(iGen, new_weight);
    }

    normalize();
    changed();
}

} // namespace zephyr::mesh::decomp