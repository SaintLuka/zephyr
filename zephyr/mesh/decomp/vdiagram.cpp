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
    : domain_(domain), actual_(false) {

    auto gen = math::Random::Uniform(domain, 13);

    weights_.resize(size, 0.0);
    coords_.resize(size);
    for (int i = 0; i < size; ++i) {
        coords_[i] = gen->next();
    }
    build();
}

VDiagram::VDiagram(const Box& domain, const std::vector<Vector3d>& gs)
    : domain_(domain), coords_(gs), actual_(false) {
    weights_.resize(gs.size(), 0.0);
    build();
}

void VDiagram::normalize() {
    double avg = 0.0;
    for (auto &w: weights_) {
        avg += w;
    }
    avg /= weights_.size();

    for (auto &w: weights_) {
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
    if (actual_) {
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
    double eps = 1.0e-3 * domain_.diameter();

    lines_.resize(n_gen);
    std::vector<double> phi(M, 0.0);
    for (int j = 0; j < M; ++j) {
        phi[j] = M_PI * (2.0 * j - M + 1.0) / (M - 1.0);
    }

    std::vector<Vector3d> inside(M, Vector3d::Zero());
    std::vector<Vector3d> outside(M, Vector3d::Zero());
    for (int iGen = 0; iGen < n_gen; ++iGen) {
        double x_c = coords_[iGen].x();
        double y_c = coords_[iGen].y();

        // Полярные углы углов области
        double phi1 = std::atan2(domain_.vmin.y() - y_c, domain_.vmin.x() - x_c);
        double phi2 = std::atan2(domain_.vmin.y() - y_c, domain_.vmax.x() - x_c);
        double phi3 = std::atan2(domain_.vmax.y() - y_c, domain_.vmax.x() - x_c);
        double phi4 = std::atan2(domain_.vmax.y() - y_c, domain_.vmin.x() - x_c);

        for(int j = 0; j < M; ++j) {
            inside[j].x() = x_c;
            inside[j].y() = y_c;

            if (phi1 <= phi[j] && phi[j] < phi2) {
                outside[j].y() = domain_.vmin.y();
                outside[j].x() = x_c + (domain_.vmin.y() - y_c) / std::tan(phi[j]);
            }
            else if (phi2 <= phi[j] && phi[j] < phi3) {
                outside[j].x() = domain_.vmax.x();
                outside[j].y() = y_c + (domain_.vmax.x() - x_c) * std::tan(phi[j]);
            }
            else if (phi3 <= phi[j] && phi[j] < phi4) {
                outside[j].y() = domain_.vmax.y();
                outside[j].x() = x_c + (domain_.vmax.y() - y_c) / std::tan(phi[j]);
            }
            else {
                outside[j].x() = domain_.vmin.x();
                outside[j].y() = y_c + (domain_.vmin.x() - x_c) * std::tan(phi[j]);
            }

            domain_.shove_in( inside[iGen]);
            domain_.shove_in(outside[iGen]);
        }

        lines_[iGen].resize(M);
        for (int j = 0; j < M; ++j) {
            double diff = domain_.diameter();
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
            lines_[iGen][j] = 0.5 * (inside[j] + outside[j]);
        }
    }

    // ========================================================================
    //  Устанавливаем центры ячеек
    // ========================================================================

    centers_.resize(n_gen);
    for (int i = 0; i < n_gen; ++i) {
        Vector3d c = Vector3d::Zero();
        for (auto &p: lines_[i]) {
            c += p;
        }
        centers_[i] = c / lines_[i].size();
    }

    // ========================================================================
    //  Устанавливаем радиус поиска
    // ========================================================================

    search_radii_.resize(n_gen);

    enum class RadiusType { AVERAGE, MINIMUM };

    RadiusType type = RadiusType::AVERAGE;

    if (type == RadiusType::AVERAGE) {
        // Среднее расстояние
        for (int iGen = 0; iGen < n_gen; ++iGen) {
            search_radii_[iGen] = 0.0;
            for (auto &p: lines_[iGen]) {
                search_radii_[iGen] += (p - coords_[iGen]).norm();
            }
            search_radii_[iGen] /= lines_[iGen].size();
        }
    }
    else if (type == RadiusType::MINIMUM) {
        // Минимальное расстояние
        for (int iGen = 0; iGen < n_gen; ++iGen) {
            search_radii_[iGen] = std::numeric_limits<double>::max();
            for (auto &p: lines_[iGen]) {
                search_radii_[iGen] = std::min(
                        search_radii_[iGen],
                        (p - coords_[iGen]).norm()
                );
            }
        }
    }

    // ========================================================================
    //  Поиск смежных ячеек
    // ========================================================================

    adjacency_.clear();
    adjacency_.resize(n_gen, {});
    for (int iGen = 0; iGen < n_gen - 1; ++iGen) {
        for (auto& v: lines_[iGen]) {
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
                adjacency_[iGen].insert(jGen);
                adjacency_[jGen].insert(iGen);
            }
        }
    }

    actual_ = true;
}

void VDiagram::changed() {
    actual_ = false;

    lines_.clear();
    adjacency_.clear();
    search_radii_.clear();

    // Обнулить цвета
    std::fill(colors_.begin(), colors_.end(), -1);
}

int VDiagram::size() const {
    return coords_.size();
}

const Vector3d& VDiagram::coords(int idx) const {
    return coords_[idx];
}

std::vector<Vector3d>& VDiagram::centers() {
    return centers_;
}

std::vector<double> VDiagram::coords_x() const {
    return np::get_x(coords_);
}

std::vector<double> VDiagram::coords_y() const {
    return np::get_y(coords_);
}

std::vector<double> VDiagram::coords_z() const {
    return np::get_z(coords_);
}

std::vector<double> VDiagram::centers_x() const {
    return np::get_x(centers_);
}

std::vector<double> VDiagram::centers_y() const {
    return np::get_y(centers_);
}

std::vector<double> VDiagram::centers_z() const {
    return np::get_z(centers_);
}

std::vector<double> VDiagram::weights() const {
    return weights_;
}

std::vector<int> VDiagram::degrees() {
    build();

    std::vector<int> deg(size());
    for (int iGen = 0; iGen < size(); ++iGen) {
        deg[iGen] = adjacency_[iGen].size();
    }
    return deg;
}

int VDiagram::chromatic_number() const {
    std::set<int> cols;
    for (auto& col: colors_) {
        cols.insert(col);
    }
    return cols.size();
}

const std::vector<int>& VDiagram::colors() {
    colors_.resize(size(), -1);
    return colors_;
}

std::vector<std::vector<double>> VDiagram::lines_x() {
    build();

    std::vector<std::vector<double>> lx(lines_.size());
    for (int i = 0; i < lines_.size(); ++i) {
        lx[i].resize(lines_[i].size());
        for (int j = 0; j < lines_[i].size(); ++j) {
            lx[i][j] = lines_[i][j].x();
        }
    }
    return lx;
}

std::vector<std::vector<double>> VDiagram::lines_y() {
    build();

    std::vector<std::vector<double>> ly(lines_.size());
    for (int i = 0; i < lines_.size(); ++i) {
        ly[i].resize(lines_[i].size());
        for (int j = 0; j < lines_[i].size(); ++j) {
            ly[i][j] = lines_[i][j].y();
        }
    }
    return ly;
}

std::vector<std::vector<double>> VDiagram::connections_x() {
    build();

    std::vector<std::vector<double>> lx;
    for (int iGen = 0; iGen < size(); ++iGen) {
        for (auto jGen: adjacency_[iGen]) {
            if (jGen > iGen) {
                std::vector<double> segment = {
                        coords_[iGen].x(), coords_[jGen].x()
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
        for (auto jGen: adjacency_[iGen]) {
            if (jGen > iGen) {
                std::vector<double> segment = {
                        coords_[iGen].y(), coords_[jGen].y()
                };
                ly.push_back(segment);
            }
        }
    }
    return ly;
}

const std::vector<double>& VDiagram::search_radii() {
    build();
    return search_radii_;
}

double VDiagram::search_radius(int iGen) const {
    return search_radii_[iGen];
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
            arr_x[iGen][i] = coords_[iGen].x() +
                             search_radii_[iGen] * cos(2 * M_PI * double(i) / M);
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
            arr_y[iGen][i] = coords_[iGen].y() +
                             search_radii_[iGen] * sin(2 * M_PI * double(i) / M);
        }
    }
    return arr_y;
}

void VDiagram::add_generator(double x, double y, double w) {
    x = std::max(domain_.vmin.x(), std::min(x, domain_.vmax.x()));
    y = std::max(domain_.vmin.y(), std::min(y, domain_.vmax.y()));
    double z = 0.0;

    coords_.emplace_back(Vector3d{x, y, z});
    weights_.push_back(w);
    normalize();
    changed();
}

void VDiagram::set_coords(int iGen, double x, double y) {
    coords_.at(iGen) = {x, y, 0.0};
    // Если здесь вызывать changed(), то ломается балансировка
}

void VDiagram::set_coords(int iGen, const Vector3d& p) {
    coords_.at(iGen) = domain_.shove_in(p);
    // Если здесь вызывать changed(), то ломается балансировка
}

void VDiagram::set_weight(int iGen, double w) {
    if (iGen < 0 || iGen >= weights_.size()) {
        throw std::out_of_range("VDiagram::set_weight");
    }
    weights_[iGen] = w;
    // Если здесь вызывать changed(), то ломается балансировка
}

double VDiagram::get_weight(int iGen) const {
    return weights_.at(iGen);
}

double VDiagram::get_coord_x(int iGen) const {
    return coords_.at(iGen).x();
}

double VDiagram::get_coord_y(int iGen) const {
    return coords_.at(iGen).y();
}

double VDiagram::get_coord_z(int iGen) const {
    return coords_.at(iGen).z();
}

Vector3d VDiagram::get_coord(int iGen) const {
    return coords_.at(iGen);
}

void VDiagram::set_coords(const std::vector<Vector3d>& coords) {
    coords_ = coords;
    weights_.resize(coords_.size());
    changed();
}

void VDiagram::set_weights(const std::vector<double>& ws) {
    weights_ = ws;
    coords_.resize(ws.size());
    normalize();
    changed();
}

void VDiagram::paint() {
    if (coords_.empty())
        return;

    colors_.resize(coords_.size(), -1);
    if (colors_[0] >= 0) {
        return;
    }

    // Раскраска графа
    build();

    int N = size();

    // Составляем порядок обхода вершин
    std::vector<int> order;

    std::map<int, std::set<int>> adj;
    for (int iGen = 0; iGen < N; ++iGen) {
        adj[iGen] = adjacency_[iGen];
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

    std::vector full_set(size(), true);

    std::ranges::reverse(order);
    for (auto iGen: order) {
        for (auto& jGen: adjacency_[iGen]) {
            if (colors_[jGen] >= 0) {
                full_set[colors_[jGen]] = false;
            }
        }

        int col = std::distance(full_set.begin(), std::find(full_set.begin(), full_set.end(), true));
        colors_[iGen] = col;

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
    return wdistance(p, coords_[iGen], weights_[iGen]);
}

double VDiagram::distance_gen(int i, int j) {
    return (coords_[i] - coords_[j]).norm();
}

inline double imb_func(double Imax) {
    // Значение дисбаланса, при котором возвращается 1/2
    // Приемлемое значение дисбаланса
    const double I0 = 1.0e-2;
    return Imax / (Imax + I0);
}

void VDiagram::balancing() {
    build();

    std::vector<double> ws(size());
    for (int iGen = 0; iGen < size(); ++iGen) {
        geom::Polygon poly(lines_[iGen]);
        poly.sort();
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

        for (auto jGen: adjacency_[iGen]) {
            if (iGen == jGen) {
                continue;
            }
            double imb = (loads[jGen] - loads[iGen]) / (loads[jGen] + loads[iGen]);

            // Проверка на NaN, возникает при loads[iGen] = loads[jGen] = 0.0
            if (imb != imb) {
                imb = 0.0;
            }

            Vector3d dir = (coords_[jGen] - coords_[iGen]).normalized();

            dr += imb * dir;
            dw += imb;

            min_imb = std::min(min_imb, std::abs(imb));
            max_imb = std::max(max_imb, std::abs(imb));
        }

        // Обновляем координаты

        // Минимальное расстояние до границы ячейки
        double DR = search_radii_[iGen];

        double theta_x = xi_x * imb_func(max_imb);

        Vector3d new_coords = coords_[iGen] + theta_x * DR * dr;

        // Сдвиг к центру масс
        new_coords = sigma * centers_[iGen] + (1.0 - sigma) * new_coords;

        coords_[iGen] = new_coords;

        if (!weighted) {
            continue;
        }

        // Предельный угол для асимптот гипербол
        const double angle_threshold = 90.0 * M_PI / 180.0;
        const double fix = std::cos(0.5 * angle_threshold);

        double min_w = -std::numeric_limits<double>::infinity();
        double max_w = +std::numeric_limits<double>::infinity();
        for (auto jGen: adjacency_[iGen]) {
            min_w = std::max(min_w, weights_[jGen] - fix * distance_gen(iGen, jGen));
            max_w = std::min(max_w, weights_[jGen] + fix * distance_gen(iGen, jGen));
        }

        if (min_w > max_w) {
            double avg = 0.5 * (min_w + max_w);
            min_w = max_w = avg;
        }

        double DW = dw > 0.0 ?
                    std::max(0.0, max_w - weights_[iGen]) :
                    std::max(0.0, weights_[iGen] - min_w);

        double theta_w = xi_w * imb_func(max_imb);

        double new_weight = weights_[iGen] + theta_w * DW * dw;

        new_weight = std::min(max_w, std::max(min_w, new_weight));

        set_weight(iGen, new_weight);
    }

    normalize();
    changed();
}

} // namespace zephyr::mesh::decomp