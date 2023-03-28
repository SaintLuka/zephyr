#include <Vector3d>
#include <iostream>
#include <algorithm>

#include <zephyr/network/mpi/network.h>

#include <zephyr/network/decomposition/ORB/blocks.h>

namespace zephyr { namespace network { namespace decomposition {

Blocks::Blocks(const std::string& type, size_t n_proc,
       int n_proc_1, int n_proc_2, int n_proc_3)
       : m_ready(false), m_map(type), m_size(n_proc) {

    switch (m_map.dimension()) {
        case 1:
            init_sizes_1D();
            break;
        case 2:
            init_sizes_2D(n_proc_1);
            break;
        case 3:
            init_sizes_3D(n_proc_1, n_proc_2);
            break;
        default:
            throw std::runtime_error("Странная декомпозиция с размерностью > 3");
    }

    init_limits();
    init_indices();
    init_coords();
    init_widths();
    init_loads();
}

Blocks::Blocks(const std::string& type, size_t n_proc, const Vector3d& sizes)
        : m_ready(false), m_map(type), m_size(n_proc) {

    auto dim = m_map.dimension();

    if (dim == 1) {
        init_sizes_1D();
    }
    else if (dim == 2) {
        Vector3d L = m_map(sizes);
        int n_proc_1 = static_cast<int>(std::round(L.x * std::sqrt(m_size / (L.x * L.y))));
        n_proc_1 = std::min(std::max(1, n_proc_1), int(m_size));
        init_sizes_2D(n_proc_1);
    }
    else if (dim == 3) {
        Vector3d L = m_map(sizes);
        int n_proc_1 = static_cast<int>(std::round(L.x * std::cbrt(m_size / (L.x * L.y * L.z))));
        n_proc_1 = std::min(std::max(1, n_proc_1), int(m_size));

        int n_proc_2 = static_cast<int>(std::round(L.y * std::cbrt(m_size / (L.x * L.y * L.z))));
        n_proc_2 = std::min(std::max(1, n_proc_2), int(m_size));

        init_sizes_3D(n_proc_1, n_proc_2);
    }
    else {
        throw std::runtime_error("Странная декомпозиция с размерностью > 3");
    }

    init_limits();
    init_indices();
    init_coords();
    init_widths();
    init_loads();
}

void Blocks::readiness() const {
    if (!m_ready) {
        std::cerr << "Block Decomposition is not prepared, call Blocks::setup_bounds first\n";
    }
}

Vector3d Blocks::mapping(const Vector3d& vec) const {
    return m_map.mapping(vec);
}

Vector3d Blocks::inverse(const Vector3d& vec) const {
    return m_map.inverse(vec);
}

// Коэффициент расширения области
const double EC = 0.01;

void Blocks::setup_bounds(Vector3d vmin, Vector3d vmax) {
    vmin = m_map(vmin);
    vmax = m_map(vmax);

    int dim = m_map.dimension();
    for (size_t d = 0; d < dim; ++d) {
        if (std::isinf(vmin[d]) || std::isnan(vmin[d]) ||
            std::isinf(vmax[d]) || std::isnan(vmax[d])) {
            throw std::runtime_error("Blocks::setup_bounds error: bad limits");
        }

        double dv = std::abs(vmax[d] - vmin[d]);
        m_vmin[d] = vmin[d] - EC * dv;
        m_vmax[d] = vmax[d] + EC * dv;
    }

    init_coords();
    init_widths();

    m_ready = true;
}

void Blocks::update_bounds(const Vector3d &vmin, const Vector3d &vmax) {
    readiness();

    int dim = m_map.type().size();
    for (size_t d = 0; d < dim; ++d) {
        if (std::isinf(vmin[d]) || std::isnan(vmin[d]) ||
            std::isinf(vmax[d]) || std::isnan(vmax[d])) {
            throw std::runtime_error("Blocks::update_bounds error: bad limits");
        }

        double dv = std::abs(vmax[d] - vmin[d]);
        m_vmin[d] = vmin[d] - EC * dv;
        m_vmax[d] = vmax[d] + EC * dv;
    }

    if (dim > 0) {
        m_x_coords.front() = m_vmin.x;
        m_x_coords.back()  = m_vmax.x;

        if (dim > 1) {
            for (int i = 0; i < m_nx; ++i) {
                m_y_coords[i].front() = m_vmin.y;
                m_y_coords[i].back()  = m_vmax.y;

                if (dim > 2) {
                    for (int j = 0; j < m_ny[i]; ++j) {
                        m_z_coords[i][j].front() = m_vmin.z;
                        m_z_coords[i][j].back()  = m_vmax.z;
                    }
                }
            }
        }
    }

    init_widths();
}

int Blocks::rank(const Vector3d& v) const {
    readiness();

    Vector3d vec = m_map.mapping(v);

    int i = 0;
    for (; i < m_nx - 1; ++i) {
        if (vec.x < m_x_coords[i + 1]) {
            break;
        }
    }

    int j = 0;
    for (; j < m_ny[i] - 1; ++j) {
        if (vec.y < m_y_coords[i][j + 1]) {
            break;
        }
    }

    int k = 0;
    for (; k < m_nz[i][j] + 1; ++k) {
        if (vec.z < m_z_coords[i][j][k + 1]) {
            break;
        }
    }

    return m_ranks[i][j][k];
}

size_t Blocks::size() const {
    return m_size;
}

size_t Blocks::size(int i) const {
    size_t res = 0;
    for (int j = 0; j < m_ny[i]; ++j) {
        res += m_nz[i][j];
    }
    return res;
}

size_t Blocks::size(int i, int j) const {
    return m_nz[i][j];
}

Vector3d Blocks::center(int rank) const {
    int i = m_i[rank];
    int j = m_j[rank];
    int k = m_k[rank];

    auto dim = m_map.dimension();

    Vector3d c = {0.0, 0.0, 0.0};
    c.x = m_nx > 1 ? (m_x_coords[i] + m_x_coords[j]) / 2.0 : 0.0;
    if (dim > 1) {
        c.y = m_ny[i] > 1 ? (m_y_coords[i][j] + m_y_coords[i][j + 1]) / 2.0 : 0.0;

        if (dim > 2) {
            c.z = m_nz[i][j] > 1 ? (m_z_coords[i][j][k] + m_z_coords[i][j][k + 1]) / 2.0 : 0.0;
        }
    }

    return m_map.inverse(c);
}

std::vector<std::vector<Vector3d>> Blocks::lines() const {
    readiness();

    std::vector<std::vector<Vector3d>> L(m_size);

    auto add_line = [this](
            std::vector<Vector3d>& line,
            Vector3d v1, Vector3d v2) {

        const int M = 50;
        for (int k = 0; k < M; ++k) {
            Vector3d v;
            if (std::isinf(v1.x) || std::isinf(v2.x)) {
                double xi1 = std::atan(v1.x);
                double xi2 = std::atan(v2.x);
                v.x = std::tan(xi1 + ((xi2 - xi1) * k) / (M - 1));
            }
            else {
                v.x = v1.x + ((v2.x - v1.x) * k) / M;
            }

            if (std::isinf(v1.y) || std::isinf(v2.y)) {
                double xi1 = std::atan(v1.y);
                double xi2 = std::atan(v2.y);
                v.y = std::tan(xi1 + ((xi2 - xi1) * k) / (M - 1));
            }
            else {
                v.y = v1.y + ((v2.y - v1.y) * k) / M;
            }
            line.push_back(m_map.inverse(v));
        }
    };

    for (size_t r = 0; r < m_size; ++r) {
        int i = m_i[r];
        int j = m_j[r];

        Vector3d p1 = {m_x_coords[i], m_y_coords[i][j], 0.0};
        Vector3d p2 = {m_x_coords[i + 1], m_y_coords[i][j], 0.0};
        Vector3d p3 = {m_x_coords[i + 1], m_y_coords[i][j + 1], 0.0};
        Vector3d p4 = {m_x_coords[i], m_y_coords[i][j + 1], 0.0};

        add_line(L[r], p1, p2);
        add_line(L[r], p2, p3);
        add_line(L[r], p3, p4);
        add_line(L[r], p4, p1);
        L[r].push_back(L[r][0]);
    }
    return L;
}

void Blocks::set_loads(const std::vector<double>& loads) {
    if (loads.size() != m_size) {
        throw std::runtime_error("Blocks::set_loads error #1");
    }
    int N = static_cast<int>(m_size);

    // Нагрузка блоков
    for (int r = 0; r < N; ++r) {
        int i = m_i[r];
        int j = m_j[r];
        int k = m_k[r];

        m_z_loads[i][j][k] = loads[r];
    }

    for (int i = 0; i < m_nx; ++i) {
        int count = 0;
        double load_x = 0.0;
        for (int j = 0; j < m_ny[i]; ++j) {
            double load_y = 0.0;
            for (int k = 0; k < m_nz[i][j]; ++k) {
                load_y += m_z_loads[i][j][k];
            }
            load_x += load_y;
            count += m_nz[i][j];
            m_y_loads[i][j] = load_y / m_nz[i][j];
        }
        m_x_loads[i] = load_x / count;
    }
}

bool Blocks::is_near(int m, int neib_rank, const Vector3d &v, double R) {
    Vector3d c = m_map(v);

    int i = m_i[neib_rank];
    int j = m_j[neib_rank];
    int k = m_k[neib_rank];

    if (m_x_coords[i] - R <= c.x && c.x < m_x_coords[i + 1] + R) {
        if (m_y_coords[i][j] - R <= c.y && c.y < m_y_coords[i][j + 1] + R) {
            if (m_z_coords[i][j][k] - R <= c.z && c.z < m_z_coords[i][j][k + 1] + R) {
                return true;
            }
        }
    }
    return false;
}

void Blocks::move_simple(
        double alpha,
        const std::vector<double>& loads,
        const std::vector<double>& widths,
        std::vector<barrier>& coords) {

    int n_blocks = int(loads.size());

    if (n_blocks < 2) {
        return;
    }

    for (int i = 1; i < n_blocks; ++i) {
        double L1 = loads[i - 1];
        double L2 = loads[i];
        double dx1 = widths[i - 1];
        double dx2 = widths[i];

        double dx = std::min(dx1, dx2);
        double I = L1 != L2 ? (L2 - L1) / (L2 + L1) : 0.0;
        coords[i] += alpha * dx * I;
    }
}

// Метод прогонки, решение для "почти" трехдиагональной СЛАУ вида
// A_i u_{i-1} + B_i u_i + C_i u_{i+1} = F_i, i = 1..N-1
// B_0 u_0 + C_0 u_1 + A_0 u_2 = F_0
// C_N u_{N-2} + A_N u_{N-1} + B_N u_N = F_N
std::vector<double> progon(std::vector<double> A, std::vector<double> B,std::vector<double>  C, std::vector<double> F) {
    auto N = A.size();

    std::vector<double> a(N, 0.0);
    std::vector<double> b(N, 0.0);

    a[1] = -C[0] / B[0];
    b[1] = F[0] / B[0];
    for (size_t i = 1; i < N - 1; ++i) {
        double Q = (A[i] * a[i] + B[i]);
        a[i + 1] = -C[i] / Q;
        b[i + 1] = (F[i] - A[i] * b[i]) / Q;
    }

    std::vector<double> u(N, 0.0);
    u[N - 1] = (F[N - 1] - A[N - 1] * b[N - 1]) / (B[N - 1] + A[N - 1] * a[N - 1]);

    for (size_t i = N - 1; i > 0; --i) {
        u[i - 1] = a[i] * u[i] + b[i];
    }

    return u;
}

inline double fit(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

inline double imb_func(double Imax) {
    // Значение дисбаланса, при котором возвращается 1/2
    // Приемлимое значение дисбаланса
    const double I0 = 1.0e-1;
    return Imax / (Imax + I0);
}

void Blocks::move_newton(
        double alpha,
        const std::vector<double>& loads,
        const std::vector<double>& widths,
        std::vector<barrier>& coords) {

    zephyr::network::mpi::Network& net = zephyr::network::mpi::main();

    int n_blocks = int(loads.size());

    if (n_blocks < 2) {
        return;
    }

    // Центры ячеек
    std::vector<double> xc(n_blocks);
    for (int i = 0; i < n_blocks; ++i) {
        xc[i] = 0.5 * (coords[i] + coords[i+1]);
    }

    // Плотность нагрузки
    std::vector<double> dens_c(n_blocks);
    for (int i = 0; i < n_blocks; ++i) {
        dens_c[i] = std::max(0.0, loads[i] / (coords[i + 1] - coords[i]));
    }

    // Плотность на границах блоков
    std::vector<double> dens(n_blocks + 1);
    for (int i = 1; i < n_blocks; ++i) {
        double nu = (xc[i] - coords[i]) / (xc[i] - xc[i - 1]);
        dens[i] = nu * dens_c[i - 1] + (1.0 - nu) * dens_c[i];
    }
    dens[0] = 0.0;
    dens[n_blocks] = 0.0;

    std::vector<double> A(n_blocks - 1, 0.0);
    std::vector<double> B(n_blocks - 1, 0.0);
    std::vector<double> C(n_blocks - 1, 0.0);
    std::vector<double> F(n_blocks - 1, 0.0);

    for (int i = 1; i < n_blocks; ++i) {
        double L1 = loads[i - 1];
        double L2 = loads[i];

        double Imb = L1 != L2 ? std::abs((L2 - L1) / (L2 + L1)) : 0.0;

        A[i - 1] = -dens[i - 1];
        B[i - 1] = 2.0 * dens[i];
        C[i - 1] = -dens[i + 1];

        F[i - 1] = (L2 - L1) * imb_func(Imb);
    }

    double err = 0.0;
    auto dx = progon(A, B, C, F);

    for (int i = 1; i < n_blocks; ++i) {
        dx[i - 1] = fit(dx[i - 1], alpha * (coords[i - 1] - coords[i]), alpha * (coords[i + 1] - coords[i]));

        //if (coords[i].dx * dx[i - 1] < 0.0) {
        //    dx[i - 1] *= std::min(std::abs(0.05 * coords[i].dx / dx[i - 1]), 1.0);
        //}
        double L1 = loads[i - 1];
        double L2 = loads[i];

        double I = L1 != L2 ? (L2 - L1) / (L2 + L1) : 0.0;

        coords[i].dx = dx[i - 1]; // * imb_func(std::abs(I));
        coords[i] += dx[i - 1];
    }
}

void Blocks::balancing_simple(const std::vector<double>& loads, double alpha) {
    readiness();

    if (loads.size() != m_size) {
        throw std::runtime_error("loads.size != decomposition.size");
    }
    set_loads(loads);

    alpha = std::max(0.0, std::min(alpha, 0.45));

    move_simple(alpha, m_x_loads, m_x_widths, m_x_coords);

    for (int i = 0; i < m_nx; ++i) {
        move_simple(alpha, m_y_loads[i], m_y_widths[i], m_y_coords[i]);
    }

    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < m_ny[i]; ++j) {
            move_simple(alpha, m_z_loads[i][j], m_z_widths[i][j], m_z_coords[i][j]);
        }
    }

    init_widths();
}

void Blocks::balancing_newton(const std::vector<double>& loads, double alpha) {
    readiness();

    if (loads.size() != m_size) {
        throw std::runtime_error("loads.size != decomposition.size");
    }
    set_loads(loads);

    alpha = std::max(0.0, std::min(alpha, 0.45));

    move_newton(alpha, m_x_loads, m_x_widths, m_x_coords);

    for (int i = 0; i < m_nx; ++i) {
        move_newton(alpha, m_y_loads[i], m_y_widths[i], m_y_coords[i]);
    }

    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < m_ny[i]; ++j) {
            move_newton(alpha, m_z_loads[i][j], m_z_widths[i][j], m_z_coords[i][j]);
        }
    }

    init_widths();
}

void Blocks::init_sizes_1D() {
    m_nx = m_size;
    m_ny.resize(m_nx);
    m_nz.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ny[i] = 1;
        m_nz[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_nz[i][j] = 1;
        }
    }
}

void Blocks::init_sizes_2D(int n_proc_1) {
    if (n_proc_1 < 1) {
        std::string message = "Для двумерной декомпозиции необходимо указать n_proc_1";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (n_proc_1 >= m_size) {
        n_proc_1 = m_size;
    }
    int n_proc_2 = m_size / n_proc_1;

    int big_count = m_size - n_proc_1 * n_proc_2;

    m_nx = n_proc_1;
    m_ny.resize(m_nx);
    m_nz.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ny[i] = n_proc_2;
    }

    int counter = 0;
    int beg = (m_nx - big_count) / 2;
    while (counter < big_count) {
        m_ny[beg + counter] += 1;
        counter += 1;
    }

    m_nz.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_nz[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_nz[i][j] = 1;
        }
    }
}

void Blocks::init_sizes_3D(int n_proc_1, int n_proc_2) {
    if (n_proc_1 < 1 || n_proc_2 < 1) {
        std::string message = "Для трехмерной декомпозиции необходимо указать n_proc_1 и n_proc_2";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    // Уменьшим значения n_proc_1, n_proc_2
    if (n_proc_1 * n_proc_2 > m_size) {
        int counter = 0;
        while (n_proc_1 * n_proc_2 > m_size) {
            if (n_proc_1 > 1) {
                n_proc_1 -= (counter + 1) % 2;
            }
            if (n_proc_2 > 1) {
                n_proc_2 -= counter % 2;
            }
            ++counter;
        }
    }

    int n_proc_12 = n_proc_1 * n_proc_2;

    int n_proc_3 = m_size / n_proc_12;

    int big_count = m_size - n_proc_12 * n_proc_3;

    m_nx = n_proc_1;
    m_ny.resize(m_nx);
    m_nz.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ny[i] = n_proc_2;
        m_nz[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_nz[i][j] = n_proc_3;
        }
    }

    struct Index {
        int i, j;

        Index(int _i, int _j) : i(_i), j(_j) {}

        double abs(const int &n_proc_1, const int &n_proc_2) const {
            double x = i - 0.5 * (n_proc_1 - 1.0);
            double y = j - 0.5 * (n_proc_2 - 1.0);
            return x * x + y * y;
        }
    };

    std::vector<Index> indices;
    indices.reserve(n_proc_12);
    for (int i = 0; i < n_proc_1; ++i) {
        for (int j = 0; j < n_proc_2; ++j) {
            indices.emplace_back(Index(i, j));
        }
    }

    // Сортируем двумерные индексы в квадрате n_proc_1 x n_proc_2
    // по убыванию от центра области ( 0.5 n_proc_1, 0.5 n_proc_2 )
    std::sort(indices.begin(), indices.end(),
              [n_proc_1, n_proc_2](const Index &idx1, const Index &idx2) {
                  return idx1.abs(n_proc_1, n_proc_2) < idx2.abs(n_proc_1, n_proc_2);
              });

    for (int counter = 0; counter < big_count; ++counter) {
        m_nz[indices[counter].i][indices[counter].j] += 1;
    }
}

void Blocks::init_limits() {
    m_vmin = m_map.limits().min;
    m_vmax = m_map.limits().max;
}

void Blocks::init_indices() {
    int counter = 0;
    m_ranks.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ranks[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_ranks[i][j].resize(m_nz[i][j]);
            for (int k = 0; k < m_nz[i][j]; ++k) {
                m_i.push_back(i);
                m_j.push_back(j);
                m_k.push_back(k);
                m_ranks[i][j][k] = counter;
                ++counter;
            }
        }
    }
    if (counter > m_size) {
        throw std::runtime_error("Blocks::init_indices error");
    }
}

void Blocks::init_coords() {
    m_x_coords.resize(m_nx + 1);
    if (m_nx == 1) {
        m_x_coords[0] = m_map.limits().min.x;
        m_x_coords[1] = m_map.limits().max.x;
    }
    else {
        // Тут небольшая поправочка на EC, поскольку используется
        // расширенная область (больше domain на долю EC)
        double dx = (m_vmax.x - m_vmin.x) * EC / (1 + 2 * EC);
        double xi = (m_vmax.x - m_vmin.x - 2 * dx) / m_size;

        m_x_coords[0] = m_vmin.x + dx;
        for (int i = 1; i < m_nx; ++i) {
            m_x_coords[i] = m_x_coords[i - 1] + xi * size(i - 1);
        }
        m_x_coords[0] = m_vmin.x;
        m_x_coords[m_nx] = m_vmax.x;
    }

    m_y_coords.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_y_coords[i].resize(m_ny[i] + 1);
        if (m_ny[i] == 1) {
            m_y_coords[i][0] = m_map.limits().min.y;
            m_y_coords[i][1] = m_map.limits().max.y;
        }
        else {
            // Тут небольшая поправочка на EC, поскольку используется
            // расширенная область (больше domain на долю EC)
            double dy = (m_vmax.y - m_vmin.y) * EC / (1 + 2 * EC);
            double xi = (m_vmax.y - m_vmin.y - 2 * dy) / size(i);

            m_y_coords[i][0] = m_vmin.y + dy;
            for (int j = 1; j < m_ny[i]; ++j) {
                m_y_coords[i][j] = m_y_coords[i][j - 1] + xi * size(i, j - 1);
            }
            m_y_coords[i][0] = m_vmin.y;
            m_y_coords[i][m_ny[i]] = m_vmax.y;
        }
    }

    m_z_coords.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_z_coords[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_z_coords[i][j].resize(m_nz[i][j] + 1);
            if (m_nz[i][j] == 1) {
                m_z_coords[i][j][0] = m_map.limits().min.z;
                m_z_coords[i][j][1] = m_map.limits().max.z;
            }
            else {
                // Тут небольшая поправочка на EC, поскольку используется
                // расширенная область (больше domain на долю EC)
                double dz = (m_vmax.z - m_vmin.z) * EC / (1 + 2 * EC);

                m_z_coords[i][j][0] = m_vmin.z;
                for (int k = 1; k < m_nz[i][j]; ++k) {
                    m_z_coords[i][j][k] = m_vmin.z + dz + ((m_vmax.z - m_vmin.z - 2 * dz) * k) / m_nz[i][j];
                }
                m_z_coords[i][j][m_nz[i][j]] = m_vmax.z;
            }
        }
    }
}

void Blocks::init_widths() {
    m_x_widths.resize(m_nx);
    m_y_widths.resize(m_nx);
    m_z_widths.resize(m_nx);

    for (int i = 0; i < m_nx; ++i) {
        m_x_widths[i] = m_x_coords[i + 1] - m_x_coords[i];
        if (m_x_coords.size() <= i + 1) {
            std::cerr << m_x_coords.size() << " " << i << " " << m_nx << "\n";
            throw std::runtime_error("Wtf 1");
        }

        m_y_widths[i].resize(m_ny[i]);
        m_z_widths[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_y_widths[i][j] = m_y_coords[i][j + 1] - m_y_coords[i][j];
            if (m_y_coords[i].size() <= j + 1) {
                std::cerr << m_y_coords[i].size() << " " << j << " " << m_ny[i] << "\n";
                throw std::runtime_error("Wtf 2");
            }

            m_z_widths[i][j].resize(m_nz[i][j]);
            for (int k = 0; k < m_nz[i][j]; ++k) {
                m_z_widths[i][j][k] = m_z_coords[i][j][k + 1] - m_z_coords[i][j][k];
                if (m_z_coords[i][j].size() <= k + 1) {
                    throw std::runtime_error("Wtf 3");
                }
            }
        }
    }
}

void Blocks::init_loads() {
    m_x_loads.resize(m_nx);
    m_y_loads.resize(m_nx);
    m_z_loads.resize(m_nx);

    for (int i = 0; i < m_nx; ++i) {
        m_x_loads[i] = 0.0;

        m_y_loads[i].resize(m_ny[i]);
        m_z_loads[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_y_loads[i][j] = 0.0;

            m_z_loads[i][j].resize(m_nz[i][j]);
            for (int k = 0; k < m_nz[i][j]; ++k) {
                m_z_loads[i][j][k] = 0.0;
            }
        }
    }
}

void Blocks::info() const {
    readiness();

    std::cout << "Информация по блочной декомпозиции.\n\n";
    std::cout << "Индексы рангов:\n";
    for (size_t r = 0; r < m_i.size(); ++r) {
        std::cout << "    Ранг " << r << "    (" <<
                  m_i[r] << ", " << m_j[r] << ", " << m_k[r] << ")\n";
    }
    std::cout << "Блоков по X:   " << m_nx << "\n";
    std::cout << "Ширина по X:   ";
    for (size_t i = 0; i < m_nx; ++i) {
        std::cout << m_x_widths[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Нагрузка по X: ";
    for (size_t i = 0; i < m_nx; ++i) {
        std::cout << m_x_loads[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Разделители:   ";
    for (size_t i = 0; i <= m_nx; ++i) {
        std::cout << m_x_coords[i] << " ";
    }
    std::cout << "\n";
    for (size_t i = 0; i < m_nx; ++i) {
        std::cout << "    X блок " << i << "\n";
        std::cout << "        Блоков по Y:   " << m_ny[i] << "\n";
        std::cout << "        Ширина по Y:   ";
        for (size_t j = 0; j < m_ny[i]; ++j) {
            std::cout << m_y_widths[i][j] << " ";
        }
        std::cout << "\n";
        std::cout << "        Нагрузка по Y: ";
        for (size_t j = 0; j < m_ny[i]; ++j) {
            std::cout << m_y_loads[i][j] << " ";
        }
        std::cout << "\n";
        std::cout << "        Разделители:   ";
        for (size_t j = 0; j <= m_ny[i]; ++j) {
            std::cout << m_y_coords[i][j] << " ";
        }
        std::cout << "\n";
        for (size_t j = 0; j < m_ny[i]; ++j) {
            std::cout << "            Y блок " << j << "\n";
            std::cout << "                Блоков по Z:   " << m_nz[i][j] << "\n";
            std::cout << "                Ширина по Z:   ";
            for (size_t k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_z_widths[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Нагрузка по Z: ";
            for (size_t k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_z_loads[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Разделители:   ";
            for (size_t k = 0; k <= m_nz[i][j]; ++k) {
                std::cout << m_z_coords[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Ранги:         ";
            for (size_t k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_ranks[i][j][k] << " ";
            }
            std::cout << "\n";
        }
    }
}

} // decomposition
} // network
} // zephyr