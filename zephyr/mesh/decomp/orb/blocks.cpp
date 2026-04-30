#include <iostream>
#include <algorithm>
#include <numeric>
#include <span>
#include <list>
#include <set>

#include <zephyr/mesh/decomp/orb/blocks.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh::decomp {
static constexpr double inf = std::numeric_limits<double>::infinity();

// Декомпозиция прямоугольника на size блоков,
// по оси x задано nx блоков
inline std::vector<int> decomp(int size, int n_x) {
    int N_YZ = size / n_x;
    int big_count = size - n_x * N_YZ;

    std::vector<int> n_y(n_x);
    for (int i = 0; i < n_x; ++i) {
        n_y[i] = N_YZ;
    }

    int counter = 0;
    int beg = (n_x - big_count) / 2;
    while (counter < big_count) {
        n_y[beg + counter] += 1;
        counter += 1;
    }

    int check_sum = std::accumulate(n_y.begin(), n_y.end(), 0);
    if (check_sum != size) {
        throw std::runtime_error("Something wrong in Blocks constructor");
    }

    return n_y;
}

// Декомпозиция прямоугольника [Lx, Ly] на size блоков
inline std::vector<int> decomp(double Lx, double Ly, int size) {
    double H = std::sqrt((Lx * Ly) / size);
    int nx = static_cast<int>(std::round(Lx / H ));
    nx = std::min(std::max(1, nx), size);
    return decomp(size, nx);
}

Blocks::Blocks() {
    m_ready = true;

    m_map = Transform();

    m_size = 1;
    m_nx = 1;
    m_ny = {1};
    m_nz = {{1}};

    m_ids = {index_t{0, 0, 0}};

    m_ranks = {{{0}}};

    m_x_coords = {  m_map.x_min(), m_map.x_max()};
    m_y_coords = {{ m_map.y_min(), m_map.y_max()}};
    m_z_coords = {{{m_map.z_min(), m_map.z_max()}}};

    m_x_widths = {  m_map.x_width()};
    m_y_widths = {{ m_map.y_width()}};
    m_z_widths = {{{m_map.z_width()}}};

    m_x_loads = {0.0};
    m_y_loads = {{0.0}};
    m_z_loads = {{{0.0}}};
}

Blocks::Blocks(const Box& domain, const std::string& type, int size)
    : m_size(size), m_map(domain, type) {

    if (type.size() == 1) {
        init_sizes_1D(m_size);
    }
    else if (type.size() == 2) {
        Vector3d L = m_map.box().sizes();
        std::vector<int> n_y = decomp(L.x(), L.y(), size);
        init_sizes_2D(n_y);
    }
    else if (type.size() == 3) {
        Vector3d L = m_map.box().sizes();
        double H = std::cbrt((L.x() * L.y() * L.z()) / m_size);
        int nx = static_cast<int>(std::round(L.x() / H ));
        nx = std::min(std::max(1, nx), m_size);

        std::vector<int> nYZ = decomp(size, nx);
        std::vector<std::vector<int>> nz(nx);
        for (int i = 0; i < nx; ++i) {
            nz[i] = decomp(L.y(), L.z(), nYZ[i]);
        }
        int check_sum = 0;
        for (auto& r: nz) {
            check_sum += std::accumulate(r.begin(), r.end(), 0);
        }
        if (check_sum != size) {
            throw std::runtime_error("Something wrong in Blocks constructor");
        }
        init_sizes_3D(nz);
    }
    else {
        throw std::runtime_error("Странная декомпозиция с размерностью > 3");
    }

    init_indices();
    init_coords();
    init_widths();
    init_loads();
    m_ready = true;
}

Blocks::Blocks(const Box& domain, const std::string& type, int size, int nx)
    : m_size(size), m_map(domain, type) {

    if (type.size() == 1) {
        if (size < 0) {
            m_size = nx;
        }
        init_sizes_1D(m_size);
    }
    else if (type.size() == 2) {
        std::vector<int> n_y = decomp(size, nx);
        init_sizes_2D(n_y);
    }
    else if (type.size() == 3) {
        std::vector<int> nYZ = decomp(size, nx);
        std::vector<std::vector<int>> nz(nx);
        for (int i = 0; i < nx; ++i) {
            Vector3d L = m_map.box().sizes();
            nz[i] = decomp(L.y(), L.z(), nYZ[i]);
        }
        int check_sum = 0;
        for (auto& r: nz) {
            check_sum += std::accumulate(r.begin(), r.end(), 0);
        }
        if (check_sum != size) {
            throw std::runtime_error("Something wrong in Blocks constructor");
        }
        init_sizes_3D(nz);
    }
    else {
        throw std::runtime_error("Странная декомпозиция с размерностью > 3");
    }

    init_indices();
    init_coords();
    init_widths();
    init_loads();
    m_ready = true;
}

Blocks::Blocks(const Box& domain, const std::string& type, int size,
               const std::vector<int>& ny)
    : m_size(size), m_map(domain, type) {

    if (type.size() == 1) {
        init_sizes_1D(m_size);
    }
    else if (type.size() == 2) {
        if (size < 0) {
            m_size = std::accumulate(ny.begin(), ny.end(), 0);
        }
        init_sizes_2D(ny);
    }
    else if (type.size() == 3) {
        Vector3d L = m_map.box().sizes();
        int nx = ny.size();
        std::vector<std::vector<int>> nz(nx);
        for (int i = 0; i < nx; ++i) {
            nz[i] = decomp(L.y(), L.z(), ny[i]);
        }

        int check_size = 0;
        for (auto& r: nz) {
            check_size += std::accumulate(r.begin(), r.end(), 0);
        }
        if (size < 0) {
            m_size = check_size;
        }
        else if (check_size != size) {
            throw std::runtime_error("Something wrong in Blocks constructor");
        }
        init_sizes_3D(nz);
    }
    else {
        throw std::runtime_error("Странная декомпозиция с размерностью > 3");
    }

    init_indices();
    init_coords();
    init_widths();
    init_loads();
    m_ready = true;
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

void Blocks::update_bounds(const Box& domain) {
    readiness();

    m_map.update_limits(domain);

    int dim = m_map.dim();
    if (dim > 0) {
        m_x_coords.front() = m_map.x_min();
        m_x_coords.back()  = m_map.x_max();

        if (dim > 1) {
            for (int i = 0; i < m_nx; ++i) {
                m_y_coords[i].front() = m_map.y_min();
                m_y_coords[i].back()  = m_map.y_max();

                if (dim > 2) {
                    for (int j = 0; j < m_ny[i]; ++j) {
                        m_z_coords[i][j].front() = m_map.z_min();
                        m_z_coords[i][j].back()  = m_map.z_max();
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
        if (vec.x() < m_x_coords[i + 1]) {
            break;
        }
    }

    int j = 0;
    for (; j < m_ny[i] - 1; ++j) {
        if (vec.y() < m_y_coords[i][j + 1]) {
            break;
        }
    }

    int k = 0;
    for (; k < m_nz[i][j] + 1; ++k) {
        if (vec.z() < m_z_coords[i][j][k + 1]) {
            break;
        }
    }

    return m_ranks[i][j][k];
}

int Blocks::size() const {
    return m_size;
}

int Blocks::size(int i) const {
    int res = 0;
    for (int j = 0; j < m_ny[i]; ++j) {
        res += m_nz[i][j];
    }
    return res;
}

int Blocks::size(int i, int j) const {
    return m_nz[i][j];
}

Vector3d Blocks::center(int rank) const {
    auto [i, j, k] = m_ids[rank];

    const auto dim = m_map.dim();

    Vector3d c = {0.0, 0.0, 0.0};
    c.x() = m_nx > 1 ? (m_x_coords[i] + m_x_coords[j]) / 2.0 : 0.0;
    if (dim > 1) {
        c.y() = m_ny[i] > 1 ? (m_y_coords[i][j] + m_y_coords[i][j + 1]) / 2.0 : 0.0;

        if (dim > 2) {
            c.z() = m_nz[i][j] > 1 ? (m_z_coords[i][j][k] + m_z_coords[i][j][k + 1]) / 2.0 : 0.0;
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
            if (std::isinf(v1.x()) || std::isinf(v2.x())) {
                double xi1 = std::atan(v1.x());
                double xi2 = std::atan(v2.x());
                v.x() = std::tan(xi1 + ((xi2 - xi1) * k) / (M - 1));
            }
            else {
                v.x() = v1.x() + ((v2.x() - v1.x()) * k) / M;
            }

            if (std::isinf(v1.y()) || std::isinf(v2.y())) {
                double xi1 = std::atan(v1.y());
                double xi2 = std::atan(v2.y());
                v.y() = std::tan(xi1 + ((xi2 - xi1) * k) / (M - 1));
            }
            else {
                v.y() = v1.y() + ((v2.y() - v1.y()) * k) / M;
            }
            line.push_back(m_map.inverse(v));
        }
    };

    for (int r = 0; r < m_size; ++r) {
        auto [i, j, k] = m_ids[r];

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

    // Нагрузка блоков
    for (int r = 0; r < m_size; ++r) {
        auto [i, j, k] = m_ids[r];
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
std::vector<double> progon(std::vector<double> A, std::vector<double> B, std::vector<double>  C, std::vector<double> F) {
    auto N = A.size();

    std::vector<double> a(N, 0.0);
    std::vector<double> b(N, 0.0);

    a[1] = -C[0] / B[0];
    b[1] = F[0] / B[0];
    for (int i = 1; i < N - 1; ++i) {
        double Q = (A[i] * a[i] + B[i]);
        a[i + 1] = -C[i] / Q;
        b[i + 1] = (F[i] - A[i] * b[i]) / Q;
    }

    std::vector<double> u(N, 0.0);
    u[N - 1] = (F[N - 1] - A[N - 1] * b[N - 1]) / (B[N - 1] + A[N - 1] * a[N - 1]);

    for (int i = N - 1; i > 0; --i) {
        u[i - 1] = a[i] * u[i] + b[i];
    }

    return u;
}

inline double fit(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

inline double imb_func(double Imax) {
    // Значение дисбаланса, при котором возвращается 1/2
    // Приемлемое значение дисбаланса
    constexpr double I0 = 0.1;
    return Imax / (Imax + I0);
}

void Blocks::move_newton(
        double alpha,
        const std::vector<double>& loads,
        const std::vector<double>& widths,
        std::vector<barrier>& coords) {

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

void Blocks::init_sizes_1D(int nx) {
    if (nx != m_size) {
        throw std::runtime_error("Blocks::init_sizes_1D error: Wrong 1D decomposition");
    }

    m_nx = nx;

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

void Blocks::init_sizes_2D(const std::vector<int>& ny) {
    int checkout_sum = std::accumulate(ny.begin(), ny.end(), 0);
    if (checkout_sum != m_size) {
        throw std::runtime_error("Blocks::init_sizes_2D error: Wrong 2D decomposition");
    }

    m_nx = ny.size();

    m_ny.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ny[i] = ny[i];
    }

    m_nz.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_nz[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_nz[i][j] = 1;
        }
    }
}

void Blocks::init_sizes_3D(const std::vector<std::vector<int>>& nz) {
    int size = 0;
    for(auto& r: nz) {
        size += std::accumulate(r.begin(), r.end(), 0);
    }
    if (size != m_size) {
        throw std::runtime_error("Blocks::init_sizes_3D error: Wrong 3D decomposition");
    }

    m_nz = nz;
    m_nx = nz.size();
    m_ny.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ny[i] = nz[i].size();
    }
}

void Blocks::init_indices() {
    int counter = 0;
    m_ranks.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_ranks[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_ranks[i][j].resize(m_nz[i][j]);
            for (int k = 0; k < m_nz[i][j]; ++k) {
                m_ids.emplace_back(i, j, k);
                m_ranks[i][j][k] = counter;
                ++counter;
            }
        }
    }
    if (counter > m_size) {
        throw std::runtime_error("Blocks::init_indices error");
    }
}

constexpr double EC = 0.01;

void Blocks::init_coords() {
    m_x_coords.resize(m_nx + 1);
    if (m_nx == 1) {
        m_x_coords[0] = m_map.x_min();
        m_x_coords[1] = m_map.x_max();
    }
    else {
        // Тут небольшая поправочка на EC, поскольку используется
        // расширенная область (больше domain на долю EC)
        double dx = m_map.x_width() * EC / (1 + 2 * EC);
        double xi = (m_map.x_width() - 2 * dx) / m_size;

        m_x_coords[0] = m_map.x_min() + dx;
        for (int i = 1; i < m_nx; ++i) {
            m_x_coords[i] = m_x_coords[i - 1] + xi * size(i - 1);
        }
        m_x_coords[0] = m_map.x_min();
        m_x_coords[m_nx] = m_map.x_max();
    }

    m_y_coords.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_y_coords[i].resize(m_ny[i] + 1);
        if (m_ny[i] == 1) {
            m_y_coords[i][0] = m_map.y_min();
            m_y_coords[i][1] = m_map.y_max();
        }
        else {
            // Тут небольшая поправочка на EC, поскольку используется
            // расширенная область (больше domain на долю EC)
            double dy = m_map.y_width() * EC / (1 + 2 * EC);
            double xi = (m_map.y_width() - 2 * dy) / size(i);

            m_y_coords[i][0] = m_map.y_min() + dy;
            for (int j = 1; j < m_ny[i]; ++j) {
                m_y_coords[i][j] = m_y_coords[i][j - 1] + xi * size(i, j - 1);
            }
            m_y_coords[i][0] = m_map.y_min();
            m_y_coords[i][m_ny[i]] = m_map.y_max();
        }
    }

    m_z_coords.resize(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        m_z_coords[i].resize(m_ny[i]);
        for (int j = 0; j < m_ny[i]; ++j) {
            m_z_coords[i][j].resize(m_nz[i][j] + 1);
            if (m_nz[i][j] == 1) {
                m_z_coords[i][j][0] = m_map.z_min();
                m_z_coords[i][j][1] = m_map.z_max();
            }
            else {
                // Тут небольшая поправочка на EC, поскольку используется
                // расширенная область (больше domain на долю EC)
                double dz = m_map.z_width() * EC / (1 + 2 * EC);

                m_z_coords[i][j][0] = m_map.z_min();
                for (int k = 1; k < m_nz[i][j]; ++k) {
                    m_z_coords[i][j][k] = m_map.z_min() + dz + ((m_map.z_width() - 2 * dz) * k) / m_nz[i][j];
                }
                m_z_coords[i][j][m_nz[i][j]] = m_map.z_max();
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
    for (int r = 0; r < m_ids.size(); ++r) {
        auto [i, j, k] = m_ids[r];
        std::cout << "    Ранг " << r << "    (" << i << ", " << j << ", " << k << ")\n";
    }
    std::cout << "Блоков по X:   " << m_nx << "\n";
    std::cout << "Ширина по X:   ";
    for (int i = 0; i < m_nx; ++i) {
        std::cout << m_x_widths[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Нагрузка по X: ";
    for (int i = 0; i < m_nx; ++i) {
        std::cout << m_x_loads[i] << " ";
    }
    std::cout << "\n";
    std::cout << "Разделители:   ";
    for (int i = 0; i <= m_nx; ++i) {
        std::cout << m_x_coords[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < m_nx; ++i) {
        std::cout << "    X блок " << i << "\n";
        std::cout << "        Блоков по Y:   " << m_ny[i] << "\n";
        std::cout << "        Ширина по Y:   ";
        for (int j = 0; j < m_ny[i]; ++j) {
            std::cout << m_y_widths[i][j] << " ";
        }
        std::cout << "\n";
        std::cout << "        Нагрузка по Y: ";
        for (int j = 0; j < m_ny[i]; ++j) {
            std::cout << m_y_loads[i][j] << " ";
        }
        std::cout << "\n";
        std::cout << "        Разделители:   ";
        for (int j = 0; j <= m_ny[i]; ++j) {
            std::cout << m_y_coords[i][j] << " ";
        }
        std::cout << "\n";
        for (int j = 0; j < m_ny[i]; ++j) {
            std::cout << "            Y блок " << j << "\n";
            std::cout << "                Блоков по Z:   " << m_nz[i][j] << "\n";
            std::cout << "                Ширина по Z:   ";
            for (int k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_z_widths[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Нагрузка по Z: ";
            for (int k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_z_loads[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Разделители:   ";
            for (int k = 0; k <= m_nz[i][j]; ++k) {
                std::cout << m_z_coords[i][j][k] << " ";
            }
            std::cout << "\n";
            std::cout << "                Ранги:         ";
            for (int k = 0; k < m_nz[i][j]; ++k) {
                std::cout << m_ranks[i][j][k] << " ";
            }
            std::cout << "\n";
        }
    }
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& values) {
    for (auto& v: values) {
        os << v << ", ";
    }
    return os;
}

// Раскладывает точки по
template <class Container>
void bucket_sort(const Container& points, std::span<const double> sections,
    std::span<std::list<Vector3d>> buckets, int coord) {

    if (points.empty()) { return; }

#ifdef ZEPHYR_ASSERTS
    double min_coord = std::numeric_limits<double>::max();
    double max_coord = std::numeric_limits<double>::lowest();
    for (const auto& p: points) {
        min_coord = std::min(min_coord, p[coord]);
        max_coord = std::max(max_coord, p[coord]);
    }
    if (min_coord < sections.front() || min_coord > sections.back()) {
        std::cout << "min_coord: " << min_coord << "\n";
        std::cout << "max_coord: " << max_coord << "\n";
        std::cout << "real secs: " << sections.front() << " " << sections.back() << "\n";
        throw std::runtime_error("wrong assumption bucket_sort #1");
    }
#endif

    if (sections.size() < 2 || sections.size() != buckets.size() + 1) {
        throw std::runtime_error("wrong assumption bucket_sort #2");
    }

    for (const auto& p: points) {
        int i = std::ranges::lower_bound(sections, p[coord]) - sections.begin();
        if (i < 1 || i >= sections.size()) {
            throw std::runtime_error("Blocks::split_exact: bad range");
        }
        buckets[i - 1].push_back(p);
    }
}

using zephyr::utils::mpi;

// Формальный тип для MPI-группы
struct group_t {
#ifdef ZEPHYR_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
#endif
    group_t() = default;
};

// Локальное число точек
inline int sum(const std::vector<int>& counts) {
    return std::accumulate(counts.begin(), counts.end(), 0);
}

inline std::vector<double> merge_sorted(const std::vector<std::vector<double>>& vectors) {
    if (vectors.empty()) return {};

    // Подсчитываем общий размер
    size_t total_size = 0;
    for (const auto& v : vectors) {
        total_size += v.size();
    }

    std::vector<double> result;
    result.reserve(total_size);

    // Индексы для каждого вектора
    std::vector<size_t> indices(vectors.size(), 0);

    // Слияние всех отсортированных векторов
    while (true) {
        double min_val = std::numeric_limits<double>::max();
        int min_idx = -1;

        for (size_t i = 0; i < vectors.size(); ++i) {
            if (indices[i] < vectors[i].size() && vectors[i][indices[i]] < min_val) {
                min_val = vectors[i][indices[i]];
                min_idx = i;
            }
        }

        if (min_idx == -1) break;

        result.push_back(min_val);
        indices[min_idx]++;
    }

    return result;
}

// На каждом MPI-процессе в группе (!) есть отсортированный массив, собрать на корневом процессе
// группы единый отсортированный массив.
inline std::vector<double> merge_sorted(const std::vector<double>& points, group_t group) {
#ifndef ZEPHYR_MPI
    return points;
#else
    // Ранг и размер группы
    int rank, size;
    MPI_Comm_rank(group.comm, &rank);
    MPI_Comm_size(group.comm, &size);

    if (size == 1) {
        return points;
    }

    constexpr int root = 0; // ранги в группе с нуля

    const int n_points = points.size();
    std::vector<int> all_sizes(size, 0);

    MPI_Gather(&n_points, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, root, group.comm);

    if (rank == root) {
        std::vector<int> offsets(size, 0);
        int total_size = 0;
        for (int i = 0; i < size; ++i) {
            offsets[i] = total_size;
            total_size += all_sizes[i];
        }

        std::vector<double> flat_data(total_size);
        MPI_Gatherv(points.data(), n_points, MPI_DOUBLE,
                    flat_data.data(), all_sizes.data(), offsets.data(), MPI_DOUBLE,
                    root, group.comm);

        std::vector<std::vector<double>> result(size);
        for (int i = 0; i < size; ++i) {
            result[i].assign(flat_data.begin() + offsets[i],
                             flat_data.begin() + offsets[i] + all_sizes[i]);
        }

        return merge_sorted(result);
    }
    else {
        MPI_Gatherv(points.data(), n_points, MPI_DOUBLE,
                    nullptr, nullptr, nullptr, MPI_DOUBLE,
                    root, group.comm);
        return {};
    }
#endif
}

class Buckets {
public:
    // points -- Полный набор точек на данном процессе
    // barriers -- Исходные разделяющие барьеры
    // n_parts -- Число подблоков разбиения для каждого исходного блока
    Buckets(const std::vector<Vector3d>& points, const std::vector<double>& barriers, const int n_parts) {
        z_assert(barriers.size() > 1, "Bad barriers size");
        const int n_blocks = static_cast<int>(barriers.size()) - 1;

        // Сечения на подблоки
        sections.reserve(n_parts * n_blocks + 1);

        // Заполним подразбиение линейно (хорошо бы проверить варианты со сгущением)
        sections.push_back(barriers[0]);
        for (int i = 0; i < n_blocks; ++i) {
            for (int k = 0; k < n_parts; ++k) {
                sections.push_back(barriers[i] + (barriers[i + 1] - barriers[i]) * (k + 1) / n_parts);
            }
        }

        // Раскидаем точки по подблокам
        buckets.resize(n_parts * n_blocks);
        bucket_sort(points, sections, buckets, 0);

        // Заполним массив с размерами
        counts.reserve(buckets.size());
        for (const auto& b: buckets) {
            counts.push_back(b.size());
        }

        if (mpi::master()) {
            std::cout << "barriers: " << barriers << "\n";
            std::cout << "sections: " << sections << "\n";
            std::cout << "counts:   " << counts << "\n";
        }

        // Число точек равно сумме counts.
        z_assert(points.size() == sum(counts), "Bad sum #195");
    }

    // Выполнить подразбитие интересных блоков на n_parts частей
    void specify(const std::vector<int>& k_stat, int n_parts) {
        // Интересные блоки
        auto subblock_idx = interesting_blocks_(k_stat);

        // Подразбить интересные блоки
        split_(subblock_idx, n_parts);

        if (mpi::master()) {
            std::cout << "sections: " << sections << "\n";
            std::cout << "counts:   " << counts << "\n";
        }
    }

    std::vector<double> approx_barriers(const std::vector<int>& k_stat) const {
        // Интересные блоки
        const auto subblock_idx_set = interesting_blocks_(k_stat);

        std::vector<int> subblock_idx;
        subblock_idx.reserve(subblock_idx_set.size());
        for (const auto& idx: subblock_idx_set) {
            subblock_idx.push_back(idx);
        }

        std::vector<double> barriers;

        barriers.reserve(subblock_idx.size());

        auto comm_summ = accumulate_counts_();

        for (int i = 0; i < subblock_idx.size(); ++i) {
            int idx = subblock_idx[i];
            int k = k_stat[i];
            // Простая медиана
            //double x = 0.5 * (sections[idx] + sections[idx + 1]);
            // линейная аппроксимация
            double x = sections[idx] + (sections[idx + 1] - sections[idx]) * (k - comm_summ[idx]) / (comm_summ[idx + 1] - comm_summ[idx]);
            barriers.push_back(x);
        }
        return barriers;
    }

    std::vector<double> exact_barriers(const std::vector<int>& k_stat) const {
        // Интересные блоки
        const auto subblock_idx_set = interesting_blocks_(k_stat);

        std::vector<int> subblock_idx;
        subblock_idx.reserve(subblock_idx_set.size());
        for (const auto& idx: subblock_idx_set) {
            subblock_idx.push_back(idx);
        }

        std::vector<double> barriers;
        barriers.reserve(subblock_idx.size());

        auto comm_summ = accumulate_counts_();

        for (int i = 0; i < subblock_idx.size(); ++i) {
            int idx = subblock_idx[i];
            std::vector<double> points;
            points.reserve(counts[idx]);
            for (const auto& p: buckets[idx]) {
                points.push_back(p.x());
            }
            std::ranges::sort(points);

            std::cout << "Subblock " << idx << "; points: ";
            for (auto& p: points) {
                std::cout << p << ", ";
            }
            std::cout << "\n";

            // TODO: MPI VERSION!!!!
            if (!mpi::single()) {
                throw std::runtime_error("OH NO NO NO");
            }
            //std::vector<double> global_array = merge_sorted(points);
            //points = std::move(global_array);

            z_assert(k_stat[i] >= comm_summ[idx] && k_stat[i] <= comm_summ[idx + 1], "Bad k-stat");

            double x;
            if (mpi::master()) {
                int loc_idx = k_stat[i] - comm_summ[idx] - 1;
                z_assert(0 <= loc_idx && loc_idx < points.size(), "Bad loc index exact ORB");
                x = (1.0 + 1.1235423e-14) * points[loc_idx];

                std::cout << k_stat[i] << " " << comm_summ[idx] << " " << comm_summ[idx + 1] << " " << comm_summ[idx + 1] - comm_summ[idx] << "; " << points.size() << "\n";
                std::cout << loc_idx << ", " << x << "\n";
                std::cout << "Bunch: ";
                for (auto& p: points) {
                    std::cout << p << ", ";
                }
                std::cout << std::endl;
            }
            mpi::broadcast(0, x);
            barriers.push_back(x);
        }

        return barriers;
    }

public:

    // Кумулятивная гистограмма по всем процессам
    std::vector<int> accumulate_counts_() const {
        std::vector<int> comm_summ(counts.size() + 1);
        comm_summ[0] = 0;
        for (int i = 0; i < counts.size(); ++i) {
            comm_summ[i + 1] = comm_summ[i] + counts[i];
        }
        // Кумулятивная сумма по подблокам со всех процессов
        comm_summ = mpi::sum(comm_summ);
        if (mpi::master()) {
            std::cout << "comm_summ: " << comm_summ << "\n";
        }

#ifdef ZEPHYR_ASSERTS
        const int n_points = mpi::sum(sum(counts));
        if (comm_summ.back() != n_points) {
            throw std::runtime_error("Blocks::split_exact: bad comm_summ");
        }
#endif
        return comm_summ;
    }

    // Индексы подблоков, куда попадает k-статистика
    // comm_summ -- Кумулятивное число точек по подблокам
    // k_stat -- Индексы k-ой статистики
    std::set<int> interesting_blocks_(const std::vector<int>& k_stat) const {
        // Кумулятивная сумма точек по подблокам
        std::vector<int> comm_summ = accumulate_counts_();

        if (mpi::master()) {
            std::cout << "comm summ: " << comm_summ << "\n";
        }

        std::set<int> subblock_idx;
        for (int i = 0; i < k_stat.size(); ++i) {
            const int m = std::ranges::lower_bound(comm_summ, k_stat[i]) - comm_summ.begin();
            if (m < 1 || m >= comm_summ.size()) {
                throw std::runtime_error("Blocks::split_exact: bad range");
            }
            subblock_idx.insert(m - 1);
        }

        if (mpi::master()) {
            std::cout << "subblock_num: ";
            for (int c: subblock_idx) {
                std::cout << c << ", ";
            }
            std::cout << "\n";
        }
        return subblock_idx;
    }

    // Подразбить интересные блоки, остальные склеить
    void split_(const std::set<int>& subblock_idx, int n_parts) {
        // Теперь сложно, надо склеить subblocks и nums
        std::vector<double> prev_sections = std::move(sections);
        std::vector<int> prev_counts = std::move(counts);
        std::vector<std::list<Vector3d>> prev_buckets = std::move(buckets);

        // Делаем склейку путем заполнения по новой
        const int int_blocks = subblock_idx.size();
        sections.clear();
        sections.reserve(int_blocks * n_parts + int_blocks + 2);

        counts.clear();
        counts.reserve(int_blocks * n_parts + int_blocks + 1);

        buckets.clear();
        buckets.reserve(int_blocks * n_parts + int_blocks + 1);

        sections.push_back(prev_sections[0]);

        bool prev_interesting = true;
        for (int i = 0; i < prev_buckets.size(); ++i) {
            // Номер нового блока для добавления
            int curr = static_cast<int>(buckets.size());

            if (subblock_idx.contains(i)) {
                // Интересный блок.
                // Добавим линейно распределенные узлы и пустые корзины
                for (int k = 0; k < n_parts; ++k) {
                    double xi = prev_sections[i] + (prev_sections[i + 1] - prev_sections[i]) * (k + 1) / n_parts;
                    sections.push_back(xi);
                    buckets.emplace_back();
                }

                // распределим точки по корзинам
                bucket_sort(prev_buckets[i], std::span(&sections[curr], n_parts + 1), std::span(&buckets[curr], n_parts), 0);
                for (int k = 0; k < n_parts; ++k) {
                    counts.push_back(static_cast<int>(buckets[curr + k].size()));
                }
                prev_interesting = true;

#ifdef ZEPHYR_ASSERTS
                // Новое разбиение в сумме дает старое число точек
                int prev_count = prev_counts[i];
                int new_count = std::accumulate(counts.data() + curr, counts.data() + curr + n_parts, 0);
                if (prev_count != new_count) {
                    throw std::runtime_error("Blocks::split_exact: bad counts");
                }
#endif
            }
            else {
                // неинтересный блок
                if (prev_interesting) {
                    // Предыдущий интересный, копируем данные
                    sections.push_back(prev_sections[i + 1]);
                    // Не переносим список точек, это же неинтересный блок
                    buckets.emplace_back();
                    // Но число точек запоминаем
                    counts.push_back(prev_counts[i]);
                }
                else {
                    // Предыдущий неинтересный, склейка с ним
                    if (buckets.empty()) {
                        throw std::runtime_error("Impossible exact ORB error");
                    }
                    sections.back() = prev_sections[i + 1];
                    counts.back() += prev_counts[i];
                }
                prev_interesting = false;
            }
        }
    }

    // Границы для гистограммы, одинаковые на всех процессах, охватывают всю область целиком
    std::vector<double> sections;

    // Число точек в каждом отсеке гистограммы. Считаются только локальные точки!
    std::vector<int> counts;

    // Список точек в каждом отсеке гистограммы. Список ОТСУТСТВУЕТ для неинтересной
    // части гистограммы, в этом случае buckets[i].size() != counts[i].
    //
    std::vector<std::list<Vector3d>> buckets;
};

void Blocks::split_exact(const std::vector<Vector3d>& points_arg) {
    std::vector<Vector3d> points(points_arg.size());
    for (int i = 0; i < points_arg.size(); ++i) {
        points[i] = m_map(points_arg[i]);
    }

    // Число точек
    const int n_points = mpi::sum(static_cast<int>(points.size()));

    // k-порядковая статистика (индексы точек для разбиения)
    std::vector<int> k_stat(m_nx - 1);
    for (int i = 0; i < m_nx - 1; ++i) {
        k_stat[i] = (i + 1) * n_points / m_nx;
    }
    if (mpi::master()) {
        std::cout << "k stat: " << k_stat << "\n";
    }

    std::vector<double> barriers(m_x_coords.size());
    for (int i = 0; i < m_x_coords.size(); ++i) {
        barriers[i] = m_x_coords[i].val;
    }

    Buckets buckets(points, barriers, 5);
    buckets.specify(k_stat, 8);
    buckets.specify(k_stat, 4);
    buckets.specify(k_stat, 4);

    //auto xnew = buckets.approx_barriers(k_stat);
    auto xnew = buckets.exact_barriers(k_stat);

    double x = buckets.sections[2];
    int less = 0;
    for (auto& p: points) {
        if (p.x() < x) {
            ++less;
        }
    }
    less = mpi::sum(less);
    auto comm_summ = buckets.accumulate_counts_();
    std::cout << "count less: " << less << "; " << comm_summ[2] << "\n";

    z_assert(xnew.size() + 2 == m_x_coords.size(), "OLOL #123");
    for (int i = 0; i < xnew.size(); ++i) {
        m_x_coords[i + 1] = xnew[i];
    }
}

} // namespace zephyr::mesh::decomp