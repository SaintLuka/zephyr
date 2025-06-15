#include <zephyr/utils/mpi.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/decomp/ORB.h>

namespace zephyr::mesh::decomp {

using namespace zephyr::utils;

ORB::ORB(Box domain, const std::string &type, int size, const ORB::params &p)
    : Decomposition(size) {
    m_newton   = p.newton;
    m_mobility = p.mobility;

    m_blocks = Blocks(domain, type, m_size);
}

ORB::ORB(Box domain, const std::string& type, int size, int nx,
    const params& p)
    : Decomposition(size) {
    m_newton   = p.newton;
    m_mobility = p.mobility;

    m_blocks = Blocks(domain, type, size, nx);
    m_size = m_blocks.size();
}

ORB::ORB(Box domain, const std::string& type, int size,
    const std::vector<int>& ny, const params& p)
    : Decomposition(size) {
    m_newton   = p.newton;
    m_mobility = p.mobility;

    m_blocks = Blocks(domain, type, size, ny);
    m_size = m_blocks.size();
}

ORB::ORB(Box domain, const utils::Json& config)
    : Decomposition(mpi::size()) {
    params p;
    m_newton = p.newton;
    m_mobility = p.mobility;

    if (config["mobility"]) {
        m_mobility = config["mobility"].as<double>();
    }
    if (config["newton"]) {
        m_newton = config["newton"].as<double>();
    }

    std::string type;
    if (config["type"]) {
        type = config["type"].as<std::string>();

        // Попытка разбить двумерную сетку по оси Z приводит к ошибке
        if (domain.is_2D()) {
            type.erase(type.find('Z'), 1);
        }
    } else {
        type = domain.is_2D() ? "XY" : "XYZ";
    }

    if (config["proc_nx"]) {
        m_blocks = Blocks(domain, type, m_size);
    } else {
        int nx = config["nx"].as<int>();
        m_blocks = Blocks(domain, type, m_size, nx);
    }
}

int ORB::rank(const Vector3d& v) const {
    return m_blocks.rank(v);
}

int ORB::rank(QCell& elem) const {
    return m_blocks.rank(elem.center());
}

int ORB::rank(AmrStorage::Item& elem) const {
    return m_blocks.rank(elem.center);
}

inline double between(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

void ORB::use_newton(bool val) {
    m_newton = val;
}

void ORB::set_mobility(double val) {
    m_mobility = between(val, 0.0, 0.45);
}

double ORB::mobility() const {
    return m_mobility;
}

bool ORB::newton() const {
    return m_newton;
}

void ORB::balancing(const std::vector<double> &w) {
    if (m_newton) {
        m_blocks.balancing_newton(w, m_mobility);
    }
    else {
        m_blocks.balancing_simple(w, m_mobility);
    }
}

} // namespace zephyr::mesh::decomp