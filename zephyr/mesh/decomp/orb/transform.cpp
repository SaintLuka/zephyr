#include <iostream>
#include <cmath>
#include <limits>

#include <algorithm>
#include <stdexcept>

#include <zephyr/mesh/decomp/orb/transform.h>

namespace zephyr::mesh::decomp {

const std::vector<std::string> Transform::available = {
        // Одномерная декартова декомпозиция
        "X", "Y", "Z",

        // Двумерная декартова декомпозиция
        "XY", "XZ", "YX", "YZ", "ZX", "ZY",

        // Трехмерная декартова декмпозиция
        "XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX",

        // Для полярной декомпозиции пока не продуманы
        // алгоритмы поиска соседей и окрестностей

        // Одномерная полярная декомпозиция
        // "R", "P",

        // Двумерная полярная декомпозиция
        // "RP", "PR",
};

Vector3d rotate(const Vector3d& v) {
    const double s = 0.0; // std::sin(1.0/129.0);
    const double c = 1.0; // std::cos(1.0/129.0);

    return {
            c * c * v.x() + c * s * v.z() + s * v.y(),
            -c * s * v.x() - s * s * v.z() + c * v.y(),
            c * v.z() - s * v.x()
    };
}

Vector3d map_xy(const Vector3d& v) { return rotate({v.x(), v.y(), v.z()}); }

Vector3d map_xz(const Vector3d& v) { return rotate({v.x(), v.z(), v.y()}); }

Vector3d map_yx(const Vector3d& v) { return rotate({v.y(), v.x(), v.z()}); }

Vector3d map_yz(const Vector3d& v) { return rotate({v.y(), v.z(), v.x()}); }

Vector3d map_zx(const Vector3d& v) { return rotate({v.z(), v.x(), v.y()}); }

Vector3d map_zy(const Vector3d& v) { return rotate({v.z(), v.y(), v.x()}); }

Vector3d map_rp(const Vector3d& v) {
    return {std::sqrt(v.x() * v.x() + v.y() * v.y()), std::atan2(v.y(), v.x()), v.z()};
}

Vector3d map_pr(const Vector3d& v) {
    return {std::atan2(v.y(), v.x()), std::sqrt(v.x() * v.x() + v.y() * v.y()), v.z()};
}

Vector3d inv_xy(const Vector3d& v) { return {v.x(), v.y(), v.z()}; }

Vector3d inv_xz(const Vector3d& v) { return {v.x(), v.z(), v.y()}; }

Vector3d inv_yx(const Vector3d& v) { return {v.y(), v.x(), v.z()}; }

Vector3d inv_yz(const Vector3d& v) { return {v.z(), v.x(), v.y()}; }

Vector3d inv_zx(const Vector3d& v) { return {v.y(), v.z(), v.x()}; }

Vector3d inv_zy(const Vector3d& v) { return {v.z(), v.y(), v.x()}; }

Vector3d inv_rp(const Vector3d& v) {
    return {v.x() * std::cos(v.y()), v.x() * std::sin(v.y()), v.z()};
}

Vector3d inv_pr(const Vector3d& v) {
    return {v.y() * std::cos(v.x()), v.y() * std::sin(v.x()), v.z()};
}

Limits::Limits() {
    const double neg_infty = -std::numeric_limits<double>::infinity();
    const double pos_infty = +std::numeric_limits<double>::infinity();

    min = {neg_infty, neg_infty, neg_infty};
    max = {pos_infty, pos_infty, pos_infty};
}

Limits::Limits(std::string type)
        : Limits() {

    if (type == "R") {
        type = "RP";
    }
    if (type == "P") {
        type = "PR";
    }
    if (type.length() > 3) {
        throw std::runtime_error("Limits::Limits error");
    }
    for (size_t i = 0; i < type.size(); ++i) {
        set(type[i], i);
    }
}

Limits::Limits(const Vector3d& vmin, const Vector3d& vmax)
        :min(vmin), max(vmax) {
}

void Limits::set(char c, int axes) {
    switch (std::tolower(c)) {
        case 'r':
            set_radius(axes);
            break;
        case 'p':
            set_angle(axes);
            break;
        default:
            set_coord(axes);
            break;
    }
}

void Limits::set_coord(int axes) {
    min[axes] = -std::numeric_limits<double>::infinity();
    max[axes] = +std::numeric_limits<double>::infinity();
}

void Limits::set_radius(int axes) {
    min[axes] = 0.0;
    max[axes] = std::numeric_limits<double>::infinity();
}

void Limits::set_angle(int axes) {
    min[axes] = -M_PI;
    max[axes] = +M_PI;
}

Transform::Transform() {
    m_type = "XY";
    m_limits = Limits();
    m_mapping = map_xy;
    m_inverse = inv_xy;
}

Transform::Transform(const std::string& type) {
    auto it = std::find(available.begin(), available.end(), type);
    size_t idx = std::distance(available.begin(), it);

    if (idx >= available.size()) {
        std::cerr << "Unknown orthogonal decomposition type " << type << ".\n";
        std::cerr << "Available Options: ";
        for (size_t i = 0; i < available.size() - 1; ++i) {
            std::cerr << "\"" << available[i] << "\", ";
        }
        std::cerr << available.back() << ".\n";
        throw std::runtime_error("Unknown orthogonal decomposition type '" + type + "'");
    }

    const std::vector<map_function> mapping = {
            map_xy, map_yx, map_zx,
            map_xy, map_xz, map_yx, map_yz, map_zx, map_zy,
            map_xy, map_xz, map_yx, map_yz, map_zx, map_zy,
            //map_rp, map_pr,
            //map_rp, map_pr
    };

    const std::vector<map_function> inverse = {
            inv_xy, inv_yx, inv_zx,
            inv_xy, inv_xz, inv_yx, inv_yz, inv_zx, inv_zy,
            inv_xy, inv_xz, inv_yx, inv_yz, inv_zx, inv_zy,
            //inv_rp, inv_pr,
            //inv_rp, inv_pr
    };

    if (mapping.size() != available.size() ||
        inverse.size() != available.size()) {
        throw std::runtime_error("Transform: something wrong");
    }

    m_type = available[idx];
    m_limits = Limits(m_type);
    m_mapping = mapping[idx];
    m_inverse = inverse[idx];
}

const std::string& Transform::type() const {
    return m_type;
}

size_t Transform::dimension() const {
    return m_type.size();
}

const Limits& Transform::limits() const {
    return m_limits;
}

Vector3d Transform::operator()(const Vector3d& vec) const {
    return m_mapping(vec);
}

Vector3d Transform::mapping(const Vector3d& vec) const {
    return m_mapping(vec);
}

Vector3d Transform::inverse(const Vector3d& vec) const {
    return m_inverse(vec);
}

void Transform::check() {
    Vector3d vec0(1.0, 2.0, 3.0);
    for (auto& type: available) {
        Transform t(type);
        Vector3d vec1 = t.mapping(vec0);
        Vector3d vec2 = t.inverse(vec1);

        double err = std::fabs(vec2.x() - vec0.x() ) +
                     std::fabs(vec2.y() - vec0.y() ) +
                     std::fabs(vec2.z() - vec0.z() );

        if (err > 1.0e-13) {
            std::cerr << "vec0: (" << vec0.x()  << ", " << vec0.y()  << ", " << vec0.z()  << ")\n";
            std::cerr << "vec1: (" << vec1.x()  << ", " << vec1.y()  << ", " << vec1.z()  << ")\n";
            std::cerr << "vec2: (" << vec2.x()  << ", " << vec2.y()  << ", " << vec2.z()  << ")\n";
            throw std::runtime_error("Wrong transformation for '" + type + "'");
        }
    }
}

} // namespace zephyr::mesh::decomp