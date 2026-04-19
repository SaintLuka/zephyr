#include <iostream>
#include <cmath>
#include <limits>

#include <algorithm>
#include <stdexcept>

#include <zephyr/mesh/decomp/orb/transform.h>

#define USE_ROTATION 1

namespace zephyr::mesh::decomp {

#if USE_ROTATION

// Небольшие синусы для углов поворота
static const double S1 = 0.00795;
static const double S2 = 0.00735;

// Соответствующие косинусы
static const double C1 = std::sqrt(1.0 - S1 * S1);
static const double C2 = std::sqrt(1.0 - S2 * S2);

// Матрица поворота (2D)
static const double R2[3][3] = {
    {+C1, S1, 0.0},
    {-S1, C1, 0.0},
    {0.0, 0.0, 1.0}
};

// Матрица поворота (3D)
static const double R3[3][3] = {
    {+C1, S1, 0.0},
    {S1*C2, C1*C2, -S2},
    {S2*S1, C1*S2, C2}
};

inline Vector3d rotate_2d(const Vector3d& v) {
    return ((const geom::Matrix3d&) R2) * v;
}
inline Vector3d rotate_3d(const Vector3d& v) {
    return ((const geom::Matrix3d&) R3) * v;
}
inline Vector3d rotate_2d_inv(const Vector3d& v) {
    return ((const geom::Matrix3d&) R2).transpose() * v;
}
inline Vector3d rotate_3d_inv(const Vector3d& v) {
    return ((const geom::Matrix3d&) R3).transpose() * v;
}
#else
inline Vector3d rotate_2d(const Vector3d& v) { return v; }
inline Vector3d rotate_3d(const Vector3d& v) { return v; }
inline Vector3d rotate_2d_inv(const Vector3d& v) { return v; }
inline Vector3d rotate_3d_inv(const Vector3d& v) { return v; }
#endif

Vector3d map_xy (const Vector3d& v) { return rotate_2d({v.x(), v.y(), v.z()}); }
Vector3d map_xyz(const Vector3d& v) { return rotate_3d({v.x(), v.y(), v.z()}); }
Vector3d map_xz (const Vector3d& v) { return rotate_2d({v.x(), v.z(), v.y()}); }
Vector3d map_xzy(const Vector3d& v) { return rotate_3d({v.x(), v.z(), v.y()}); }
Vector3d map_yx (const Vector3d& v) { return rotate_2d({v.y(), v.x(), v.z()}); }
Vector3d map_yxz(const Vector3d& v) { return rotate_3d({v.y(), v.x(), v.z()}); }
Vector3d map_yz (const Vector3d& v) { return rotate_2d({v.y(), v.z(), v.x()}); }
Vector3d map_yzx(const Vector3d& v) { return rotate_3d({v.y(), v.z(), v.x()}); }
Vector3d map_zx (const Vector3d& v) { return rotate_2d({v.z(), v.x(), v.y()}); }
Vector3d map_zxy(const Vector3d& v) { return rotate_3d({v.z(), v.x(), v.y()}); }
Vector3d map_zy (const Vector3d& v) { return rotate_2d({v.z(), v.y(), v.x()}); }
Vector3d map_zyx(const Vector3d& v) { return rotate_3d({v.z(), v.y(), v.x()}); }


Vector3d inv_xy (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.x(), w.y(), w.z()}; }
Vector3d inv_xyz(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.x(), w.y(), w.z()}; }
Vector3d inv_xz (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.x(), w.z(), w.y()}; }
Vector3d inv_xzy(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.x(), w.z(), w.y()}; }
Vector3d inv_yx (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.y(), w.x(), w.z()}; }
Vector3d inv_yxz(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.y(), w.x(), w.z()}; }
Vector3d inv_yz (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.z(), w.x(), w.y()}; }
Vector3d inv_yzx(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.z(), w.x(), w.y()}; }
Vector3d inv_zx (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.y(), w.z(), w.x()}; }
Vector3d inv_zxy(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.y(), w.z(), w.x()}; }
Vector3d inv_zy (const Vector3d& v) { Vector3d w = rotate_2d_inv(v); return {w.z(), w.y(), w.x()}; }
Vector3d inv_zyx(const Vector3d& v) { Vector3d w = rotate_3d_inv(v); return {w.z(), w.y(), w.x()}; }

Vector3d map_rp(const Vector3d& v) {
    return {std::sqrt(v.x() * v.x() + v.y() * v.y()), std::atan2(v.y(), v.x()), v.z()};
}
Vector3d map_pr(const Vector3d& v) {
    return {std::atan2(v.y(), v.x()), std::sqrt(v.x() * v.x() + v.y() * v.y()), v.z()};
}
Vector3d inv_rp(const Vector3d& v) {
    return {v.x() * std::cos(v.y()), v.x() * std::sin(v.y()), v.z()};
}
Vector3d inv_pr(const Vector3d& v) {
    return {v.y() * std::cos(v.x()), v.y() * std::sin(v.x()), v.z()};
}

Limits::Limits() {
    constexpr double neg_infty = -std::numeric_limits<double>::infinity();
    constexpr double pos_infty = +std::numeric_limits<double>::infinity();

    min = {neg_infty, neg_infty, neg_infty};
    max = {pos_infty, pos_infty, pos_infty};
}

Limits::Limits(std::string type)
        : Limits() {
    if (type.size() != 2 && type.size() != 3) {
        throw std::runtime_error("Limits::Limits error");
    }
    for (int i = 0; i < type.size(); ++i) {
        set(type[i], i);
    }
}

Limits::Limits(const Vector3d& vmin, const Vector3d& vmax)
        :min(vmin), max(vmax) {
}

void Limits::set(char c, int axes) {
    switch (c) {
        case 'R':
            set_radius(axes);
            break;
        case 'P':
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

void Transform::init_(const std::string &type) {
    auto available = available_types();
    const auto idx = std::ranges::find(available, type) - available.begin();

    if (idx < 0 || idx >= available.size()) {
        std::cerr << "Unknown orthogonal decomposition type " << type << ".\n";
        std::cerr << "Available Options: ";
        for (size_t i = 0; i < available.size() - 1; ++i) {
            std::cerr << "\"" << available[i] << "\", ";
        }
        std::cerr << available.back() << ".\n";
        throw std::runtime_error("Unknown orthogonal decomposition type '" + type + "'");
    }

    const std::vector<map_function> mapping = {
        map_xy,  map_xz,  map_yx,  map_yz,  map_zx,  map_zy,
        map_xyz, map_xzy, map_yxz, map_yzx, map_zxy, map_zyx,
        map_rp, map_pr
    };

    const std::vector<map_function> inverse = {
        inv_xy,  inv_xz,  inv_yx,  inv_yz,  inv_zx,  inv_zy,
        inv_xyz, inv_xzy, inv_yxz, inv_yzx, inv_zxy, inv_zyx,
        inv_rp, inv_pr
    };

    if (mapping.size() != available.size() ||
        inverse.size() != available.size()) {
        throw std::runtime_error("Transform: something wrong");
        }

    m_type = available[idx];
    m_mapping = mapping[idx];
    m_inverse = inverse[idx];
}

void Transform::update_limits(const Box& domain) {
    if (m_type == "RP" || m_type == "PR") {
        throw std::runtime_error("Transform::init_box_ polar not supported");
    }

    // Строим box, который захватывает углы m_mapping(domain)
    if (domain.is_2D()) {
        m_box = Box::Empty(2);
        for (double x: {domain.vmin.x(), domain.vmax.x()}) {
            for (double y: {domain.vmin.y(), domain.vmax.y()}) {
                Vector3d v = {x, y, 0.0};
                m_box.capture(m_mapping(v));
            }
        }
        m_box.extend(0.01, 0.01);
        m_box.vmin.z() = -std::numeric_limits<double>::infinity();
        m_box.vmax.z() = +std::numeric_limits<double>::infinity();
    }
    else {
        m_box = Box::Empty(3);
        for (double x: {domain.vmin.x(), domain.vmax.x()}) {
            for (double y: {domain.vmin.y(), domain.vmax.y()}) {
                for (double z: {domain.vmin.z(), domain.vmax.z()}) {
                    Vector3d v = {x, y, z};
                    m_box.capture(m_mapping(v));
                }
            }
        }
        m_box.extend(0.01, 0.01, 0.01);
    }
}

//Transform::Transform(const std::string& type) {
//    init_(type);
//}

std::string full_type(const Box& domain, const std::string& type) {
    if (domain.is_2D()) {
        // Нормальные двумерные отображения
        if (type == "XY" || type == "YX" || type == "RP" || type == "PR") {
            return type;
        }

        // Одномерные дополняем до двумерных
        if (type == "X") { return "XY"; }
        if (type == "Y") { return "YX"; }
        if (type == "R") { return "RP"; }
        if (type == "P") { return "PR"; }

        throw std::runtime_error("Bad decomposition type '" + type + "' for 2D domain");
    }
    else {
        // Нормальные трёхмерные отображения
        if (type == "XYZ" || type == "XZY" || type == "YXZ" ||
            type == "YZX" || type == "ZXY" || type == "ZYX") {
            return type;
            }

        // Одномерные и двумерные дополняем до трёхмерных
        if (type == "X" || type == "XY") { return "XYZ"; }
        if (type == "Y" || type == "YZ") { return "YZX"; }
        if (type == "Z" || type == "ZX") { return "ZXY"; }

        if (type == "YX") { return "YXZ"; }
        if (type == "ZY") { return "ZYX"; }
        if (type == "XZ") { return "XZY"; }

        throw std::runtime_error("Bad decomposition type '" + type + "' for 3D domain");
    }
}

Transform::Transform() {
    init_("XYZ");
    m_box = Box::Infinite(3);
}

Transform::Transform(const Box& domain, const std::string& type) {
    init_(full_type(domain, type));
    update_limits(domain);
}

const std::string& Transform::type() const {
    return m_type;
}

int Transform::dim() const {
    return m_type.size();
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

std::vector<std::string> Transform::available_types() {
    return {
        // Двумерная декартова декомпозиция
        "XY", "XZ", "YX", "YZ", "ZX", "ZY",

        // Трехмерная декартова декомпозиция
        "XYZ", "XZY", "YXZ", "YZX", "ZXY", "ZYX",

        // Двумерная полярная декомпозиция
        "RP", "PR",
    };
}

void Transform::check() {
    Box domain({0.0, 0.0, 0.0}, {3.0, 3.0, 3.0});
    Vector3d vec0 = {1.232, 0.823, 2.123};
    for (auto& type: available_types()) {
        Transform t(domain, type);
        Vector3d vec1 = t.mapping(vec0);
        Vector3d vec2 = t.inverse(vec1);

        if ((vec2 - vec0).norm() > 1.0e-13) {
            std::cerr << "vec0: {" << vec0.x()  << ", " << vec0.y()  << ", " << vec0.z()  << "}\n";
            std::cerr << "vec1: {" << vec1.x()  << ", " << vec1.y()  << ", " << vec1.z()  << "}\n";
            std::cerr << "vec2: {" << vec2.x()  << ", " << vec2.y()  << ", " << vec2.z()  << "}\n";
            throw std::runtime_error("Wrong transformation for '" + type + "'");
        }
    }
}

} // namespace zephyr::mesh::decomp