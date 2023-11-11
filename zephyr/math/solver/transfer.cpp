#include <zephyr/math/solver/transfer.h>

#include <zephyr/geom/geom.h>
#include <zephyr/geom/polygon.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/geom/primitives/basic_face.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;

static const Transfer::State U = Transfer::datatype();

Transfer::State Transfer::datatype() {
    return {};
}

Transfer::Transfer() {
    m_CFL = 0.5;
    m_dt = 1.0e+300;
}

double Transfer::CFL() const {
    return m_CFL;
}

void Transfer::set_CFL(double C) {
    m_CFL = std::max(0.0, std::min(C, 1.0));
}

double Transfer::dt() const {
    return m_dt;
}

Vector3d Transfer::velocity(const Vector3d& c) const {
    return Vector3d::UnitX();
}

double Transfer::compute_dt(EuCell &cell) {
    double max_area = 0.0;
    for (auto &face: cell.faces()) {
        max_area = std::max(max_area, face.area());
    }
    double dx = cell.volume() / max_area;
    return m_CFL * dx / velocity(cell.center()).norm();
}

double Transfer::compute_dt(EuMesh& mesh) {
    double dt = std::numeric_limits<double>::max();
    for (auto &cell: mesh) {
        dt = std::min(dt, compute_dt(cell));
    }
    return dt;
}

// Поток из пустой ячейки в ячейку с объемной долей a и шириной h
// verts -- скорость потока, dt -- шаг интегрирования
double flux_1D_0(double a, double h, double vs, double dt) {
    return std::min(0.0, dt * vs + (1.0 - a) * h);
}

// Поток из ячейки с объемной долей a и шириной h в полную ячейку
// verts -- скорость потока, dt -- шаг интегрирования
double flux_1D_1(double a, double h, double vs, double dt) {
    return std::min(a * h, dt * vs);
}

// Поток между ячейками с объемными долями a1, a2 вдоль нормали от a1 к a2
// Предполагается, что a1 < a2.
// h1, h2 -- длина ячеек вдоль нормали к грани
// as -- объемная доля на грани
// verts -- скорость вдоль нормали грани
// dt -- шаг интегрирования по времени
double flux_2D_less(double a1, double a2, double h1, double h2, double as, double vs, double dt) {
#if 0
    double F_in = 0.0;
    double F_ex = 0.0;
    if (as < 1.0) {
        F_ex = flux_1D_0((a2 - as) / (1.0 - as), h2, verts, dt);
    }
    if (as > 0.0) {
        F_in = flux_1D_1(a1 / as, h1, verts, dt);
    }
    return (1.0 - as) * F_ex + as * F_in;
#else
    double F_in = std::min(a1 * h1, as * dt * vs);
    double F_ex = std::min(0.0, (1.0 - as) * dt * vs + (1.0 - a2) * h2);
    return F_in + F_ex;
#endif
}

double flux_2D(double a1, double a2, double h1, double h2, double as, double vs, double dt) {
    if (a1 <= a2) {
        return +flux_2D_less(a1, a2, h1, h2, as, +vs, dt);
    }
    else {
        return -flux_2D_less(a2, a1, h2, h1, as, -vs, dt);
    }
}

double face_fraction(double a1, double a2) {
    double a_min = std::min(a1, a2);
    double a_max = std::max(a1, a2);
    double a_sig = a_min / (1.0 - (a_max - a_min));

    if (std::isnan(a_sig)) {
        // Случай a_min = 0, a_max = 1
        a_sig = 0.5; // ??
    }

    return std::max(a_min, std::min(a_sig, a_max));
}

void Transfer::fluxes_CRP(EuCell &cell, Direction dir) {
    auto &zc = cell(U);
    double a1 = zc.u1;

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        double a2 = neib(U).u1;

        double vs = velocity(face.center()).dot(face.normal());

        double as = face_fraction(a1, a2);

        double h1 = cell.volume() / face.area();
        double h2 = neib.volume() / face.area();

        double F = flux_2D(a1, a2, h1, h2, as, vs, m_dt);

        fluxes += F * face.area();
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

// Предполагаем verts > 0.0
double flux_VOF(double a, Vector3d& p, Vector3d& n, const AmrCell& cell, const AmrFace& face, double vs, double dt) {
    // Типа upwind
    Vector3d v1 = cell.vertices[face.vertices[0]];
    Vector3d v2 = cell.vertices[face.vertices[1]];
    Vector3d v3 = v1 - dt * vs * face.normal;
    Vector3d v4 = v2 - dt * vs * face.normal;

    // alpha - объемная доля, n - нормаль, face.area - площадь грани, xi = verts * dt.

    if (a < 1.0e-12) {
        return 0.0;
    } else {
        PolyQuad poly(v1, v2, v4, v3);
        return poly.clip_area(p, n);
    }
}

void Transfer::fluxes_VOF(EuCell &cell, Direction dir) {
    auto &zc = cell(U);

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();

        double vs = velocity(face.center()).dot(face.normal());

        // Типа upwind
        double F;
        if (vs > 0.0) {
            F = +flux_VOF(zc.u1, zc.p, zc.n, cell.geom(), face.geom(), vs, m_dt);
        }
        else {
            auto zn = neib(U);
            F = -flux_VOF(zn.u1, zn.p, zn.n, cell.geom(), face.geom(), vs, m_dt);
        }

        fluxes += F;
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

// Находит оптимальное деление грани a_sig, при котором поток совпадает с потоком по FOV
double best_face_fraction(double a1, double a2, double h1, double h2, double vs, double dt, double F_VOF) {
    double a_min = std::min(a1, a2);
    double a_max = std::max(a1, a2);

    // Ищем a_sig: func(a_sig) = 0
    auto func = [&](double a_sig) -> double {
        return flux_2D(a1, a2, h1, h2, a_sig, vs, dt) - F_VOF;
    };

    double f_min = func(a_min);
    double f_max = func(a_max);

    if (f_min * f_max >= 0.0) {
        // Метод дихотомии не применим
        // Возвращаем одно из крайних значений
        if (std::abs(f_min) < std::abs(f_max)) {
            return a_min;
        } else {
            return a_max;
        }
    }

    while (a_max - a_min > 1.0e-4) {
        double a_avg = 0.5 * (a_min + a_max);
        double f_avg = func(a_avg);

        if (f_avg * f_min < 0.0) {
            a_max = a_avg;
            f_max = f_avg;
        } else {
            a_min = a_avg;
            f_min = f_avg;
        }
    }

    double a_sig = 0.5 * (a_min + a_max);
    //double a_sig = a_min - (a_max - a_min) * f_min / (f_max - f_min);
    return a_sig;
}

void Transfer::fluxes_MIX(EuCell &cell, Direction dir) {
    auto &zc = cell(U);

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        auto zn = neib(U);

        double vs = velocity(face.center()).dot(face.normal());

        double F_VOF;
        if (vs > 0.0) {
            F_VOF = +flux_VOF(zc.u1, zc.p, zc.n, cell.geom(), face.geom(), vs, m_dt);
        }
        else {
            F_VOF = -flux_VOF(zn.u1, zn.p, zn.n, cell.geom(), face.geom(), vs, m_dt);
        }
        F_VOF /= face.area();

        double h1 = cell.volume() / face.area();
        double h2 = neib.volume() / face.area();

        double a1 = zc.u1;
        double a2 = zn.u1;

        // Хочу найти as, при котором flux_2D дает F_VOF
        //double a_sig = face_fraction(a1, a2);
        double a_sig = best_face_fraction(a1, a2, h1, h2, vs, m_dt, F_VOF);

        double F = flux_2D(a1, a2, h1, h2, a_sig, vs, m_dt);

        fluxes += F * face.area();
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

void Transfer::update(EuMesh &mesh, int ver) {
    switch (ver) {
        case 1:
            update_ver1(mesh);
            break;
        case 2:
            update_ver2(mesh);
            break;
        case 3:
            update_ver3(mesh);
            break;
        default:
            throw std::runtime_error("Unknown update version");
    }
}

void Transfer::update_ver1(EuMesh& mesh) {
    m_dt = compute_dt(mesh);

    // Считаем потоки
    static int counter = 0;
    for (auto cell: mesh) {
        if (counter % 2 == 0) {
            fluxes_CRP(cell, Direction::ANY);
        }
        else {
            fluxes_CRP(cell, Direction::ANY);
        }
    }
    ++counter;

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }
}

void Transfer::update_ver2(EuMesh& mesh) {
    // Определяем dt
    m_dt = compute_dt(mesh);

    compute_normals(mesh, 3);
    find_sections(mesh);

    // Считаем потоки
    static int counter = 0;
    for (auto cell: mesh) {
        if (counter % 2 == 0) {
            fluxes_VOF(cell, Direction::ANY);
        }
        else {
            fluxes_VOF(cell, Direction::ANY);
        }
    }
    ++counter;

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }

    compute_normals(mesh, 3);
    find_sections(mesh);
}

void Transfer::update_ver3(EuMesh& mesh) {
    // Определяем dt
    m_dt = compute_dt(mesh);

    compute_normals(mesh, 3);
    find_sections(mesh);

    // Считаем потоки
    static int counter = 0;
    for (auto cell: mesh) {
        if (counter % 2 == 0) {
            fluxes_MIX(cell, Direction::X);
        }
        else {
            fluxes_MIX(cell, Direction::Y);
        }
    }
    ++counter;

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }

    compute_normals(mesh, 3);
    find_sections(mesh);
}

void Transfer::compute_normal(EuCell& cell) {
    double uc = cell(U).u1;

    double n_x = 0.0;
    double n_y = 0.0;

    for (auto &face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();

        double un = neib(U).u1;

        Vector3d S = 0.5 * face.normal() * face.area();
        n_x -= (uc + un) * S.x();
        n_y -= (uc + un) * S.y();
    }

    cell(U).n.x() = n_x / cell.volume();
    cell(U).n.y() = n_y / cell.volume();

    cell(U).n.normalize();
}

void Transfer::smooth_normal(EuCell& cell) {
    double n_x = 4 * cell(U).n.x();
    double n_y = 4 * cell(U).n.y();

    for (auto &face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        n_x += neib(U).n.x();
        n_y += neib(U).n.y();
    }

    cell(U).n.x() = n_x;
    cell(U).n.y() = n_y;

    cell(U).n.normalize();
}

void Transfer::find_section(EuCell &cell) {
    auto poly = cell.geom().polygon();
    cell(U).p = poly.find_section(cell(U).n, cell(U).u1);
}

void Transfer::compute_normals(EuMesh& mesh, int smoothing) {
    for (auto cell: mesh) {
        compute_normal(cell);
    }

    for (int i = 0; i < smoothing; ++i) {
        for (auto cell: mesh) {
            smooth_normal(cell);
        }
    }
}

void Transfer::find_sections(EuMesh& mesh) {
    for (auto cell: mesh) {
        find_section(cell);
    }
}

void Transfer::set_flags(EuMesh& mesh) {
    for (auto cell: mesh) {
        double min_val = cell(U).u1;
        double max_val = cell(U).u1;

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            min_val = std::min(min_val, face.neib()(U).u1);
            max_val = std::max(max_val, face.neib()(U).u1);
        }

        if (max_val - min_val > 0.0) {
            cell.set_flag(1);
        }
        else {
            cell.set_flag(-1);
        }
    }
}

Distributor Transfer::distributor() const {
    using mesh::Children;

    Distributor distr;

    distr.split = [](AmrStorage::Item &parent, Children &children) {
        for (auto &child: children) {
            Vector3d dr = parent.center - child.center;
            child(U).u1 = parent(U).u1;
        }
    };

    distr.merge = [](Children &children, AmrStorage::Item& parent) {
        double sum = 0.0;
        for (auto &child: children) {
            sum += child(U).u1 * child.volume();
        }
        parent(U).u1 = sum / parent.volume();
    };

    return distr;
}

AmrStorage Transfer::body(EuMesh& mesh) {
    double eps = 1.0e-5;
    int count = 0;
    for (auto cell: mesh) {
        if (cell(U).u1 < eps) {
            continue;
        }
        ++count;
    }

    AmrStorage cells(U, count);

    count = 0;
    for (auto cell: mesh) {
        if (cell(U).u1 < eps) {
            continue;
        }

        auto poly = cell.geom().polygon();
        auto part = poly.clip(cell(U).p, cell(U).n);
        cells[count] = AmrCell(part);
        ++count;
    }

    return cells;
}

} // namespace zephyr::math