#include <zephyr/math/solver/transfer.h>

#include <zephyr/geom/polygon.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/geom/primitives/bface.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace mesh;

static const Transfer::State U = Transfer::datatype();

Transfer::State Transfer::datatype() {
    return {};
}

Transfer::Transfer() {
    m_dt  = 1.0e+300;
    m_CFL = 0.5;
    m_ver = 1;
    m_dir = Direction::ANY;
}

double Transfer::CFL() const {
    return m_CFL;
}

void Transfer::set_CFL(double C) {
    m_CFL = std::max(0.0, std::min(C, 1.0));
}

void Transfer::set_version(int ver) {
    m_ver = ver;
}

void Transfer::dir_splitting(bool flag) {
    if (flag) {
        m_dir = Direction::X;
    }
    else {
        m_dir = Direction::ANY;
    }
}

double Transfer::dt() const {
    switch (m_dir) {
        case Direction::X:
            return 0.0001 * m_dt;
        case Direction::Y:
            return 0.9999 * m_dt;
        default:
            return m_dt;
    }
}

Vector3d Transfer::velocity(const Vector3d &c) const {
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

double Transfer::compute_dt(EuMesh &mesh) {
    double dt = std::numeric_limits<double>::max();
    for (auto &cell: mesh) {
        dt = std::min(dt, compute_dt(cell));
    }
    return dt;
}

// Функция знака
inline double sign(double x) {
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

// Помещает x внутрь отрезка [x_min, x_max]
inline double between(double x_min, double x, double x_max) {
    return std::max(x_min, std::min(x, x_max));
}

inline std::tuple<double, double> minmax(double a1, double a2) {
    return std::make_tuple(std::min(a1, a2), std::max(a1, a2));
}

// a1, a2 -- объемные доли
// S -- площадь грани
// V1, V2 -- объемы ячеек
// as -- разбиение грани
// vn -- скорость
// dt -- шаг интегрирования
double flux_2D(double a1, double a2, double S, double V1, double V2, double as, double vn, double dt) {
    double a = vn > 0.0 ? a1 : a2;
    double V = vn > 0.0 ? V1 : V2;

    double gamma = dt * std::abs(vn) * S;

    double F_min = std::max(0.0, gamma - (1.0 - a) * V);
    double F_max = a * V;
    return sign(vn) * between(F_min, gamma * as, F_max);
}

double face_fraction(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);

    // Формула Серёжкина
    double a_sig = a_min / (1.0 - (a_max - a_min));

    // Случай a_min = 0, a_max = 1
    if (std::isnan(a_sig)) {
        a_sig = 0.5;
    }

    return between(a_min, a_sig, a_max);
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

        double S = face.area();
        double V1 = cell.volume();
        double V2 = neib.volume();

        fluxes += flux_2D(a1, a2, S, V1, V2, as, vs, m_dt);
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

// V1, V2 -- скорость в узлах грани
// fn -- нормаль к грани
// Предполагаем (V1 + V2).dot(fn) > 0.0
double flux_VOF(
        double a, Vector3d &p, Vector3d &n,
        EuCell &cell, EuFace &face,
        const Vector3d& V1, const Vector3d& V2,
        double dt, const Vector3d& fn) {

    // Нормальная скорость к грани
    Vector3d Vc = 0.5 * (V1 + V2);
    Vector3d Vn = Vc.dot(fn) * fn;
    Vector3d Vt = Vc - Vn;

    const auto &v1 = face.vs(0);
    const auto &v2 = face.vs(1);

    Line line1 = {v1, v1 - V1};
    Line line2 = {v2, v2 - V2};

    double xi = dt * Vn.norm();

    auto poly1 = cell.polygon();
    auto poly2 = poly1.clip(face.center() - xi * fn, -fn);
    auto poly3 = poly2.clip(line1.center(), line1.normal(face.center()));
    auto poly4 = poly3.clip(line2.center(), line2.normal(face.center()));

    auto& poly = poly4;

    if (a < 1.0e-6) {
        // Маленькую часть отправляем по нормали
        if (Vn.squaredNorm() > Vt.squaredNorm())
            return a * poly.area();
        else {
            return 0.0;
        }
    } else if (a > 1.0 - 1.0e-6) {
        // От полной ячейки отрезаем весь кусок
        return poly.area();
    }
    else {
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

        Vector3d fn = face.normal();
        Vector3d V1 = velocity(face.vs(0));
        Vector3d V2 = velocity(face.vs(1));

        // Расчет с расщеплением
        if (m_dir != Direction::ANY) {
            V1 = V1.dot(fn) * fn;
            V2 = V2.dot(fn) * fn;
        }

        // Типа upwind
        double F;
        if ((V1 + V2).dot(face.normal()) > 0.0) {
            F = +flux_VOF(zc.u1, zc.p, zc.n, cell, face, V1, V2, m_dt, fn);
        }
        else {
            auto zn = neib(U);
            F = -flux_VOF(zn.u1, zn.p, zn.n, neib, face, V1, V2, m_dt, -fn);
        }

        fluxes += F;
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

// Находит оптимальное деление грани a_sig, при котором поток максимально близок к F_VOF
double best_face_fraction(double a1, double a2, double S, double vs, double dt, double F_VOF) {
    // Да, удивительно, но вот так.
    if (vs == 0.0) {
        return 0.5 * (a1 + a2);
    }

    auto[a_min, a_max] = minmax(a1, a2);

    // F_VOF и vs имеют один знак
    assert(F_VOF * vs >= 0.0);

    return between(a_min, F_VOF / (dt * vs * S), a_max);
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

        const auto &fn = face.normal();
        Vector3d V1 = velocity(face.vs(0));
        Vector3d V2 = velocity(face.vs(1));

        // Расчет с расщеплением
        if (m_dir != Direction::ANY) {
            V1 = V1.dot(fn) * fn;
            V2 = V2.dot(fn) * fn;
        }

        double vs = 0.5 * (V1 + V2).dot(fn);

        double F_VOF;
        if (vs > 0.0) {
            F_VOF = +flux_VOF(zc.u1, zc.p, zc.n, cell, face, V1, V2, m_dt, fn);
        } else {
            F_VOF = -flux_VOF(zn.u1, zn.p, zn.n, neib, face, V1, V2, m_dt, -fn);
        }

        double a1 = zc.u1;
        double a2 = zn.u1;

        double S = face.area();
        double vol1 = cell.volume();
        double vol2 = neib.volume();

        // Хочу найти as, при котором flux_2D дает F_VOF
        //double a_sig = face_fraction(a1, a2);
        double a_sig = best_face_fraction(a1, a2, S, vs, m_dt, F_VOF);

        double F_CRP = flux_2D(a1, a2, S, vol1, vol2, a_sig, vs, m_dt);

        fluxes += F_CRP;
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

void Transfer::update(EuMesh &mesh) {
    switch (m_ver) {
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

void Transfer::update_dir() {
    // Флаг Direction::ANY не трогаем,
    // X и Y меняем местами

    if (m_dir == Direction::X) {
        m_dir = Direction::Y;
    }
    else if (m_dir == Direction::Y) {
        m_dir = Direction::X;
    }
}

void Transfer::update_ver1(EuMesh& mesh) {
    m_dt = compute_dt(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_CRP(cell, m_dir);
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }

    update_dir();
}

void Transfer::update_ver2(EuMesh& mesh) {
    // Определяем dt
    m_dt = compute_dt(mesh);

    compute_normals(mesh, 3);
    find_sections(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_VOF(cell, m_dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }

    compute_normals(mesh, 3);
    find_sections(mesh);

    update_dir();
}

void Transfer::update_ver3(EuMesh& mesh) {
    // Определяем dt
    m_dt = compute_dt(mesh);

    compute_normals(mesh, 3);
    find_sections(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_MIX(cell, m_dir);
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = std::max(0.0, std::min(cell(U).u2, 1.0));
        cell(U).u2 = 0.0;
    }

    compute_normals(mesh, 3);
    find_sections(mesh);

    update_dir();
}

// Производная по теореме Гаусса--Остроградского
void Transfer::compute_normal(EuCell& cell) {
    double uc = cell(U).u1;

    if (uc < 1.0e-6 || uc > 1.0 - 1.0e-6) {
        cell(U).n = Vector3d::Zero();
        return;
    }

    Vector3d n = Vector3d::Zero();
    for (auto &face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        double un = face.neib(U).u1;

        Vector3d S = 0.5 * face.normal() * face.area();
        n -= (uc + un) * S;
    }

    cell(U).n = n.normalized();
}

void Transfer::smooth_normal(EuCell& cell) {
    Vector3d n2 = Vector3d::Zero();

    if (cell(U).n.squaredNorm() == 0.0) {
        cell(U).n2 = n2;
        return;
    }

    int counter = 0;
    for (auto &face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        Vector3d nn = face.neib(U).n;
        if (nn.squaredNorm() > 0.0) {
            n2 += nn;
            ++counter;
        }
    }

    if (counter == 0) {
        cell(U).n2 = cell(U).n;
        int counter = 0;
        for (auto &face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            Vector3d nn = face.neib(U).n;
            if (nn.squaredNorm() > 0.0) {
                n2 += nn;
                ++counter;
            }
        }
        n2 += counter * cell(U).n;

    }
    else {
        n2 += counter * cell(U).n;
        cell(U).n2 = n2.normalized();
    }

}

void Transfer::update_normal(EuCell& cell) {
    cell(U).n = cell(U).n2;
}

void Transfer::find_section(EuCell &cell) {
    if (cell(U).zero_normal()) {
        cell(U).p = cell.center();
    } else {
        auto poly = cell.polygon();
        cell(U).p = poly.find_section(cell(U).n, cell(U).u1);
    }
}

void Transfer::compute_normals(EuMesh& mesh, int smoothing) {
    for (auto cell: mesh) {
        compute_normal(cell);
    }

    for (int i = 0; i < smoothing; ++i) {
        for (auto &cell: mesh) {
            smooth_normal(cell);
        }
        for (auto &cell: mesh) {
            update_normal(cell);
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
    int count = 0;
    for (auto cell: mesh) {
        if (cell(U).u1 < 1.0e-6) {
            continue;
        }
        ++count;
    }

    AmrStorage cells(U, count);

    count = 0;
    for (auto cell: mesh) {
        if (cell(U).u1 < 1.0e-6) {
            continue;
        }

        if (cell(U).u1 < 1.0 - 1.0e-6) {
            auto poly = cell.polygon();
            auto part = poly.clip(cell(U).p, cell(U).n);
            cells[count] = AmrCell(part);
        }
        else {
            cells[count] = cell.geom();
        }
        ++count;
    }

    return cells;
}

AmrStorage Transfer::scheme(EuMesh& mesh) {
    int count = 0;

    for (auto cell: mesh) {
        if (cell(U).u1 < 1.0e-12) {
            continue;
        }

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            const auto& V1 = velocity(face.vs(0));
            const auto& V2 = velocity(face.vs(1));

            // Нормаль к грани
            auto &fn = face.normal();

            if ((V1 + V2).dot(fn) < 0.0) {
                continue;
            }

            ++count;
        }
    }

    AmrStorage cells(U, count);

    count = 0;
    for (auto cell: mesh) {
        if (cell(U).u1 < 1.0e-12) {
            continue;
        }

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            const auto& V1 = velocity(face.vs(0));
            const auto& V2 = velocity(face.vs(1));

            // Нормаль к грани
            auto &fn = face.normal();

            if ((V1 + V2).dot(fn) < 0.0) {
                continue;
            }

            const auto& v1 = face.vs(0);
            const auto& v2 = face.vs(1);

            Line seg1 = {v1, v1 - V1};
            Line seg2 = {v2, v2 - V2};

            auto poly1 = cell.polygon();
            auto poly2 = poly1.clip(face.center() - 0.5 * m_dt * (V1 + V2).dot(fn) * fn, - fn);
            auto poly3 = poly2.clip(seg1.center(), seg1.normal(face.center()));
            auto poly4 = poly3.clip(seg2.center(), seg2.normal(face.center()));

            // Типа upwind
            double xi1 = V2.dot(fn) / V1.dot(fn);
            double xi2 = V1.dot(fn) / V2.dot(fn);

            if (V1.dot(fn) == 0.0) {
                if (V2.dot(fn) == 0.0) {
                    return 0.0;
                }

                xi1 = -1.0;
            }
            if (V2.dot(fn) == 0.0) {
                xi2 = -1.0;
            }

            Vector3d v3 = v1 - 0.5 * m_dt * (1.0 + xi1) * V1;
            Vector3d v4 = v2 - 0.5 * m_dt * (1.0 + xi2) * V2;

            //PolyQuad poly(v1, v2, v4, v3);
            //std::cout << poly5.size() << "\n";
            cells[count] = AmrCell(poly4);

            ++count;
        }
    }

    return cells;
}

} // namespace zephyr::math