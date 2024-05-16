#include <zephyr/math/solver/transfer.h>

#include <zephyr/geom/primitives/bface.h>
#include <zephyr/geom/polygon.h>
#include <zephyr/geom/intersection.h>
#include <zephyr/math/cfd/face_extra.h>

#include <zephyr/math/calc/weno.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;
using namespace mesh;

inline bool CRP_type(Transfer::Method m) {
    return m == Transfer::Method::CRP_V3 ||
           m == Transfer::Method::CRP_V5 ||
           m == Transfer::Method::CRP_SE ||
           m == Transfer::Method::CRP_N1 ||
           m == Transfer::Method::CRP_N2;
}

inline bool VOF_type(Transfer::Method m) {
    return m == Transfer::Method::VOF ||
           m == Transfer::Method::VOF_CRP;
}

inline bool MUSCL_type(Transfer::Method m) {
    return m == Transfer::Method::MUSCLd ||
           m == Transfer::Method::MUSCLn ||
           m == Transfer::Method::MUSCLd_CRP ||
           m == Transfer::Method::MUSCLn_CRP;
}

inline bool WENO_type(Transfer::Method m) {
    return m == Transfer::Method::WENO ||
           m == Transfer::Method::WENO_CRP;
}

static const Transfer::State U = Transfer::datatype();

Transfer::State Transfer::datatype() {
    return {};
}

Transfer::Transfer()
    : interface{offsetof(State, u1),
                offsetof(State, n),
                offsetof(State, p)} {

    m_dt  = 1.0e+300;
    m_CFL = 0.5;
    m_method = Method::VOF;
}

double Transfer::CFL() const {
    return m_CFL;
}

void Transfer::set_CFL(double C) {
    m_CFL = std::max(0.0, std::min(C, 1.0));
}

void Transfer::Transfer::set_method(Transfer::Method method) {
    m_method = method;
}

double Transfer::get_dt() const {
    return m_dt;
}

void Transfer::set_dt(double dt) {
    m_dt = dt;
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

//inline double sqr(double x) { return x * x; }

// Функция знака
inline double sign(double x) {
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

// Функция Хевисайда
inline double heav(double x) {
    return x > 0.0 ? 1.0 : 0.0;
}

// гладкая функция знака
inline double sign_s(double x) {
    const double x0 = 0.1;
    double xi = std::abs(x / x0);
    return sign(x) * (1.0 - heav(1.0 - xi) * sqr(1.0 - xi));
}

// Гладкая функция Хевисайда
inline double heav_s(double x) {
    const double x0 = 0.05;
    double xi = std::abs(x / x0);
    return heav(x) * (1.0 - heav(1.0 - xi) * sqr(1.0 - xi));
}

// Помещает x внутрь отрезка [x_min, x_max]
inline double between(double x_min, double x, double x_max) {
    return std::max(x_min, std::min(x, x_max));
}

inline std::tuple<double, double> minmax(double a1, double a2) {
    return std::make_tuple(std::min(a1, a2), std::max(a1, a2));
}

inline bool is_zero(const Vector3d& p) {
    return p.x() == 0.0 && p.y() == 0.0 && p.z() == 0.0;
}

inline double cross(const Vector3d& v1, const Vector3d& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
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

double face_fraction_v3(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);
    double a_avg = 0.5 * (a1 + a2);
    double a_s = a_avg + 0.5 * sign_s(2 * a_avg - 1) * (a_max - a_min);
    return between(a_min, a_s, a_max);
}

double face_fraction_v5(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);
    double a_avg = 0.5 * (a1 + a2);
    double delta = 0.5 * std::abs(a2 - a1);

    double L = 0.5 * (1.0 - sign_s(2.0 * a_avg - 1.0));
    double G = 0.5 * (1.0 + sign_s(2.0 * a_avg - 1.0));

    double a_s = L * a_avg * heav_s(a_avg - delta) +
                 G * (1.0 - (1.0 - a_avg) * heav_s(1.0 - a_avg - delta));

    return between(a_min, a_s, a_max);
}

// Формула Серёжкина
double face_fraction_s(double a1, double a2) {
    auto [a_min, a_max] = minmax(a1, a2);

    double a_sig = a_min / (1.0 - (a_max - a_min));

    // Случай a_min = 0, a_max = 1
    if (std::isnan(a_sig)) {
        a_sig = 0.5;
    }

    return between(a_min, a_sig, a_max);
}

// n1, n2 --- нормали к интерфейсу
// fn --- нормаль к грани
// vn --- нормальная компонента скорости
double face_fraction_n1(double a1, double a2, const Vector3d& n1, const Vector3d& n2,
                        const Vector3d& fn, double vn) {
    auto [a_min, a_max] = minmax(a1, a2);

    double a_ser = face_fraction_s(a1, a2);
    double a_up = vn > 0.0 ? a1 : a2;

    double cos = fn.dot(vn > 0.0 ? n1 : n2);
    double xi = std::abs(cos);

    double a_sig = xi * a_ser + (1.0 - xi) * a_up;

    return between(a_min, a_sig, a_max);
}

double face_fraction_n2(EuCell& cell, EuCell& neib, EuFace& face, double vn) {
    // Отрезок - грань
    obj::segment seg{
            .v1 = face.vs(0),
            .v2 = face.vs(1)
    };

    // Реконструкция в ячейке
    obj::plane plane{
            .p = vn > 0.0 ? cell(U).p : neib(U).p,
            .n = vn > 0.0 ? cell(U).n : neib(U).n
    };

    bool in1 = plane.under(seg.v1);
    bool in2 = plane.under(seg.v2);

    double a_sig;
    if (in1 && in2) {
        a_sig = 1.0;
    }
    else if (!in1 && !in2) {
        a_sig = 0.0;
    }
    else {
        Vector3d in = geom::intersection2D::find_fast(plane, seg);

        if (in1) {
            a_sig = (in - seg.v1).norm() / seg.length();
        }
        else {
            a_sig = (in - seg.v2).norm() / seg.length();
        }
    }

    auto [a_min, a_max] = minmax(cell(U).u1, neib(U).u1);

    return between(a_min, a_sig, a_max);
}

// Находит оптимальное деление грани a_sig, при котором поток максимально близок к Flux
double best_face_fraction(double a1, double a2, double S, double vn, double dt, double Flux) {
    // Да, удивительно, но вот так.
    if (vn == 0.0) {
        return 0.5 * (a1 + a2);
    }

    auto[a_min, a_max] = minmax(a1, a2);

    // Flux и vn имеют один знак
    assert(Flux * vn >= 0.0);

    return between(a_min, Flux / (dt * vn * S), a_max);
}

double flux_CRP(EuCell& cell, EuCell& neib, EuFace& face, double vn, double dt, double Flux) {
    double a1 = cell(U).u1;
    double a2 = neib(U).u1;

    double S = face.area();
    double vol1 = cell.volume();
    double vol2 = neib.volume();

    // Хочу найти a_sig, при котором flux_2D дает Flux
    double a_sig = best_face_fraction(a1, a2, S, vn, dt, Flux);

    return flux_2D(a1, a2, S, vol1, vol2, a_sig, vn, dt);
}

void Transfer::fluxes_CRP(EuCell &cell, Direction dir) {
    auto &zc = cell(U);
    double a1 = zc.u1;
    Vector3d n1 = cell(U).n;

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        double a2 = neib(U).u1;
        Vector3d n2 = neib(U).n;

        Vector3d fn = face.normal();
        double vn = velocity(face.center()).dot(fn);

        double a_sig;
        switch (m_method) {
            case Method::CRP_V3:
                a_sig = face_fraction_v3(a1, a2);
                break;
            case Method::CRP_V5:
                a_sig = face_fraction_v5(a1, a2);
                break;
            case Method::CRP_N1:
                a_sig = face_fraction_n1(a1, a2, n1, n2, fn, vn);
                break;
            case Method::CRP_N2:
                a_sig = face_fraction_n2(cell, neib, face, vn);
                break;
            default:
                a_sig = face_fraction_s(a1, a2);
                break;
        }

        double S = face.area();
        double V1 = cell.volume();
        double V2 = neib.volume();

        fluxes += flux_2D(a1, a2, S, V1, V2, a_sig, vn, m_dt);
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

    if (a < 1.0e-8) {
        // Маленькую часть отправляем по нормали
        if (Vn.squaredNorm() > Vt.squaredNorm())
            return a * poly.area();
        else {
            return 0.0;
        }
    } else if (a > 1.0 - 1.0e-8) {
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
        double vn = 0.5 * (V1 + V2).dot(fn);

        // Расчет с расщеплением
        if (dir != Direction::ANY) {
            V1 = V1.dot(fn) * fn;
            V2 = V2.dot(fn) * fn;
        }

        // Типа upwind
        double Flux;
        if (vn > 0.0) {
            Flux = +flux_VOF(zc.u1, zc.p, zc.n, cell, face, V1, V2, m_dt, fn);
        }
        else {
            auto zn = neib(U);
            Flux = -flux_VOF(zn.u1, zn.p, zn.n, neib, face, V1, V2, m_dt, -fn);
        }

        // CRP поправка
        if (m_method == Method::VOF_CRP) {
            Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
        }

        fluxes += Flux;
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

void Transfer::fluxes_MUSCL(EuCell &cell, Direction dir) {
    auto &zc = cell(U);

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        auto &zn = neib(U);

        auto& fn = face.normal();
        double vn = velocity(face.center()).dot(fn);

        auto fe = FaceExtra::ATvL(
                zc.u1, zc.du_dx, zc.du_dy, 0.0,
                zn.u1, zn.du_dx, zn.du_dy, 0.0,
                cell.center(), neib.center(), face.center());

        double a_sig = vn > 0.0 ? fe.m(zc.u1) : fe.p(zn.u1);
        double Flux = a_sig * vn * m_dt * face.area();

        // CRP поправка
        if (m_method == Method::MUSCLd_CRP ||
            m_method == Method::MUSCLn_CRP) {
            Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
        }

        fluxes += Flux;
    }

    zc.u2 = zc.u1 - fluxes / cell.volume();
}

void Transfer::compute_slopes(EuMesh& mesh) {
    if (m_method == Method::MUSCLn ||
        m_method == Method::MUSCLn_CRP) {
        for (auto cell: mesh) {
            double du_dx = 0.0;
            double du_dy = 0.0;

            // Реконструкция в ячейке
            obj::plane plane{
                    .p = cell(U).p,
                    .n = cell(U).n
            };

            for (auto face: cell.faces()) {
                // Отрезок - грань
                obj::segment seg{
                        .v1 = face.vs(0),
                        .v2 = face.vs(1)
                };

                bool in1 = plane.under(seg.v1);
                bool in2 = plane.under(seg.v2);

                double a_sig;
                if (in1 && in2) {
                    a_sig = 1.0;
                } else if (!in1 && !in2) {
                    a_sig = 0.0;
                } else {
                    Vector3d in = geom::intersection2D::find_fast(plane, seg);

                    if (in1) {
                        a_sig = (in - seg.v1).norm() / seg.length();
                    } else {
                        a_sig = (in - seg.v2).norm() / seg.length();
                    }
                }

                du_dx += a_sig * face.area() * face.normal().x();
                du_dy += a_sig * face.area() * face.normal().y();
            }
            cell(U).du_dx = du_dx / cell.volume();
            cell(U).du_dy = du_dy / cell.volume();
        }

        return;
    }

    for (auto cell: mesh) {
        double du_dx = 0.0;
        double du_dy = 0.0;

        double uc = cell(U).u1;
        for (auto face: cell.faces()) {
            double un = face.neib(U).u1;

            du_dx += 0.5 * (uc + un) * face.area() * face.normal().x();
            du_dy += 0.5 * (uc + un) * face.area() * face.normal().y();
        }

        cell(U).du_dx = du_dx / cell.volume();
        cell(U).du_dy = du_dy / cell.volume();
    }
}

void Transfer::update(EuMesh &mesh, Direction dir) {
    if (CRP_type(m_method)) {
        update_CRP(mesh, dir);
    }
    else if (VOF_type(m_method)) {
        update_VOF(mesh, dir);
    }
    else if (MUSCL_type(m_method)) {
        update_MUSCL(mesh, dir);
    }
    else if (WENO_type(m_method)) {
        update_WENO(mesh, dir);
    }
    else {
        update_other(mesh, dir);
    }
}

void Transfer::update_CRP(EuMesh& mesh, Direction dir) {
    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_CRP(cell, dir);
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = between(0.0, cell(U).u2, 1.0);
        cell(U).u2 = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_VOF(EuMesh& mesh, Direction dir) {
    for (auto cell: mesh) {
        if (cell(U).u1 < 0.0) {
            throw std::runtime_error("WTF MAN 1?");
        }
    }
    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_VOF(cell, dir);
    }
    for (auto cell: mesh) {
        if (cell(U).u1 < 0.0) {
            throw std::runtime_error("WTF MAN 2?");
        }
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell(U).u1 = between(0.0, cell(U).u2, 1.0);
        cell(U).u2 = 0.0;
    }
    for (auto cell: mesh) {
        if (cell(U).u1 < 0.0) {
            throw std::runtime_error("WTF MAN 3?");
        }
    }

    update_interface(mesh);
    for (auto cell: mesh) {
        if (cell(U).u1 < 0.0) {
            throw std::runtime_error("WTF MAN 4?");
        }
    }
}

void Transfer::update_MUSCL(EuMesh& mesh, Direction dir) {
    compute_slopes(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_MUSCL(cell, dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell(U).u1 = between(0.0, cell(U).u2, 1.0);
        cell(U).u2 = 0.0;
    }

    update_interface(mesh);

}

void Transfer::update_WENO(EuMesh& mesh, Direction dir) {
    using zephyr::math::WENO5;

    // Считаем потоки
    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);

            double fluxes = 0.0;

            // LEFT
            if (dir == Direction::X || dir == Direction::ANY) {
                auto face = cell.face(Side::L);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int I = vn > 0.0 ? i : i - 1;

                WENO5 weno {
                    mesh(I - 2, j)(U).u1,
                    mesh(I - 1, j)(U).u1,
                    mesh(I + 0, j)(U).u1,
                    mesh(I + 1, j)(U).u1,
                    mesh(I + 2, j)(U).u1,
                };
                double a_sig = vn > 0.0 ? weno.m() : weno.p();
                a_sig = between(0.0, a_sig, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
                }
                fluxes += Flux;
            }

            // RIGHT
            if (dir == Direction::X || dir == Direction::ANY) {
                auto face = cell.face(Side::R);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int I = vn > 0.0 ? i : i + 1;

                WENO5 weno {
                        mesh(I - 2, j)(U).u1,
                        mesh(I - 1, j)(U).u1,
                        mesh(I + 0, j)(U).u1,
                        mesh(I + 1, j)(U).u1,
                        mesh(I + 2, j)(U).u1,
                };

                double a_sig = vn > 0.0 ? weno.p() : weno.m();
                a_sig = between(0.0, a_sig, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
                }
                fluxes += Flux;
            }

            // BOTTOM
            if (dir == Direction::Y || dir == Direction::ANY) {
                auto face = cell.face(Side::B);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int J = vn > 0.0 ? j : j - 1;

                WENO5 weno {
                        mesh(i, J - 2)(U).u1,
                        mesh(i, J - 1)(U).u1,
                        mesh(i, J + 0)(U).u1,
                        mesh(i, J + 1)(U).u1,
                        mesh(i, J + 2)(U).u1,
                };

                double a_sig = vn > 0.0 ? weno.m() : weno.p();
                a_sig = between(0.0, a_sig, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
                }
                fluxes += Flux;
            }

            // TOP
            if (dir == Direction::Y || dir == Direction::ANY) {
                auto face = cell.face(Side::T);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int J = vn > 0.0 ? j : j + 1;

                WENO5 weno {
                        mesh(i, J - 2)(U).u1,
                        mesh(i, J - 1)(U).u1,
                        mesh(i, J + 0)(U).u1,
                        mesh(i, J + 1)(U).u1,
                        mesh(i, J + 2)(U).u1,
                };
                double a_sig = vn > 0.0 ? weno.p() : weno.m();
                a_sig = between(0.0, a_sig, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux);
                }
                fluxes += Flux;
            }

            cell(U).u2 = cell(U).u1 - fluxes / cell.volume();
        }
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(U).u1 = between(0.0, cell(U).u2, 1.0);
        cell(U).u2 = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_other(EuMesh& mesh, Direction dir) {

}

void Transfer::update_interface(EuMesh& mesh, int smoothing) {
    interface.update(mesh, smoothing);
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
    return interface.body(mesh);
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