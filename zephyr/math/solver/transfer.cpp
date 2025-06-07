#include <zephyr/math/solver/transfer.h>
#include <zephyr/mesh/primitives/bface.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/intersection.h>
#include <zephyr/geom/sections.h>

#include <zephyr/math/funcs.h>
#include <zephyr/math/calc/weno.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>

#include <zephyr/geom/sections.h>

namespace zephyr::math {

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
           m == Transfer::Method::MUSCLn_CRP ||
           m == Transfer::Method::MUSCL_MC ||
           m == Transfer::Method::MUSCL_MC_CRP;
}

inline bool WENO_type(Transfer::Method m) {
    return m == Transfer::Method::WENO ||
           m == Transfer::Method::WENO_CRP;
}

Transfer::Transfer() {
    m_dt  = 1.0e+300;
    m_CFL = 0.5;
    m_method = Method::VOF;
    m_limiter = "MC";
}

void Transfer::add_types(SoaMesh& mesh) {
    data.u1 = mesh.add<double>("u1");
    data.u2 = mesh.add<double>("u2");
    data.n = mesh.add<Vector3d>("n");
    data.p = mesh.add<Vector3d>("p");
    data.du_dx = mesh.add<double>("du/dx");
    data.du_dy = mesh.add<double>("du/dy");

    interface = InterfaceRecovery(data.u1, data.n, data.p);
}

double Transfer::CFL() const {
    return m_CFL;
}

void Transfer::set_CFL(double C) {
    m_CFL = std::max(0.0, std::min(C, 1.0));
}

Transfer::Method Transfer::method() const {
    return m_method;
}

void Transfer::set_method(Transfer::Method method) {
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

double Transfer::compute_dt(QCell &cell) {
    double max_area = 0.0;
    for (auto &face: cell.faces()) {
        max_area = std::max(max_area, face.area());
    }
    double dx = cell.volume() / max_area;
    return dx / velocity(cell.center()).norm();
}

double Transfer::compute_dt(SoaMesh &mesh) {
    double tau = std::numeric_limits<double>::max();

    for (auto &cell: mesh) {
        tau = std::min(tau, compute_dt(cell));
    }
    return m_CFL * tau;
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
    return sign(vn) * between(gamma * as, F_min, F_max);
}

// n1, n2 --- нормали к интерфейсу
// fn --- нормаль к грани
// vn --- нормальная компонента скорости
double face_fraction_n1(
        double a1, double a2, const Vector3d& n1, const Vector3d& n2,
        const Vector3d& fn, double vn, double V1, double V2, double S, double dt) {
    double a = vn > 0.0 ? a1 : a2;
    double V = vn > 0.0 ? V1 : V2;
    double cos = vn > 0.0 ? n1.dot(fn) : -n2.dot(fn);

    double C = dt * std::abs(vn) * S / V;

    double a_sig = average_flux(a, cos, C);
    return between(a_sig, a1, a2);
}

double face_fraction_n2(QCell& cell, QCell& neib, QFace& face, double vn, const Transfer::State& data) {
    // Отрезок - грань
    obj::segment seg{
            .v1 = face.vs(0),
            .v2 = face.vs(1)
    };

    // Реконструкция в ячейке
    obj::plane plane{
            .p = vn > 0.0 ? cell(data.p) : neib(data.p),
            .n = vn > 0.0 ? cell(data.n) : neib(data.n)
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

    auto [a_min, a_max] = minmax(cell(data.u1), neib(data.u1));

    return between(a_sig, a_min, a_max);
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

    return between(Flux / (dt * vn * S), a_min, a_max);
}

double flux_CRP(QCell& cell, QCell& neib, QFace& face, double vn, double dt, double Flux, const Transfer::State& data) {
    double a1 = cell(data.u1);
    double a2 = neib(data.u1);

    double S = face.area();
    double vol1 = cell.volume();
    double vol2 = neib.volume();

    // Хочу найти a_sig, при котором flux_2D дает Flux
    double a_sig = best_face_fraction(a1, a2, S, vn, dt, Flux);

    return flux_2D(a1, a2, S, vol1, vol2, a_sig, vn, dt);
}

void Transfer::fluxes_CRP(QCell &cell, Direction dir) {
    double a1 = cell(data.u1);
    Vector3d n1 = cell(data.n);

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        double a2 = neib(data.u1);
        Vector3d n2 = neib(data.n);

        Vector3d fn = face.normal();
        double vn = velocity(face.center()).dot(fn);

        double S = face.area();
        double V1 = cell.volume();
        double V2 = neib.volume();

        double a_sig;
        switch (m_method) {
            case Method::CRP_V3:
                a_sig = face_fraction_v3(a1, a2);
                break;
            case Method::CRP_V5:
                a_sig = face_fraction_v5(a1, a2);
                break;
            case Method::CRP_N1:
                a_sig = face_fraction_n1(a1, a2, n1, n2, fn, vn, V1, V2, S, m_dt);
                break;
            case Method::CRP_N2:
                a_sig = face_fraction_n2(cell, neib, face, vn, data);
                break;
            default:
                a_sig = face_fraction_s(a1, a2);
                break;
        }

        fluxes += flux_2D(a1, a2, S, V1, V2, a_sig, vn, m_dt);
    }

    cell(data.u2) = cell(data.u1) - fluxes / cell.volume();
}

// V1, V2 -- скорость в узлах грани
// fn -- нормаль к грани
// Предполагаем (V1 + V2).dot(fn) > 0.0
double flux_VOF(QCell &cell, QFace &face,
        const Vector3d& V1, const Vector3d& V2,
        double dt, const Vector3d& fn, const Transfer::State& data) {

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

    if (cell(data.u1) < 1.0e-8) {
        // Маленькую часть отправляем по нормали
        if (Vn.squaredNorm() > Vt.squaredNorm())
            return cell(data.u1) * poly.area();
        else {
            return 0.0;
        }
    } else if (cell(data.u1) > 1.0 - 1.0e-8) {
        // От полной ячейки отрезаем весь кусок
        return poly.area();
    }
    else {
        return poly.clip_area(cell(data.p), cell(data.n));
    }
}

void Transfer::fluxes_VOF(QCell &cell, Direction dir) {
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
            // Для четырехугольных ячеек
            // double C = std::abs(vn) * m_dt * face.area() / cell.volume();
            // Flux = vn * m_dt * face.area() * a_sigma_vof(zc.u1, fn.dot(zc.n), C);

            Flux = +flux_VOF(cell, face, V1, V2, m_dt, fn, data);
        }
        else {
            // Для четырехугольных ячеек
            // double C = std::abs(vn) * m_dt * face.area() / neib.volume();
            // Flux = vn * m_dt * face.area() * a_sigma_vof(zn.u1, -fn.dot(zn.n), C);

            Flux = -flux_VOF(neib, face, V1, V2, m_dt, -fn, data);
        }

        // CRP поправка
        if (m_method == Method::VOF_CRP) {
            Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
        }

        fluxes += Flux;
    }

    cell(data.u2) = cell(data.u1) - fluxes / cell.volume();
}

void Transfer::fluxes_MUSCL(QCell &cell, Direction dir) {
    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();

        auto fn = face.normal();
        double vn = velocity(face.center()).dot(fn);

        double a_sig = NAN;
        if (m_method == Method::MUSCL_MC || m_method == Method::MUSCL_MC_CRP)
        {
            auto fe = FaceExtra::Direct(
                    cell(data.u1), cell(data.du_dx), cell(data.du_dy), 0.0,
                    neib(data.u1), neib(data.du_dx), neib(data.du_dy), 0.0,
                    cell.center(), neib.center(), face.center());

            a_sig = vn > 0.0 ? fe.m(cell(data.u1)) : fe.p(neib(data.u1));
            a_sig = between(a_sig, cell(data.u1), neib(data.u1));
        }
        else {
            auto fe = FaceExtra::ATvL(
                    cell(data.u1), cell(data.du_dx), cell(data.du_dy), 0.0,
                    neib(data.u1), neib(data.du_dx), neib(data.du_dy), 0.0,
                    cell.center(), neib.center(), face.center());

            a_sig = vn > 0.0 ? fe.m(cell(data.u1)) : fe.p(neib(data.u1));
        }

        double Flux = a_sig * vn * m_dt * face.area();

        // CRP поправка
        if (m_method == Method::MUSCLd_CRP ||
            m_method == Method::MUSCLn_CRP ||
            m_method == Method::MUSCL_MC_CRP) {
            Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
        }

        fluxes += Flux;
    }

    cell(data.u2) = cell(data.u1) - fluxes / cell.volume();
}

void Transfer::compute_slopes(SoaMesh& mesh) {
    if (m_method == Method::MUSCLn ||
        m_method == Method::MUSCLn_CRP) {
        for (auto cell: mesh) {
            double grad_x = 0.0;
            double grad_y = 0.0;

            // Реконструкция в ячейке
            obj::plane plane{
                    .p = cell(data.p),
                    .n = cell(data.n)
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

                grad_x += a_sig * face.area() * face.normal().x();
                grad_y += a_sig * face.area() * face.normal().y();
            }
            cell(data.du_dx) = grad_x / cell.volume();
            cell(data.du_dy) = grad_y / cell.volume();
        }

        return;
    }

    auto u1 = data.u1;
    auto get_state = [u1](QCell& cell) -> double {
        return cell(u1);
    };
    auto boundary_value = [](double u, const Vector3d& n, Boundary b) -> double {
        return u;
    };


    for (auto cell: mesh) {
        auto grad = gradient::LSM<double>(cell, get_state, boundary_value);
        cell(data.du_dx) = grad.x;
        cell(data.du_dy) = grad.y;

        if (m_method == Method::MUSCL_MC || m_method == Method::MUSCL_MC_CRP) {
            auto lim_grad = gradient::limiting<double>(cell, m_limiter,
                    grad, get_state, boundary_value);

            cell(data.du_dx) = lim_grad.x;
            cell(data.du_dy) = lim_grad.y;
        }
    }
}

void Transfer::update(SoaMesh &mesh, Direction dir) {
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
        throw std::runtime_error("Unknown solver method");
    }
}

void Transfer::update_CRP(SoaMesh& mesh, Direction dir) {
    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_CRP(cell, dir);
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(data.u1) =  between(cell(data.u2), 0.0, 1.0);
        cell(data.u2) = 0.0;
    }

    // Без сглаживаний, чисто для реконструкции
    update_interface(mesh, 0);
}

void Transfer::update_VOF(SoaMesh& mesh, Direction dir) {
    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_VOF(cell, dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell(data.u1) =  between(cell(data.u2), 0.0, 1.0);
        cell(data.u2) = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_MUSCL(SoaMesh& mesh, Direction dir) {
    compute_slopes(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_MUSCL(cell, dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell(data.u1) = between(cell(data.u2), 0.0, 1.0);
        cell(data.u2) = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_WENO(SoaMesh& mesh, Direction dir) {
    // Считаем потоки
    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);

            double fluxes = 0.0;

            // LEFT
            if (dir == Direction::X || dir == Direction::ANY) {
                auto face = cell.face(Side3D::L);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int I = vn > 0.0 ? i : i - 1;

                WENO5 weno {
                    mesh(I - 2, j)(data.u1),
                    mesh(I - 1, j)(data.u1),
                    mesh(I + 0, j)(data.u1),
                    mesh(I + 1, j)(data.u1),
                    mesh(I + 2, j)(data.u1),
                };
                double a_sig = vn > 0.0 ? weno.m() : weno.p();
                a_sig =  between(a_sig, 0.0, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
                }
                fluxes += Flux;
            }

            // RIGHT
            if (dir == Direction::X || dir == Direction::ANY) {
                auto face = cell.face(Side3D::R);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int I = vn > 0.0 ? i : i + 1;

                WENO5 weno {
                        mesh(I - 2, j)(data.u1),
                        mesh(I - 1, j)(data.u1),
                        mesh(I + 0, j)(data.u1),
                        mesh(I + 1, j)(data.u1),
                        mesh(I + 2, j)(data.u1),
                };

                double a_sig = vn > 0.0 ? weno.p() : weno.m();
                a_sig =  between(a_sig, 0.0, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
                }
                fluxes += Flux;
            }

            // BOTTOM
            if (dir == Direction::Y || dir == Direction::ANY) {
                auto face = cell.face(Side3D::B);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int J = vn > 0.0 ? j : j - 1;

                WENO5 weno {
                        mesh(i, J - 2)(data.u1),
                        mesh(i, J - 1)(data.u1),
                        mesh(i, J + 0)(data.u1),
                        mesh(i, J + 1)(data.u1),
                        mesh(i, J + 2)(data.u1),
                };

                double a_sig = vn > 0.0 ? weno.m() : weno.p();
                a_sig =  between(a_sig, 0.0, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
                }
                fluxes += Flux;
            }

            // TOP
            if (dir == Direction::Y || dir == Direction::ANY) {
                auto face = cell.face(Side3D::T);
                auto neib = face.neib();
                double vn = velocity(face.center()).dot(face.normal());

                int J = vn > 0.0 ? j : j + 1;

                WENO5 weno {
                        mesh(i, J - 2)(data.u1),
                        mesh(i, J - 1)(data.u1),
                        mesh(i, J + 0)(data.u1),
                        mesh(i, J + 1)(data.u1),
                        mesh(i, J + 2)(data.u1),
                };
                double a_sig = vn > 0.0 ? weno.p() : weno.m();
                a_sig =  between(a_sig, 0.0, 1.0);
                double Flux = a_sig * vn * m_dt * face.area();

                // CRP поправка
                if (m_method == Method::WENO_CRP) {
                    Flux = flux_CRP(cell, neib, face, vn, m_dt, Flux, data);
                }
                fluxes += Flux;
            }

            cell(data.u2) = cell(data.u1) - fluxes / cell.volume();
        }
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell(data.u1) = between(cell(data.u2), 0.0, 1.0);
        cell(data.u2) = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_interface(SoaMesh& mesh, int smoothing) {
    interface.update(mesh, smoothing);
}

void Transfer::set_flags(SoaMesh& mesh) {
    for (auto cell: mesh) {
        double min_val = cell(data.u1);
        double max_val = cell(data.u1);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            min_val = std::min(min_val, face.neib(data.u1));
            max_val = std::max(max_val, face.neib(data.u1));
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
    Distributor distr;

    auto u1 = data.u1;

    distr.split_soa = [u1](QCell& parent, SoaChildren &children) {
        for (auto child: children) {
            Vector3d dr = parent.center() - child.center();
            child(u1) = parent(u1);
        }
    };

    distr.merge_soa = [u1](SoaChildren &children, QCell& parent) {
        double sum = 0.0;
        for (auto child: children) {
            sum += child(u1) * child.volume();
        }
        parent(u1) = sum / parent.volume();
    };

    return distr;
}

SoaMesh Transfer::body(SoaMesh& mesh) {
    return interface.body(mesh);
}

} // namespace zephyr::math