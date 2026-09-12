#include <zephyr/math/solver/transfer.h>
#include <zephyr/geom/intersection.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/geom.h>

#include <zephyr/math/funcs.h>
#include <zephyr/math/calc/weno.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>

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
    m_dim = 2;
    m_method = Method::VOF;
    m_limiter = "MC";
    m_max_dt = std::numeric_limits<double>::max();
}

void Transfer::set_dim(int dim) {
    m_dim = std::max(2, std::min(dim, 3));
}

Transfer::State Transfer::add_types(Mesh& mesh) {
    data.u1 = mesh.add<double>("u1");
    data.u2 = mesh.add<double>("u2");
    data.n = mesh.add<Vector3d>("n");
    data.p = mesh.add<double>("p");
    data.grad = mesh.add<Vector3d>("grad");
    return data;
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

void Transfer::set_plic_type(Plic::Type type) {
    if (!data.u1) {
        throw std::runtime_error("Call 'add_types' before set PLIC type");
    }

    m_plic = Plic(m_dim, true, type,
        [&u=data.u1](const Cell& cell, int idx) -> double {
            return cell[u];
        });
}

double Transfer::get_dt() const {
    return m_dt;
}

void Transfer::set_dt(double dt) {
    m_dt = std::min(dt, m_max_dt);
}

void Transfer::set_max_dt(double dt) {
    m_max_dt = dt;
}

Vector3d Transfer::velocity(const Vector3d &c) const {
    return Vector3d::UnitX();
}

double Transfer::compute_dt(Cell &cell) const {
    double max_area = 0.0;
    for (auto &face: cell.faces()) {
        max_area = std::max(max_area, face.area());
    }
    double dx = cell.volume() / max_area;
    return dx / velocity(cell.center()).norm();
}

double Transfer::compute_dt(Mesh &mesh) const {
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

// Точная доля отсечения от грани для заданной плоскости
double face_fraction_n1(Cell& cell, Cell& neib, Face& face, double vn, const Transfer::State& data) {
    // Реконструкция в ячейке
    obj::plane plane{
        .p = vn > 0.0 ? cell[data.p] + cell.center().dot(cell[data.n]) : neib[data.p] + neib.center().dot(neib[data.n]),
        .n = vn > 0.0 ? cell[data.n] : neib[data.n]
    };

    double a_sig;
    if (cell.dim() == 2) {
        // Грань это отрезок
        obj::segment seg{
            .v1 = face.vs(0),
            .v2 = face.vs(1)
        };
        a_sig = intersection2D::edge_fraction(seg, plane);
    }
    else {
        throw std::runtime_error("intersection3D::quad fraction");
    }
    auto [a_min, a_max] = sorted(cell[data.u1], neib[data.u1]);
    return between(a_sig, a_min, a_max);
}

// n1, n2 --- нормали к интерфейсу
// fn --- нормаль к грани
// vn --- нормальная компонента скорости
double face_fraction_n2(
        double a1, double a2, const Vector3d& n1, const Vector3d& n2,
        const Vector3d& fn, double vn, double V1, double V2, double S, double dt) {
    double a = vn > 0.0 ? a1 : a2;
    double V = vn > 0.0 ? V1 : V2;
    double cos = vn > 0.0 ? n1.dot(fn) : -n2.dot(fn);

    double C = dt * std::abs(vn) * S / V;

    double a_sig = average_flux(a, cos, C);
    return between(a_sig, a1, a2);
}

double face_fraction_n2_3D(
        double a1, double a2, const Vector3d& n1, const Vector3d& n2,
        const Vector3d& fn, double vn, double V1, double V2, double S, double dt) {
    double a = vn > 0.0 ? a1 : a2;
    double V = vn > 0.0 ? V1 : V2;
    Vector3d n = vn > 0.0 ? n1 : n2;

    double C = dt * std::abs(vn) * S / V;

    double a_sig = average_flux(a, n, vn > 0.0 ? fn : -fn, C);
    return between(a_sig, a1, a2);
}

// Находит оптимальное деление грани a_sig, при котором поток максимально близок к Flux
double best_face_fraction(double a1, double a2, double S, double vn, double dt, double Flux) {
    // Да, удивительно, но вот так.
    if (vn == 0.0) {
        return 0.5 * (a1 + a2);
    }

    auto[a_min, a_max] = sorted(a1, a2);

    // Flux и vn имеют один знак
    z_assert(Flux * vn >= 0.0, "Strange flux");

    return between(Flux / (dt * vn * S), a_min, a_max);
}

double flux_CRP(Cell& cell, Cell& neib, Face& face, double vn, double dt, double Flux, const Transfer::State& data) {
    double a1 = cell[data.u1];
    double a2 = neib[data.u1];

    double S = face.area();
    double vol1 = cell.volume();
    double vol2 = neib.volume();

    // Хочу найти a_sig, при котором flux_2D дает Flux
    double a_sig = best_face_fraction(a1, a2, S, vn, dt, Flux);

    return flux_2D(a1, a2, S, vol1, vol2, a_sig, vn, dt);
}

void Transfer::fluxes_CRP(Cell &cell, Direction dir) const {
    double a1 = cell[data.u1];
    Vector3d n1 = cell[data.n];

    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();
        double a2 = neib[data.u1];
        Vector3d n2 = neib[data.n];

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
                a_sig = face_fraction_n1(cell, neib, face, vn, data);
                break;
            case Method::CRP_N2:
                if (cell.dim() == 2) {
                    a_sig = face_fraction_n2(a1, a2, n1, n2, fn, vn, V1, V2, S, m_dt);
                }
                else {
                    a_sig = face_fraction_n2_3D(a1, a2, n1, n2, fn, vn, V1, V2, S, m_dt);
                }
                break;
            default:
                a_sig = face_fraction_s(a1, a2);
                break;
        }

        fluxes += flux_2D(a1, a2, S, V1, V2, a_sig, vn, m_dt);
    }

    cell[data.u2] = cell[data.u1] - fluxes / cell.volume();
}

// V1, V2 -- скорость в узлах грани
// fn -- нормаль к грани
// Предполагаем (V1 + V2).dot(fn) > 0.0
double flux_VOF(Cell &cell, Face &face,
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

    if (cell[data.u1] < 1.0e-8) {
        // Маленькую часть отправляем по нормали
        if (Vn.squaredNorm() > Vt.squaredNorm())
            return cell[data.u1] * poly.area();
        else {
            return 0.0;
        }
    } else if (cell[data.u1] > 1.0 - 1.0e-8) {
        // От полной ячейки отрезаем весь кусок
        return poly.area();
    }
    else {
        Vector3d P = cell.center() + cell[data.p] * cell[data.n];
        return poly.clip_area(P, cell[data.n]);
    }
}

void Transfer::fluxes_VOF(Cell &cell, Direction dir) const {
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

    cell[data.u2] = cell[data.u1] - fluxes / cell.volume();
}

void Transfer::fluxes_MUSCL(Cell &cell, Direction dir) const {
    double fluxes = 0.0;
    for (auto &face: cell.faces(dir)) {
        if (face.is_boundary()) {
            continue;
        }

        auto neib = face.neib();

        auto fn = face.normal();
        double vn = velocity(face.center()).dot(fn);

        bool interesting = (cell[data.u1] > 0.01 && cell[data.u1] < 0.99) || (neib[data.u1] > 0.01 && neib[data.u1] < 0.99);

        double a_sig = NAN;
        if (m_method == Method::MUSCL_MC || m_method == Method::MUSCL_MC_CRP)
        {
            if (interesting) {
                //std::cout << cell[data.u1] << " " <<  cell[data.du_dx] << " " << cell[data.du_dy] << "\n";
                //std::cout << neib[data.u1] << " " <<  neib[data.du_dx] << " " << neib[data.du_dy] << "\n";
            }
            auto fe = FaceExtra::Direct(
                    cell[data.u1], cell[data.grad].x(), cell[data.grad].y(), 0.0,
                    neib[data.u1], neib[data.grad].x(), neib[data.grad].y(), 0.0,
                    cell.center(), neib.center(), face.center());

            a_sig = vn > 0.0 ? fe.m(cell[data.u1]) : fe.p(neib[data.u1]);
            a_sig = between(a_sig, cell[data.u1], neib[data.u1]);
        }
        else {
            auto fe = FaceExtra::ATvL(
                    cell[data.u1], cell[data.grad].x(), cell[data.grad].y(), 0.0,
                    neib[data.u1], neib[data.grad].x(), neib[data.grad].y(), 0.0,
                    cell.center(), neib.center(), face.center());

            a_sig = vn > 0.0 ? fe.m(cell[data.u1]) : fe.p(neib[data.u1]);
        }

        if (interesting) {
            //std::cout << "asig: " << a_sig << "\n";
            //a_sig = vn > 0.0 ? cell[data.u1] : neib[data.u1];
            //std::cout << "asig: " <<a_sig << "\n\n";
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

    cell[data.u2] = cell[data.u1] - fluxes / cell.volume();
}

void Transfer::compute_slopes(Mesh& mesh) const {
    if (m_method == Method::MUSCLn ||
        m_method == Method::MUSCLn_CRP) {
        for (auto cell: mesh) {
            // Реконструкция в ячейке
            obj::plane plane{
                    .p = cell[data.p] + cell.center().dot(cell[data.n]),
                    .n = cell[data.n]
            };

            Vector3d grad = Vector3d::Zero();
            for (auto face: cell.faces()) {
                // Отрезок - грань
                obj::segment seg{
                        .v1 = face.vs(0),
                        .v2 = face.vs(1)
                };
                double a_sig = intersection2D::edge_fraction(seg, plane);
                grad += a_sig * face.area_n();
            }
            cell[data.grad] = grad / cell.volume();
        }
        return;
    }

    auto u1 = data.u1;
    auto get_state = [u1](Cell& cell) -> double {
        return cell[u1];
    };
    auto boundary_value = [](double u, const Vector3d& n, Boundary b) -> double {
        return u;
    };

    for (auto cell: mesh) {
        auto grad = gradient::LSM<double>(cell, get_state, boundary_value);
        cell[data.grad] = {grad.x, grad.y, 0.0};

        if (m_method == Method::MUSCL_MC || m_method == Method::MUSCL_MC_CRP) {
            auto lim_grad = gradient::limiting<double>(cell, m_limiter,
                    grad, get_state, boundary_value);

            cell[data.grad] = {lim_grad.x, lim_grad.y, 0.0};
        }
    }
}

void Transfer::update(Mesh &mesh, Direction dir) {
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

void Transfer::update_CRP(Mesh& mesh, Direction dir) const {
    // Считаем потоки
    mesh.for_each([&](Cell& cell) {
        fluxes_CRP(cell, dir);
    });

    // Обновляем слои
    mesh.for_each([this](Cell& cell) {
        cell[data.u1] = between(cell[data.u2], 0.0, 1.0);
        cell[data.u2] = 0.0;
    });

    // Без сглаживаний, чисто для реконструкции
    update_interface(mesh, 0);
}

void Transfer::update_VOF(Mesh& mesh, Direction dir) {
    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_VOF(cell, dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell[data.u1] =  between(cell[data.u2], 0.0, 1.0);
        cell[data.u2] = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_MUSCL(Mesh& mesh, Direction dir) const {
    compute_slopes(mesh);

    // Считаем потоки
    for (auto cell: mesh) {
        fluxes_MUSCL(cell, dir);
    }

    // Обновляем слои
    for (auto& cell: mesh) {
        cell[data.u1] = between(cell[data.u2], 0.0, 1.0);
        cell[data.u2] = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_WENO(Mesh& mesh, Direction dir) const {
    if (mesh.dim() == 3) {
        throw std::runtime_error("NO WENO");
    }

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
                    mesh(I - 2, j)[data.u1],
                    mesh(I - 1, j)[data.u1],
                    mesh(I + 0, j)[data.u1],
                    mesh(I + 1, j)[data.u1],
                    mesh(I + 2, j)[data.u1],
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
                        mesh(I - 2, j)[data.u1],
                        mesh(I - 1, j)[data.u1],
                        mesh(I + 0, j)[data.u1],
                        mesh(I + 1, j)[data.u1],
                        mesh(I + 2, j)[data.u1],
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
                        mesh(i, J - 2)[data.u1],
                        mesh(i, J - 1)[data.u1],
                        mesh(i, J + 0)[data.u1],
                        mesh(i, J + 1)[data.u1],
                        mesh(i, J + 2)[data.u1],
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
                        mesh(i, J - 2)[data.u1],
                        mesh(i, J - 1)[data.u1],
                        mesh(i, J + 0)[data.u1],
                        mesh(i, J + 1)[data.u1],
                        mesh(i, J + 2)[data.u1],
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

            cell[data.u2] = cell[data.u1] - fluxes / cell.volume();
        }
    }

    // Обновляем слои
    for (auto cell: mesh) {
        cell[data.u1] = between(cell[data.u2], 0.0, 1.0);
        cell[data.u2] = 0.0;
    }

    update_interface(mesh);
}

void Transfer::update_interface(Mesh& mesh, int smoothing) const {
    mesh.for_each([this](Cell& cell) {
        auto [p, n] = m_plic.plane(cell, 0);
        cell[data.p] = p;
        cell[data.n] = n;
    });
}

void Transfer::set_flags(Mesh& mesh) const {
    mesh.for_each([this](Cell& cell) {
        double min_val = cell[data.u1];
        double max_val = cell[data.u1];

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            min_val = std::min(min_val, face.neib(data.u1));
            max_val = std::max(max_val, face.neib(data.u1));
        }

        if (max_val < 1.0e-13 || min_val > 1.0 - 1.0e-13) {
            cell.set_flag(-1);
        }
        else if (max_val - min_val > 1.0e-13) {
            cell.set_flag(1);
        }
        else {
            cell.set_flag(-1);
        }
    });
}

Distributor Transfer::distributor() const {
    Distributor distr;

    auto u1 = data.u1;

    distr.split = [u1](const Cell& parent, Children &children) {
        for (auto child: children) {
            child[u1] = parent[u1];
        }
    };

    distr.merge = [u1](const Children &children, Cell& parent) {
        double sum = 0.0;
        for (auto child: children) {
            sum += child[u1] * child.volume();
        }
        parent[u1] = sum / parent.volume();
        if (parent[u1] < 1.0e-13) {
            parent[u1] = 0.0;
        }
        if (parent[u1] > 1.0 - 1.0e-13) {
            parent[u1] = 1.0;
        }
    };

    return distr;
}

Mesh Transfer::body(Mesh& mesh) const {
    auto empty_cell = [this](Cell& cell) -> bool {
        return cell[data.u1] <= 1.0e-12 || (cell[data.u1] < 0.5 && cell[data.n].isZero());
    };

    using Eigen::Vector3i;

    Vector3i count = mesh.sum([&empty_cell](Cell& cell) -> Vector3i {
        if (empty_cell(cell)) {
            return {0, 0, 0};
        }
        return {1, cell.face_count() + 1, cell.node_count() + 2};
    }, Vector3i{0, 0, 0});

    int n_cells = count[0];
    int n_faces = count[1];
    int n_nodes = count[2];

    Mesh clipped = Mesh::PolySet(mesh.dim());
    clipped.locals().reserve(n_cells, n_faces, n_nodes);

    if (mesh.dim() == 2) {
        for (auto cell: mesh) {
            if (empty_cell(cell)) {
                continue;
            }

            if (cell[data.u1] > 1.0 - 1.0e-12) {
                clipped.push_back(cell.polygon());
                continue;
            }

            Vector3d P = cell.center() + cell[data.p] * cell[data.n];
            if (cell[data.n].isZero()) {
                double d = 0.5 * std::sqrt(cell[data.u1] * cell.volume());
                Polygon poly = {
                    P + Vector3d{-d, -d, 0.0},
                    P + Vector3d{+d, -d, 0.0},
                    P + Vector3d{+d, +d, 0.0},
                    P + Vector3d{-d, +d, 0.0},
                };
                clipped.push_back(poly);
            }
            else {
                auto poly = cell.polygon();
                auto part = poly.clip(P, cell[data.n]);
                clipped.push_back(part);
            }
        }
        return clipped;
    }
    else {
        for (auto& cell: mesh) {
            if (empty_cell(cell)) {
                continue;
            }
            if (cell[data.u1] > 1.0 - 1.0e-12) {
                clipped.push_back(cell.polyhedron());
                continue;
            }

            Vector3d P = cell.center() + cell[data.p] * cell[data.n];
            auto poly = cell.polyhedron();
            auto clip = poly.clip(P, cell[data.n]);
            if (!clip.empty()) {
                clipped.push_back(clip);
            }
        }
    }
    return clipped;
}

} // namespace zephyr::math