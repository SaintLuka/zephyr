#include <zephyr/math/solver/free_boundary.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>

#include "zephyr/geom/curves/interpolant.h"

namespace zephyr::math {

using namespace geom;
using namespace smf;

using utils::threads;
using utils::mpi;

double FreeBoundary::bound_pos(double t) const {
    return t;
    // return std::sin(t);
}

double FreeBoundary::bound_vel(double t) const {
    return (bound_pos(t + m_dt) - bound_pos(t)) / m_dt;
}

FreeBoundary::FreeBoundary(Eos::Ptr eos) : m_eos(eos) {
    m_CFL = 0.5;
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

FreeBoundary::Parts FreeBoundary::add_types(EuMesh& mesh) {
    part.init = mesh.add<PState>("init");
    part.next = mesh.add<PState>("next");
    part.alpha = mesh.add<double>("alpha");
    part.a_next = mesh.add<double>("alpha(next)");
    return part;
}

void FreeBoundary::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

double FreeBoundary::CFL() const {
    return m_CFL;
}

double FreeBoundary::dt() const {
    return m_dt;
}

double FreeBoundary::current_time() const {
    return curr_time;
}

void FreeBoundary::set_max_dt(double dt) {
    m_max_dt = dt;
}

void FreeBoundary::set_curr_time(double t) {
    curr_time = t;
}

namespace{
    PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
        if (flag != Boundary::WALL) {
            return zc;
        }

        PState zn(zc);
        Vector3d Vn = normal * zc.velocity.dot(normal);
        zn.velocity = zc.velocity - 2.0 * Vn;

        return zn;
    } 
}


void FreeBoundary::update(EuMesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    mesh.sync(part.init);
    mesh.sync(part.alpha);
    fluxes(mesh);

    // Обновляем слои
    swap(mesh);
}

void FreeBoundary::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        //double c = m_eos->sound_speed_rP(cell(part.density), cell(part.pressure));
        //return cell.incircle_diameter() / (cell(part.velocity).norm() + c);
        double c = m_eos->sound_speed_rP(cell[part.init].density, cell[part.init].pressure);
        return cell.incircle_diameter() / (cell[part.init].velocity.norm() + c);
    });

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void FreeBoundary::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux, F{};
        double Falpha = 0;
        for (auto &face: cell.faces(mesh::Direction::X)) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = face.neib(part.init);
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Численный поток на грани
            // Flux loc_flux = HLLC::calc_flux(zm, zp, *m_eos);
            // loc_flux.to_global(normal);

            double FalphaL, FalphaR;
            Flux Fi{};

            // a.dot(b); a.norm(); a.normalize(); b = a.normalized();

            double u = bound_vel(current_time());
            Vector3d normal_normalized = normal.normalized();
            Vector3d u_vector = Vector3d(1, 0, 0);
            u_vector = u_vector.normalized();
            double dot_vectors = u_vector.x() * normal_normalized.x() + u_vector.y() * normal_normalized.y() + u_vector.z() * normal_normalized.z();
            double ui = u * dot_vectors;

            double alphaL = cell(part.alpha);
            double alphaR = face.neib(part.alpha);

            std::tuple<double, double, Flux> alpha_fluxes;
            if (alphaL == alphaR) {
                alpha_fluxes = FreeBoundary::flux(zm, zp, *m_eos, alphaL, alphaR, ui, dt(), cell.hx());
            } else {
                alpha_fluxes = FreeBoundary::flux(zm, zp, *m_eos, alphaL, alphaR, ui, dt(), cell.hx());
            }
            std::tie(FalphaL, FalphaR, Fi) = alpha_fluxes;
            Falpha += FalphaL * face.area();
            Falpha += FalphaR * face.area();

            Fi.to_global(normal);
            F.arr() += Fi.arr() * face.area();

            // flux.arr() += loc_flux.arr() * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        // q_c.arr() -= (m_dt / cell.volume()) * flux.arr();

        cell[part.a_next] = cell[part.alpha] - (Falpha / cell.volume());
        if (std::abs(cell[part.a_next] - 1.0) < std::numeric_limits<double>::epsilon()) {
            cell[part.a_next] = 1.0;
        } else if (std::abs(cell[part.a_next]) < std::numeric_limits<double>::epsilon()) {
            cell[part.a_next] = 0.0;
        }

        // Новое значение примитивных переменных
        q_c.arr() *= (1 - cell[part.alpha]);
        q_c.arr() -= (dt() / cell.volume()) * F.arr();
        if (std::abs(cell[part.a_next] - 1.0) < std::numeric_limits<double>::epsilon()) {
            QState nan_q_c(true);
            q_c.arr() = nan_q_c.arr();
        } else {
            q_c.arr() /= (1 - cell[part.a_next]);
        }
        cell[part.next] = PState(q_c, *m_eos);
    });
}

std::tuple<double, double, Flux> FreeBoundary::flux(const PState &zL, const PState &zR, const Eos &eos,
        double alphaL, double alphaR, double u, double dt, double cellx){
    double FalphaL = 0;
    double FalphaR = 0;
    Flux F{};
    if(std::fabs(alphaL) < std::numeric_limits<double>::epsilon() &&
        std::fabs(alphaR) < std::numeric_limits<double>::epsilon()){
        F = HLLC::calc_flux(zL, zR, eos);
        return {FalphaL, FalphaR, F};
    }

    if(std::fabs(alphaL - 1.0) < std::numeric_limits<double>::epsilon() &&
        std::fabs(alphaR - 1.0) < std::numeric_limits<double>::epsilon()){
        return {FalphaL, FalphaR, F};
    }

    // Для положительной скорости
    const double &P_L = zL.pressure;
    if (std::fabs(alphaL) < std::numeric_limits<double>::epsilon() &&
        std::fabs(alphaR - 1.0) < std::numeric_limits<double>::epsilon() && u < 0) {
        FalphaL = u * dt;
        Flux F_add{0.0, {P_L, 0.0, 0.0}, P_L * u};
        F.arr() += F_add.arr();
    } else if ((alphaL > 0 && alphaL < 1) &&
        std::fabs(alphaR - 1.0) < std::numeric_limits<double>::epsilon() && u < 0) {
        if (-u * dt <= cellx - cellx * alphaL) {
            FalphaL = u * dt;
        } else {
            FalphaL = -(cellx - cellx * alphaL);
        }
        // F = 0 - потока газа через границу нет
    } else if (std::fabs(alphaL) < std::numeric_limits<double>::epsilon() &&
        (alphaR > 0 && alphaR < 1) && u < 0) {
        if (-u * dt > cellx - cellx * alphaR) {
            FalphaL = -(-u * dt - (cellx - cellx * alphaR));
        }
        F = crp_flux_inverse(zL, zR, eos, cellx * alphaR, u, dt);
        Flux F_add{0.0, {-P_L, 0.0, 0.0}, P_L * u};
        F.arr() += F_add.arr();
    } else if ((alphaL > 0 && alphaL < 1) &&
        std::fabs(alphaR) < std::numeric_limits<double>::epsilon() && u > 0) {
        F = crp_flux_classic(zL, zR, eos, cellx * alphaL, u, dt);
        Flux F_add{0.0, {P_L, 0.0, 0.0}, P_L * u};
        F.arr() += F_add.arr();
    }

    // Для отрицательной скорости
    if (std::fabs(alphaL - 1.0) < std::numeric_limits<double>::epsilon() &&
        std::fabs(alphaR) < std::numeric_limits<double>::epsilon() && u < 0) {
        FalphaL = -u * dt;
    } else if (std::fabs(alphaL - 1.0) < std::numeric_limits<double>::epsilon() &&
        (alphaR > 0 && alphaR < 1) && u < 0) {
        if (-u * dt > cellx * alphaR) {
            FalphaL = -u * dt - (cellx * alphaR);
        }
    } else if ((alphaL > 0 && alphaL < 1) &&
        std::fabs(alphaR) < std::numeric_limits<double>::epsilon() && u < 0) {
        if (-u * dt <= cellx * alphaL) {
            FalphaL = -u * dt;
        } else {
            FalphaL = cellx * alphaL;
        }
    }


    return {FalphaL, FalphaR, F};
}

// Точка в плоскости (x, t)
struct Point {
    double x, t;
};

// Характеристика (x, t) -- точка, S -- скорость
struct Char {
    double x, t, S;

    // Пересечение характеристики с прямой x = 0, время t.
    double edge_t() const {
        double ts = (S * t - x) / S;
        if (ts <= t) {
            return std::numeric_limits<double>::infinity();
        } else {
            return ts;
        }
    }

    // Пересечение характеристик, находятся только пересечения во времени дальше,
    // чем время обеих точек исходных характеристик. В обратном случае считается,
    // что характеристики пересекаются на бесконечности t = +inf.
    Point cross(const Char &c) const {
        const double inf = std::numeric_limits<double>::infinity();

        if (S == c.S) {
            return {.x = 0.5 * (x + c.x), .t = inf};
        }

        Point p = {.x = NAN, .t = NAN};
        p.t = (c.x - x + S * t - c.S * c.t) / (S - c.S);

        if (p.t <= std::max(t, c.t)) {
            return {.x = 0.5 * (x + c.x), .t = inf};
        }

        p.x = 0.5 * (x + c.x + S * (p.t - t) + c.S * (p.t - c.t));

        return p;
    }
};

Flux FreeBoundary::crp_flux_classic(const PState& zL, const PState& zR,
                                    const Eos& eos, double delta, double u, double dt) {
    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = u};

    auto[S_0L, S_0C, S_0R, Q_s0L, F_s0L, Q_s0R, F_s0R] = HLLC::wave_config(eos, zL, eos, zR);
    Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};

    // Первое взаимодействие
    Point O1 = C_0.cross(C_0L);
    if (S_0C >= 0.0) {
        const Flux& F0 = F_s0L;

        if (O1.t >= dt) {
            return F0;
        }

        PState z_s0L(Q_s0L, eos);

        // Характеристики из точки O1
        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config_u_R(u, eos, z_s0L);
        Char C_1L = {.x = O1.x, .t = O1.t, S_1L};
        Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

        double tau1 = C_1R.edge_t() / dt;
        if (tau1 >= 1.0) {
            return F0;
        }

        const Flux& F1 = F_s1R;

        double tau2 = C_1L.edge_t() / dt;
        if (tau2 >= 1.0) {
            return tau1 * F0.arr() + (1.0 - tau1) * F1.arr();
        }

        const Flux& F2{};

        return tau1 * F0.arr() + (tau2 - tau1) * F1.arr() + (1.0 - tau2) * F2.arr();
    }else {
        const Flux& F0 = F_s0R;

        if (O1.t >= dt) {
            return F0;
        }

        PState z_s0L(Q_s0L, eos);

        Char C_0C = {.x = 0.0, .t = 0.0, .S = S_0C};
        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config_u_R(u, eos, z_s0L);
        Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

        Point O2 = C_1R.cross(C_0C);

        if (O2.t > dt) {
            return F0;
        }

        PState z_s1R(Q_s1R, eos);
        PState z_s0R(Q_s0R, eos);
        // auto[S_2L, S_2C, S_2R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(eos, z_s1R, eos, z_s0R);
        // Char C_2C = {.x = O2.x, .t = O2.t, .S = S_2C};
        // Char C_2R = {.x = O2.x, .t = O2.t, .S = S_2R};
        //
        // double tau1 = C_2R.edge_t() / dt;
        // if (tau1 >= 1.0) {
        //     return F0;
        // }
        //
        // const Flux& F1 = F_s2R;
        //
        // double tau2 = C_2C.edge_t() / dt;
        // if (tau2 >= 1.0) {
        //     return tau1 * F0.arr() + (1.0 - tau1) * F1.arr();
        // }
        //
        // const Flux& F2 = F_s2L;
        //
        // return tau1 * F0.arr() + (tau2 - tau1) * F1.arr() + (1.0 - tau2) * F2.arr();
        auto[S_2L, S_2R, Q_s, F_s] = HLL::wave_config(eos, z_s1R, z_s0R);
        Char C_2L = {.x = O2.x, .t = O2.t, .S = S_2L};
        Char C_2R = {.x = O2.x, .t = O2.t, .S = S_2R};

        double tau1 = C_2R.edge_t() / dt;
        if (tau1 >= 1.0) {
            return F0;
        }

        const Flux& F1 = F_s;

        double tau2 = C_2L.edge_t() / dt;
        if (tau2 >= 1.0) {
            return tau1 * F0.arr() + (1.0 - tau1) * F1.arr();
        }

        const Flux& F2{};

        return tau1 * F0.arr() + (tau2 - tau1) * F1.arr() + (1.0 - tau2) * F2.arr();
    }
}

Flux FreeBoundary::crp_flux_inverse(const PState& zL, const PState& zR,
                                    const Eos& eos, double delta, double u, double dt) {
    auto& zL_i = const_cast<PState&>(zL);
    auto& zR_i = const_cast<PState&>(zR);

    zL_i.inverse();
    zR_i.inverse();

    auto flux = crp_flux_classic(zR_i, zL_i, eos, delta, -u, dt);
    flux.inverse();

    zL_i.inverse();
    zR_i.inverse();

    return flux;
}

void FreeBoundary::swap(EuMesh &mesh) const {
    mesh.swap(part.init, part.next);
    mesh.swap(part.alpha, part.a_next);
}

} // namespace zephyr::math