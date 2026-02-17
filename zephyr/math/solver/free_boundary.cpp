#include <zephyr/math/solver/free_boundary.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>

#include "zephyr/geom/curves/interpolant.h"
#include "zephyr/math/funcs.h"

namespace zephyr::math {

using namespace geom;
using namespace smf;

using utils::threads;
using utils::mpi;

constexpr double eps = 1.0e-15;

double FreeBoundary::bound_pos(double t) const {
    return 0.5 *t;
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

        bool interesting = cell.index() == 250;

        // Переменная для потока
        Flux flux = Flux::Zero();
        double F_a = 0;
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

            // Вектор скорости
            //Vector3d vel = {bound_vel(current_time()), 0.0, 0.0};
            Vector3d vel = {2.0, 0.0, 0.0};

            // Нормальная составляющая
            double ui = vel.dot(normal);

            double aL = cell[part.alpha];
            double aR = face.neib(part.alpha);
            if (interesting) {
            }

            // Численный поток на грани
            auto [FaL, FL, FaR, FR] = calc_flux(zm, zp, *m_eos, aL, aR, ui, dt(), cell.hx());
            FL.to_global(normal);

            F_a += FaL * face.area();
            flux.arr() += FL.arr() * face.area();

            if (interesting) {
                std::cout << "qc: " << z_c << "\n";
                std::cout << "zm: " << zm << "\n";
                std::cout << "zp: " << zp << "\n";
                std::cout << "Fal, Far, FL: " << FaL << " " << FaR << " " << FL << "\n";
            }
        }

        cell[part.a_next] = cell[part.alpha] - (dt() / cell.volume()) * F_a;

        // Ячейка без газа
        if (cell[part.a_next] > 1.0 - eps) {
            cell[part.a_next] = 1.0;
            cell[part.next] = PState::NaN();
            return;
        }

        if (cell[part.a_next] < eps) {
            cell[part.a_next] = 0.0;
        }

        // Консервативный вектор в ячейке
        QState q_c = QState::Zero();

        // Новое значение примитивных переменных
        if (cell[part.alpha] < 1.0) {
            q_c = QState(z_c);
            q_c.arr() *= (1.0 - cell[part.alpha]);
        }

        if (interesting) {
            std::cout << "alpha: " << cell[part.alpha] << " " << cell[part.a_next] << "\n";
            std::cout << "q_c: " << q_c << "\n";
        }

        q_c.arr() -= (dt() / cell.volume()) * flux.arr();
        q_c.arr() /= (1.0 - cell[part.a_next]);

        cell[part.next] = PState(q_c, *m_eos);
    });
}

std::tuple<double, Flux, double, Flux> FreeBoundary::calc_flux(
    PState zL, PState zR, const Eos &eos,
    double alphaL, double alphaR, double u, double dt, double cellx) {
    // Грань между двумя газовыми ячейками
    if(alphaL == 0.0 && alphaR == 0.0) {
        Flux F = HLLC::calc_flux(zL, zR, eos);
        return {0.0, F, 0.0, F};
    }

    // Грань между двумя пустыми ячейками
    if (alphaL == 1.0 && alphaR == 1.0) {
        Flux F = Flux::Zero();
        return {0.0, F, 0.0, F};
    }

    // Для положительной скорости
    double P_L = zL.pressure;

    // "компенсационный поток"
    // u - скорость поршня, P - давление на поршне
    auto G = [](double u, double P) -> Flux {
        return Flux(0.0, {P, 0.0, 0.0}, u * P);
    };

    double u_m = std::min(u, 0.0);  // u_минус
    double u_p = std::max(u, 0.0);  // u_плюс

    // Всё забыли. Сводим к задаче, где слева поршень, а справа газ.
    bool inversed = false;
    if (alphaL < alphaR) {
        inversed = true;
        u *= -1.0;
        std::swap(alphaL, alphaR);
        std::swap(zL, zR);
        zL.inverse();
        zR.inverse();
    }

    // теперь alphaL > alphaR

    double Fa_L{NAN};
    double Fa_R{NAN};
    Flux F_L = Flux::NaN();
    Flux F_R = Flux::NaN();

    // Слева только поршень, справа есть газ
    if (alphaL == 1.0) {
        auto wc = HLLC::wave_config_u_R(u, eos, zR);
        PState z(wc.QsR, eos); // состояние на поршне

        // Сжатие газа
        if (u >= 0.0) {
            Fa_L = 0.0;
            Fa_R = u;
            F_L = Flux::Zero();
            F_R = G(u, z.P());
        }
        else {
            // Разрежение газа
            // Справа газ и поршень (0 < alphaR < 1)
            double delta = alphaR * cellx;

            double tau1 = std::min(-delta/ (u * dt), 1.0);
            double tau2 = 1.0 - tau1;

            Fa_L = -tau2 * u;
            Fa_R = -tau1 * u;

            F_L = tau2 * wc.FsR.arr();
            F_L.momentum.x() -= tau2 * z.P();
            F_L.energy       -= tau2 * z.P() * u;

            F_R = tau2 * wc.FsR.arr();
            F_R.momentum.x() -= tau1 * z.P();
            F_R.energy       -= tau1 * z.P() * u;
        }
    }
    else {
        throw std::runtime_error("Enough");
    }

    if (!inversed) {
        return {Fa_L, F_L, Fa_R, F_R};
    }
    else {
        Fa_L *= -1.0;
        Fa_R *= -1.0;
        F_L.inverse();
        F_R.inverse();
        return {Fa_R, F_R, Fa_L, F_L};
    }

    throw std::runtime_error("Enough");

        Flux F{};
        double FalphaL = 0;
        double FalphaR = 0;

    if (std::min(alphaL, alphaR) == 0.0 && std::max(alphaL, alphaR) == 1.0) {
        if (u < 0.0) {
            FalphaL = u * dt;
            Flux F_add{0.0, {P_L, 0.0, 0.0}, P_L * u};
            F.arr() += F_add.arr();
        }
        else {

        }
    }



    else if ((0.0 < alphaL && alphaL < 1.0) && alphaR == 1.0 && u < 0.0) {
        if (-u * dt <= cellx - cellx * alphaL) {
            FalphaL = u * dt;
        } else {
            FalphaL = -(cellx - cellx * alphaL);
        }
        // F = 0 - потока газа через границу нет
    }
    else if (alphaL == 0.0 && (0.0 < alphaR && alphaR < 1.0) && u < 0.0) {
        if (-u * dt > cellx - cellx * alphaR) {
            FalphaL = -(-u * dt - (cellx - cellx * alphaR));
        }
        F = crp_flux_inverse(zL, zR, eos, cellx * alphaR, u, dt);
        Flux F_add{0.0, {-P_L, 0.0, 0.0}, P_L * u};
        F.arr() += F_add.arr();
    }
    else if ((0.0 < alphaL && alphaL < 1.0) &&
        std::fabs(alphaR) < std::numeric_limits<double>::epsilon() && u > 0) {
        F = crp_flux_classic(zL, zR, eos, cellx * alphaL, u, dt);
        Flux F_add{0.0, {P_L, 0.0, 0.0}, P_L * u};
        F.arr() += F_add.arr();
    }

    // Для отрицательной скорости
    if (alphaL == 1.0 && alphaR == 0.0 && u < 0) {
        FalphaL = -u * dt;
    }
    else if (alphaL == 1.0 && (0.0 < alphaR && alphaR < 1.0) && u < 0) {
        if (-u * dt > cellx * alphaR) {
            FalphaL = -u * dt - (cellx * alphaR);
        }
    }
    else if ((0.0 < alphaL && alphaL < 1.0) && alphaR == 0.0 && u < 0) {
        if (-u * dt <= cellx * alphaL) {
            FalphaL = -u * dt;
        } else {
            FalphaL = cellx * alphaL;
        }
    }

    return {FalphaL, F, FalphaR, F};
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