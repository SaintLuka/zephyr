#include <iostream>
#include <vector>

#include <zephyr/phys/tests/sedov3D.h>

#include <zephyr/math/calc/integrals.h>
#include <zephyr/utils/matplotlib.h>

namespace plt = zephyr::utils::matplotlib;

namespace zephyr::phys {

inline double sqr(double x) {
    return x * x;
}

using ScalarFunc = std::function<double(double)>;


SedovBlast3D::SedovBlast3D(const params& p) {
    gamma = p.gamma;
    rho0  = p.rho0;
    E     = p.E;
    
    eos = IdealGas::create(gamma, 1.0);
    
    // Параметры решения
    v1 = -(13.0 * gamma * gamma - 7.0 * gamma + 12.0) / ((3.0 * gamma - 1.0) * (2.0 * gamma + 1.0));
    v2 = 5.0 * (gamma - 1.0) / (2.0 * gamma + 1.0);
    v3 = 3.0 / (2.0 * gamma + 1.0);
    v4 = -v1 / (2.0 - gamma);
    v5 = -2.0 / (2.0 - gamma);

    // Ищем коэффициент beta
    auto Integrand = [&](double xi) -> double {
        double V = V_xi(xi);
        double Z = Z_xi(xi, V);
        double G = G_xi(xi, V);

        double res = std::abs(G * (0.5 * V * V + Z / (gamma * (gamma - 1.0)))) * sqr(sqr(xi));
        return std::isnan(res) ? 0.0 : res;
    };

    double K = 16.0 * M_PI / 25.0;
    double b_int = zephyr::math::integrate(Integrand, 0.0, 1.0, 10000);
    beta = std::pow(K * b_int , -0.2);
    
    finish = time_by_radius(1.0);
    
    // TODO: Построить сплайны для V(xi), Z(xi), G(xi).
}

double SedovBlast3D::density(const Vector3d &r, double t) const {
    SedovBlast3D* st = (SedovBlast3D*)(this);
    return st->rho(r.norm(), t);
}

double SedovBlast3D::pressure(const Vector3d &r, double t) const {
    SedovBlast3D* st = (SedovBlast3D*)(this);
    return std::max(0.001, st->p(r.norm(), t));
}

Vector3d SedovBlast3D::velocity(const Vector3d &r, double t) const {
    SedovBlast3D* st = (SedovBlast3D*)(this);
    return r / r.norm() * st->v(r.norm(), t);
}

double SedovBlast3D::energy(const Vector3d &r, double t) const {
    return eos->energy_rP(density(r, t), pressure(r, t));
}


#define TOLERANCE 1.e-9

/*check to be sure our initial guesses are valid*/
int check_initial_values(ScalarFunc func, double x_min, double x_max) {
    double y_min = func(x_min);
    double y_max = func(x_max);

    // printf("CIV x_min %15.14e y_min %15.14e x_max %15.14e y_max %15.14e\n",x_min,y_min,x_max,y_max);

    /*check for zero crossing*/
    if (y_min * y_max >= 0.0) {
        printf("No zero crossing found in range [%15.14e, %15.14e].\n", x_min, x_max);
        printf("f(%15.14e) = %15.14e, f(%15.14e) = %15.14e\n", x_min, y_min, x_max, y_max);
        printf("Change initial values.\n");
        throw std::runtime_error("Sedov Blast error: bisection failed");
    }

    /*If x_min is a root, then return a flag == 1*/
    if (std::abs(y_min) < TOLERANCE)
        return 1;

    /*If x_max is a root, then return a flag == 2*/
    if (std::abs(y_max) < TOLERANCE)
        return 2;

    /*if we reach this point, the bracket is valid.*/
    return 0;
}

/* function that uses bisection search to find a root */
double bisection_root_finding(ScalarFunc func, double x_min_start, double x_max_start) {
    double x_min = x_min_start;  /* minimum x in bracket */
    double x_mid;              /* 1/2 between x_min and x_max */
    double x_max = x_max_start;  /* maximum x in bracket */
    double y_min;              /* f at x_min */
    double y_mid;              /* f at x_mid */
    double y_max;              /* f at x_max */

    int flag = 0;          /* a flag for triggering */
    int iter = 0;          /* current number of iterations */
    int imax = 100000;          /* maximum number of iterations */

    /* check our initial bracket to make sure it's valid */
    flag = check_initial_values(func, x_min, x_max);
    if (flag) {
        /* we made a lucky guess and we're done */
        if (flag == 1)        /*x_min is the root*/
            return x_min;
        if (flag == 2)     /*x_max is the root*/
            return x_max;
    }

    /*instead, if flag==0, we continue on our way to find the root*/

    /* loop until we find the root or we reach the maximum number of iterations*/
    while (!flag) {
        /* find value midway between bracket values */
        x_mid = 0.5 * (x_min + x_max);

        /* determine the function value at this point */
        y_mid = func(x_mid);

        /* if the absolute value of f(x) is < tolerance, then exit */
        if (fabs(y_mid) < TOLERANCE || fabs(x_max - x_min) < 1.0e-14) {
            /* trigger an exit */
            flag = 1;

        } else {

            /*x_mid is not a root*/

            y_min = func(x_min);
            y_max = func(x_max);

            /* print guesses */
            //print_guess(x_min,y_min);
            //print_guess(x_mid,y_mid);
            //print_guess(x_max,y_max);

            //printf("xmin %15.14e ymin %15.14e xmid %15.14e ymid %15.14e xmax %15.14e ymax %15.14e\n",x_min,y_min,x_mid,y_mid,x_max,y_max);

            /* if the product of the function at the midpoint and at one of the endpoints */
            /* is greater than zero, replace this end point */
            //if( f(x_min,params)*f(x_mid,params) > 0)
            if (y_min * y_mid > 0) {
                /*replace x_min with x_mid */
                x_min = x_mid;
            } else {
                /*replace x_max with x_mid */
                x_max = x_mid;
            }


            /* count the iteration */
            iter++;

            /* if we have exceeded the max number of iterations, exit*/
            if (iter >= imax) {
                printf("Exceeded %d iterations, exiting (x=%f, f=%f)...\n", imax, x_mid, y_mid);
                fflush(stdout);
                exit(-1);
            }
        }
    }

    return x_mid;
}

double SedovBlast3D::eta_V(double V) const {
    double A = 1.0 / sqr(0.5 * (gamma + 1.0) * V);
    double B = std::pow((gamma + 1.0) / (7.0 - gamma) * (5.0 - (3.0 * gamma - 1.0) * V), v1);
    double C = std::pow((gamma + 1.0) / (gamma - 1.0) * (gamma * V - 1.0), v2);

    return A * B * C;
}

double SedovBlast3D::V_xi(double xi) const {
    return V_eta(std::pow(xi, 5));
}

double SedovBlast3D::V_eta(double eta) const {
    double vmin = 1.0 / gamma;
    double vmax = 2.0 / (gamma + 1.0);

    auto f_V_eta = [&](double V) -> double {
        return eta_V(V) - eta;
    };

    double eps = 1.0e-8;
    double eta_ref = eta;
    if (eta < eps) {
        eta = eps;
        // TODO: вынести метод бисекции в методы?
        // TODO: Причесать функцию, сейчас V(eta) странный вид имеет
        double near_vmin = bisection_root_finding(f_V_eta, vmin, vmax);
        return vmin + (near_vmin - vmin) * eta_ref / eps;

    } else if (eta > 1.0 - eps) {
        eta = 1.0 - eps;
        double near_vmax = bisection_root_finding(f_V_eta, vmin, vmax);
        return vmax + (vmax - near_vmax) * (eta_ref - 1.0) / eps;
    } else {
        return bisection_root_finding(f_V_eta, vmin, vmax);
    }
}

double SedovBlast3D::Z_xi(double xi, double V) const {
    if (std::isnan(V)) V = V_xi(xi);

    return gamma * (gamma - 1.0) * (1.0 - V) * sqr(V) / (2.0 * (gamma * V - 1.0));
}

double SedovBlast3D::G_xi(double xi, double V) const {
    if (std::isnan(V)) V = V_xi(xi);

    double A = (gamma + 1.0) / (gamma - 1.0);
    double B = std::pow(A * (gamma * V - 1.0), v3);
    double C = std::pow((gamma + 1.0) / (7.0 - gamma) * (5.0 - (3.0 * gamma - 1.0) * V), v4);
    double D = std::pow(A * (1.0 - V), v5);

    return A * B * C * D;
}

double SedovBlast3D::r_shock(double t) const {
    return beta * std::pow(E * t * t / rho0, 0.2);
}

double SedovBlast3D::time_by_radius(double r) const {
    return std::sqrt(rho0 / E * std::pow(r / beta, 5));
}

double SedovBlast3D::c_squared(double r, double t) const {
    double xi = r / r_shock(t);
    return sqr(0.4 * r / t) * Z_xi(xi);
}

double SedovBlast3D::rho(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? rho0 * G_xi(xi) : rho0;
}

double SedovBlast3D::v(double r, double t) const {
    double xi = r / r_shock(t);
    return xi < 1.0 ? (0.4 * r / t) * V_xi(xi) : 0.0;
}

double SedovBlast3D::p(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? rho(r, t) * c_squared(r, t) / gamma : 0.0;
}

} // namespace zephyr::phys