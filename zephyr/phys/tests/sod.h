#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace phys {

class SodTest {
public:

    SodTest() {
        m_eos = IdealGas(1.4);

        rL = 1.0; rR = 0.125;
        pL = 1.0; pR = 0.1;
        uL = 0.0; uR = 0.0;
        eL = m_eos.energy_rp(rL, pL);
        eR = m_eos.energy_rp(rR, pR);

        x_jump = 0.5;
    }

    void inverse() {
        std::swap(rL, rR);
        std::swap(uL, uR);
        std::swap(pL, pR);
        std::swap(eL, eR);
        uL *= -1.0;
        uR *= -1.0;
    }

    IdealGas &eos() { return m_eos; }

    const IdealGas &eos() const { return m_eos; }

    double xmin() const { return 0.0; }

    double xmax() const { return 1.0; }

    double max_time() const { return 0.2; }


    double density(const double &x) const { return x < x_jump ? rL : rR; }

    Vector3d velocity(const double &x) const { return {x < x_jump ? uL : uR, 0.0, 0.0}; }

    double pressure(const double &x) const { return x < x_jump ? pL : pR; }

    double energy(const double &x) const { return x < x_jump ? eL : eR; }


    double density(const Vector3d &r) const { return density(r.x()); }

    Vector3d velocity(const Vector3d &r) const { return velocity(r.x()); }

    double pressure(const Vector3d &r) const { return pressure(r.x()); }

    double energy(const Vector3d &r) const { return energy(r.x()); }

private:
    IdealGas m_eos;
    double x_jump;
    double rL, rR;
    double uL, uR;
    double pL, pR;
    double eL, eR;
};

}
}
