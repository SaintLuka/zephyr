#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace phys {

class SodTest {
public:

    SodTest() {
        m_eos = IdealGas(1.4);

        rL = 1.0; rR = 0.125;
        pL = 1.0; pL = 0.1;
        uL = 0.0; uR = 0.0;
        eL = m_eos.energy_rp(rL, pL);
        eR = m_eos.energy_rp(rR, pR);

        x_jump = 0.5;
    }

    IdealGas &eos() { return m_eos; }

    const IdealGas &eos() const { return m_eos; }


    double xmin() const { return 0.0; }

    double xmax() const { return 1.0; }

    double max_time() const { return 0.2; }


    double density(const Vector3d &r) const {
        return r.x() < x_jump ? rL : rR;
    }

    double pressure(const Vector3d &r) const {
        return r.x() < x_jump ? pL : pR;
    }

    Vector3d velocity(const Vector3d &r) const {
        return r.x() < x_jump ? Vector3d(uL, 0.0, 0.0) : Vector3d(uR, 0.0, 0.0);
    }

    double energy(const Vector3d &r) const {
        return r.x() < x_jump ? eL : eR;
    }

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
