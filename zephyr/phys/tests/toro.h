#pragma once

#include <geom/types/vector.h>
#include <physics/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

class ToroTest {
public:

    ToroTest() : m_eos("Air") {}

    IdealGas &eos() { return m_eos; }


    double xmin() const { return -5.0; }

    double xmax() const { return 5.0; }

    double max_time() const { return 2.0; }


    double density(const geom::Vector &r) const {
        return r.x() < x_jump ? 3.857143 : 1.0 + 0.2 * std::sin(5.0 * r.x());
    }

    double pressure(const geom::Vector &r) const {
        return r.x() < x_jump ? 10.33333 : 1.0;
    }

    geom::Vector velocity(const geom::Vector &r) const {
        return r.x() < x_jump ? geom::Vector(2.629369, 0.0, 0.0) : geom::Vector::Zero();
    }

    double energy(const geom::Vector &r) const {
        return m_eos.energy_rp(density(r), pressure(r));
    }

private:
    IdealGas m_eos;
    const double x_jump = -4.0;
};

}
}
