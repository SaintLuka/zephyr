#pragma once

namespace zephyr { namespace phys {

struct dRdE {
    double value;
    double dR;
    double dE;

    dRdE(double value, double dR = 0.0, double dE = 0.0)
        : value(value), dR(dR), dE(dE) { }

    operator double() const {
        return value;
    };
};

/*
struct dRdP {
    double value;
    double dR;
    double dP;

    dRdP(double value, double dR = 0.0, double dP = 0.0)
            : value(value), dR(dR), dP(dP) { }

    operator double() const {
        return value;
    };
};

struct dRdT {
    double value;
    double dR;
    double dT;

    dRdT(double value, double dR = 0.0, double dT = 0.0)
            : value(value), dR(dR), dT(dT) { }

    operator double() const {
        return value;
    };
};

struct dPdT {
    double value;
    double dP;
    double dT;

    dPdT(double value, double dP = 0.0, double dT = 0.0)
            : value(value), dP(dP), dT(dT) { }

    operator double() const {
        return value;
    };
};
 */

}
}