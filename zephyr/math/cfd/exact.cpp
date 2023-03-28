#include <zephyr/math/cfd/exact.h>

namespace zephyr { namespace math {

using namespace zephyr::phys;

ExactRiemannSolver::ExactRiemannSolver(
        const PState &zL, const PState &zR,
        const StiffenedGas &gas) :
        ExactRiemannSolver(zL, zR, gas, gas) {
}

ExactRiemannSolver::ExactRiemannSolver(
        const PState &zL, const PState &zR,
        const StiffenedGas &gasL, const StiffenedGas &gasR) :
        gasL(gasL), gasR(gasR),
        rL(zL.density), pL(zL.pressure), eL(zL.energy),
        uL(zL.velocity.x()), vL(zL.velocity.y()), wL(zL.velocity.z()),
        rR(zR.density), pR(zR.pressure), eR(zR.energy),
        uR(zR.velocity.x()), vR(zR.velocity.y()), wR(zR.velocity.z()) {
}

void ExactRiemannSolver::set_left(double density, Vector3d velocity, double pressure,
                                  double energy, const StiffenedGas &gas) {
    gasL = gas;
    rL = density;
    uL = velocity.x();
    vL = velocity.y();
    wL = velocity.z();
    pL = pressure;
    eL = energy;
}

void ExactRiemannSolver::set_right(double density, Vector3d velocity, double pressure,
                                   double energy, const StiffenedGas &gas) {
    gasR = gas;
    rR = density;
    uR = velocity.x();
    vR = velocity.y();
    wR = velocity.z();
    pR = pressure;
    eR = energy;
}


void ExactRiemannSolver::compute() {

}


int ExactRiemannSolver::material(double x, double t) const {

}

double ExactRiemannSolver::density(double x, double t) const {

}

Vector3d ExactRiemannSolver::velocity(double x, double t) const {

}

double ExactRiemannSolver::pressure(double x, double t) const {

}

double ExactRiemannSolver::energy(double x, double t) const {

}

}
}