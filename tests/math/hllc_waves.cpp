#include <zephyr/math/cfd/flux/hllc.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/tests/ivp.h>
#include <zephyr/phys/matter/materials.h>

using namespace zephyr::math;
using namespace zephyr::phys;

class Test1D : public IVP{
    double density(const Vector3d &r) const final{return 0.0;};

    Vector3d velocity(const Vector3d &r) const final{return r;};

    double pressure(const Vector3d &r) const final{return 0.0;};

    double max_time() const final{return 0.0;};
};

class Test_HLLC_waves : public Test1D {
public:
    double gamma_1 = 1.4;

    double rL = 1.0; double rR = 1.0;
    double uL = -2.0; double uR = 2.0;
    double pL = 0.4; double pR = 0.4;

    Test_HLLC_waves() {
        m_materials += IdealGas::create(gamma_1);

        Eos::Ptr eos_L = m_materials[0].eos();
        Eos::Ptr eos_R = m_materials.single() ? eos_L : m_materials[1].eos();

        double eL = eos_L->energy_rP(rL, pL);
        double eR = eos_R->energy_rP(rR, pR);

        smf::PState zL(rL, {uL, 0.0, 0.0}, pL, eL);
        smf::PState zR(rR, {uR, 0.0, 0.0}, pR, eR);

        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(*eos_L, zL, *eos_R, zR);
        std::cout << "S_1L = " << S_1L << std::endl;
        std::cout << "S_1C = " << S_1C << std::endl;
        std::cout << "S_1R = " << S_1R << std::endl;
        std::cout << "Q_s1L = " << Q_s1L << std::endl;
        std::cout << "F_s1L = " << F_s1L << std::endl;
        std::cout << "Q_s1R = " << Q_s1R << std::endl;
        std::cout << "F_s1R = " << F_s1R << std::endl;

        auto[S_2L, S_2C, S_2R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config_u_R(S_1C, *eos_R, zR);
        std::cout << "S_2L = " << S_2L << std::endl;
        std::cout << "S_2C = " << S_2C << std::endl;
        std::cout << "S_2R = " << S_2R << std::endl;
        std::cout << "Q_s2L = " << Q_s2L << std::endl;
        std::cout << "F_s2L = " << F_s2L << std::endl;
        std::cout << "Q_s2R = " << Q_s2R << std::endl;
        std::cout << "F_s2R = " << F_s2R << std::endl;

        auto[S_3L, S_3C, S_3R, Q_s3L, F_s3L, Q_s3R, F_s3R] = HLLC::wave_config_u_L(S_1C, *eos_L, zL);
        std::cout << "S_3L = " << S_3L << std::endl;
        std::cout << "S_3C = " << S_3C << std::endl;
        std::cout << "S_3R = " << S_3R << std::endl;
        std::cout << "Q_s3L = " << Q_s3L << std::endl;
        std::cout << "F_s3L = " << F_s3L << std::endl;
        std::cout << "Q_s3R = " << Q_s3R << std::endl;
        std::cout << "F_s3R = " << F_s3R << std::endl;
    }
};

int main(){
    Test_HLLC_waves test;

    return 0;
}