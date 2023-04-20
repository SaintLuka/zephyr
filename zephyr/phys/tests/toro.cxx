#include <zephyr/phys/tests/toro.h>

namespace zephyr { namespace phys {

ToroTest::ToroTest(int num)
        : eos(1.4) {

    switch (num) {
        case 1:
            rL = 1.0;
            uL = 0.75;
            pL = 1.0;
            rR = 0.125;
            uR = 0.0;
            pR = 0.1;
            x_jump = 0.3;
            finish = 0.2;
            break;
        case 2:
            rL = 1.0;
            uL = -2.0;
            pL = 0.4;
            rR = 1.0;
            uR = 2.0;
            pR = 0.4;
            x_jump = 0.5;
            finish = 0.15;
            break;
        case 3:
            rL = 1.0;
            uL = 0.0;
            pL = 1000.0;
            rR = 1.0;
            uR = 0.0;
            pR = 0.01;
            x_jump = 0.5;
            finish = 0.012;
            break;
        case 4:
            rL = 5.99924;
            uL = 19.5975;
            pL = 460.894;
            rR = 5.99242;
            uR = -6.19633;
            pR = 46.0950;
            x_jump = 0.4;
            finish = 0.035;
            break;
        case 5:
            rL = 1.0;
            uL = -19.59745;
            pL = 1000.0;
            rR = 1.0;
            uR = -19.59745;
            pR = 0.01;
            x_jump = 0.8;
            finish = 0.012;
            break;
        case 6:
            rL = 1.4;
            uL = 0.0;
            pL = 1.0;
            rR = 1.0;
            uR = 0.0;
            pR = 1.0;
            x_jump = 0.5;
            finish = 2.0;
            break;
        case 7:
            rL = 1.4;
            uL = 0.1;
            pL = 1.0;
            rR = 1.0;
            uR = 0.1;
            pR = 1.0;
            x_jump = 0.5;
            finish = 2.0;
            break;
        default:
            throw std::runtime_error("Unknown Toro test (num > 7)");
    }

    eL = eos.energy_rp(rL, pL);
    eR = eos.energy_rp(rR, pR);
}

void ToroTest::inverse() {
    std::swap(rL, rR);
    std::swap(uL, uR);
    std::swap(pL, pR);
    std::swap(eL, eR);
    uL *= -1.0;
    uR *= -1.0;
}

} // phys
} // zephyr