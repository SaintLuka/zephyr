#include <zephyr/phys/tests/test_2D.h>

namespace zephyr::phys {

// ============================================================================
//                         Riemann 2D Tests
// ============================================================================

Riemann2D::Riemann2D(int test_case) {
    m_materials += IdealGas::create(1.4);

    switch (test_case){
        case 1: {
            r1 = 1; r2 = 0.5197;  r3 = 0.1072;  r4 = 0.2579;
            u1 = 0; u2 = -0.7259; u3 = -0.7259; u4 = 0;
            v1 = 0; v2 = 0;       v3 = -1.4045; v4 = -1.4045;
            p1 = 1; p2 = 0.4;     p3 = 0.0439;  p4 = 0.15;
            break;
        }
        case 2:
            r1 = 1; r2 = 0.4;     r3 = 1;       r4 = 0.5197;
            u1 = 0; u2 = -0.7259; u3 = -0.7259; u4 = 0;
            v1 = 0; v2 = 0;       v3 = -0.7259; v4 = -0.7259;
            p1 = 1; p2 = 0.5197;  p3 = 1;       p4 = 0.4;
            break;
        case 3:
            r1 = 1.5; r2 = 0.5323; r3 = 0.138; r4 = 0.5323;
            u1 = 0;   u2= 1.206;   u3 = 1.206; u4 = 0;
            v1 = 0;   v2 = 0;      v3 = 1.206; v4 = 1.206;
            p1 = 1.5; p2 = 0.3;    p3 = 0.029; p4 = 0.3;
            break;
        case 4:
            r1 = 1.1; r2 = 0.5065; r3 = 1.1;    r4 = 0.5065;
            u1 = 0;   u2= 0.8939;  u3 = 0.8939; u4 = 0;
            v1 = 0;   v2 = 0;      v3 = 0.8939; v4 = 0.8939;
            p1 = 1.1; p2 = 0.35;   p3 = 1.1;    p4 = 0.35;
            break;
        case 5:
            r1 = 1;     r2 = 2;     r3 = 1;    r4 = 3;
            u1 = -0.75; u2 = -0.75; u3 = 0.75; u4 = 0.75;
            v1 = -0.5;  v2 = 0;     v3 = 0.5;  v4 = -0.5;
            p1 = 1;     p2 = 1;     p3 = 1;    p4 = 1;
            break;
        case 6:
            r1 = 1;     r2 = 2;     r3 = 1;     r4 = 3;
            u1 = 0.75;  u2 = 0.75;  u3 = -0.75; u4 = -0.75;
            v1 = -0.5;  v2 = 0.5;   v3 = 0.5;   v4 = -0.5;
            p1 = 1;     p2 = 1;     p3 = 1;     p4 = 1;
            finish = 0.3;
            break;
        case 7:
            r1 = 0.5197; r2 = 1.0;     r3 = 0.8; r4 = 1.0;
            p1 = 0.4;    p2 = 1.0;     p3 = 1.0; p4 = 1.0;
            u1 = 0.1;    u2 = -0.6259; u3 = 0.1; u4 = 0.1;
            v1 = 0.1;    v2 = 0.1;     v3 = 0.1; v4 = -0.6259;
            finish = 0.25;
            break;
        case 8:
            r1 = 1.0; r2 = 2.0;  r3 = 1.039;   r4 = 0.5197;
            p1 = 1.0; p2 = 1.0;  p3 = 0.4;     p4 = 0.4;
            u1 = 0.0; u2 = 0.0;  u3 = 0.0;     u4 = 0.0;
            v1 = 0.3; v2 = -0.3; v3 = -0.8133; v4 = -0.259;
            break;
        default:
            throw std::runtime_error("Unknown Test");
    }

    y_jump = 0.5;
    x_jump = 0.5;
}

double Riemann2D::density(const Vector3d &vec) const {
    if (vec.x() >= x_jump && vec.y() >= y_jump)
        return r1;
    if (vec.x() <= x_jump && vec.y() >= y_jump)
        return r2;
    if (vec.x() <= x_jump && vec.y() <= y_jump)
        return r3;
    if (vec.x() >= x_jump && vec.y() <= y_jump)
        return r4;
    return NAN;
}

Vector3d Riemann2D::velocity(const Vector3d &vec) const {
    if (vec.x() >= x_jump && vec.y() >= y_jump)
        return {u1, v1, 0.0};
    if (vec.x() <= x_jump && vec.y() >= y_jump)
        return {u2, v2, 0.0};
    if (vec.x() <= x_jump && vec.y() <= y_jump)
        return {u3, v3, 0.0};
    if (vec.x() >= x_jump && vec.y() <= y_jump)
        return {u4, v4, 0.0};
    return {NAN, NAN, NAN};
}

double Riemann2D::pressure(const Vector3d &vec) const {
    if (vec.x() >= x_jump && vec.y() >= y_jump)
        return p1;
    if (vec.x() <= x_jump && vec.y() >= y_jump)
        return p2;
    if (vec.x() <= x_jump && vec.y() <= y_jump)
        return p3;
    if (vec.x() >= x_jump && vec.y() <= y_jump)
        return p4;
    return NAN;
}

// ============================================================================
//                      Richtmyer - Meshkov Instability
// ============================================================================

RichtmyerMeshkov::RichtmyerMeshkov(double Ms, bool multimat) {
    auto air = IdealGas::create("Air");
    m_materials += air;
    if (!multimat) {
        m_materials += air;
    }
    else {
        m_materials += IdealGas::create("SF6");
    }

    double gamma = air->cast<IdealGas>().gamma;

    // Невозмущенная область
    p0 = 1.0;
    r0 = 0.01;
    u0 = 0.0;

    // Перед фронтом УВ
    pR = 1.0;
    rR = 1.0;
    uR = 0.0;

    // За фронтом УВ
    pL = pR * (2 * gamma * Ms * Ms - gamma + 1) / (gamma + 1);
    rL = rR * (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms);
    uL = 2.0 / Ms * std::sqrt(gamma * pR / rR) * (Ms * Ms - 1) / (gamma + 1);

    x_jump = 0.5;
    x_contact = 1.0;
    amplitude = 0.02;
    y_width = 1.0;
    n_peaks = 3;

    finish = 3.0;
}

int RichtmyerMeshkov::region(const Vector3d& r) const {
    if (r.x() < x_jump) {
        return 0;
    }
    // длина волны возмущения
    double L = y_width / n_peaks;
    // s in [0, 1]
    double s = r.y() / L - std::floor(r.y() / L);
    // пила
    double x = x_contact + amplitude * (std::abs(s - 0.5) - 0.5);
    return r.x() < x ? 1 : 2;
}

int RichtmyerMeshkov::index(const Vector3d &r) const {
    if (m_materials.single()) {
        return 0;
    }
    switch (region(r)) {
        case 0:  return 0;
        case 1:  return 0;
        default: return 1;
    }
}

double RichtmyerMeshkov::density(const Vector3d &r) const {
    switch (region(r)) {
        case 0:  return rL;
        case 1:  return rR;
        default: return r0;
    }
}

Vector3d RichtmyerMeshkov::velocity(const Vector3d &r) const {
    switch (region(r)) {
        case 0:  return {uL, 0.0, 0.0};
        case 1:  return {uR, 0.0, 0.0};
        default: return {u0, 0.0, 0.0};
    }
}

double RichtmyerMeshkov::pressure(const Vector3d &r) const {
    switch (region(r)) {
        case 0:  return pL;
        case 1:  return pR;
        default: return p0;
    }
}

} // namespace zephyr::phys