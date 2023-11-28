#include <iostream>

#include <zephyr/geom/polygon.h>

using namespace zephyr::geom;

void plot(Polygon& poly) {
    std::cout << "[";
    for (int i = 0; i < poly.size() - 1; ++i) {
        std::cout << poly[i].x() << ", ";
    }
    std::cout << poly[poly.size() - 1].x() << "], [";
    for (int i = 0; i < poly.size() - 1; ++i) {
        std::cout << poly[i].y() << ", ";
    }
    std::cout << poly[poly.size() - 1].y() << "]\n";
}


// Найти наилучшее отсечение
// a     -- объемная доля in [0, 1]
// n_phi -- угол нормали сечения in [0, п]
// v     -- модуль скорости [0, 0.5]
// v_phi -- угол скорости [п/2, п/2]
double best_alpha(double a, double n_phi, double v, double v_phi) {
}

void plot(double a, double n_phi, double v, double v_phi) {
    Vector3d fn = {1.0, 0.0, 0.0};

    Vector3d c1 = {0.0, 0.0, 0.0};
    Vector3d c2 = {1.0, 0.0, 0.0};
    Vector3d c3 = {1.0, 1.0, 0.0};
    Vector3d c4 = {0.0, 1.0, 0.0};

    PolyQuad cell(c1, c2, c3, c4);


    Vector3d V = {v * std::cos(v_phi), v * std::sin(v_phi), 0.0};

    Line seg1 = {c2, c2 - V};
    Line seg2 = {c3, c3 - V};

    auto poly2 = cell.clip(c2 - (V).dot(fn) * fn, -fn);
    auto poly3 = poly2.clip(seg1.center(), seg1.normal(0.5 * (c2 + c3)));
    auto poly4 = poly3.clip(seg2.center(), seg2.normal(0.5 * (c2 + c3)));

    Vector3d n = {std::cos(n_phi), std::sin(n_phi), 0.0};
    Vector3d p = cell.find_section(n, a);

    auto body = cell.clip(p, n);
    auto move1 = poly4.clip(p, n);

    double F_VOF = move1.area();

    std::cout << "F_VOF: " << F_VOF << "\n";

    double xi = v * std::cos(v_phi);
    double F_min = std::max(0.0, xi - (1.0 - a));
    double F_max = a;

    double a_sig;
    if (F_VOF < F_min) {
        a_sig = 0.0;
    } else if (F_VOF < F_max) {
        a_sig = F_VOF / xi;
    } else {
        a_sig = a;
    }

    double delta = a_sig < a ? (1 - a) / (1 - a_sig) : a / a_sig;

    PolygonS<6> crp_body;
    if (a < a_sig) {
        std::cout << "here1\n";
        crp_body.resize(4);
        crp_body[0] = {1.0 - delta, 0.0, 0.0};
        crp_body[1] = {1.0,0.0, 0.0};
        crp_body[2] = {1.0, a_sig, 0.0};
        crp_body[3] = {1.0 - delta, a_sig, 0.0};
    }
    else {
        std::cout << "here2\n";
        crp_body.resize(6);
        crp_body[0] = {0.0, 0.0, 0.0};
        crp_body[1] = {1.0, 0.0, 0.0};
        crp_body[2] = {1.0, a_sig, 0.0};
        crp_body[3] = {1.0 - delta, a_sig, 0.0};
        crp_body[4] = {1.0 - delta, 1.0, 0.0};
        crp_body[5] = {0.0, 1.0, 0.0};
    }

    //auto move2 = crp_body.clip(p, n);

    plot(poly4);
    plot(body);
    plot(move1);

    plot(crp_body);
}

int main() {


    plot(0.4, 0.4 * M_PI, 0.3, -0.1*M_PI);



}