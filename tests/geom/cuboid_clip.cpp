
#if 0
double volfrac(double Pz, Vector3d n, double a, double b, double c) {
    std::vector<double> kek={std::abs(a * n.x()), std::abs(b * n.y()), std::abs(c * n.z())};
    std::sort(kek.begin(), kek.end());

    double xi = kek[0];
    double eta = kek[1];
    double chi = kek[2];

    std::cout << "  n: " << n << "\n";
    std::cout << "  xi : " << xi << "\n";
    std::cout << "  eta: " << eta << "\n";
    std::cout << "  chi: " << chi << "\n";

    double p = Pz + 0.5 * (xi + eta + chi);

    bool inv = p > 0.5 * (xi + eta + chi);

    if (inv) {
        std::cout << "  Inversed case\n";
        p = xi + eta + chi - p;
    }

    double vol = NAN;
    if (p < xi) {
        std::cout << "  case 1\n";
        vol = std::pow(p, 3) / (6.0 * xi * eta * chi);
    }
    else if (p < eta) {
        std::cout << "  case 2\n";
        vol = (3.0 * p * (p - xi) + xi * xi) / (6.0 * eta * chi);
    }
    else if (p < std::min(xi + eta, chi)) {
        std::cout << "  case 3\n";

        vol = (std::pow(p, 3) - std::pow(p - xi, 3) - std::pow(p - eta, 3)) / (6.0 * xi * eta * chi);
    }
    else {
        if (xi + eta < chi) {
            std::cout << "  case 4.1\n";
            vol = (2.0 * p - xi - eta) / (2.0 * chi);
        } else {
            std::cout << "  case 4.2\n";

            vol = (std::pow(p, 3) - std::pow(p - xi, 3) - std::pow(p - eta, 3) - std::pow(p - chi, 3)) / (6.0 * xi * eta * chi);
        }
    }

    if (inv) {
        vol = 1.0 - vol;
    }

    return vol;

}
#endif

int main() {
    return 0;
}