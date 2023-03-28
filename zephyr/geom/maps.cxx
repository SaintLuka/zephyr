#include <vector>
#include <functional>
#include <iomanip>

#include <zephyr/geom/maps.h>

namespace zephyr { namespace geom {


Mapping1D::Mapping1D(const ShortList1D& vs)
    : v1(vs[0]), v2(vs[1]) {
    vc = 0.5 * (v1 + v2);
}

Mapping1D::Mapping1D(const LargeList1D& vs)
    : v1(vs[0]), vc(vs[1]), v2(vs[2]) {

}

Vector3d Mapping1D::get(const Vector3d& v1, const Vector3d& vc, const Vector3d& v2, double xi) {
    return vc + 0.5 * (v2 - v1) * xi + 0.5 * (v2 - 2 * vc + v1) * xi * xi;
}

Vector3d Mapping1D::tangent(const Vector3d& v1, const Vector3d& vc, const Vector3d& v2, double xi) {
    return 0.5 * (v2 - v1) + (v2 - 2.0 * vc + v1) * xi;
}

Vector3d Mapping1D::operator()(double xi) const {
    return get(v1, vc, v2, xi);
}

Vector3d Mapping1D::normal(double xi) const {
    Vector3d a = tangent(v1, vc, v2, xi);
    if (std::abs(a.z()) <= std::abs(a.x()) && std::abs(a.z()) <= std::abs(a.y())) {
        return Vector3d(+a.y(), -a.x(), 0.0).normalized();
    }
    else if (std::abs(a.y()) <= std::abs(a.x()) && std::abs(a.y()) <= std::abs(a.z())) {
        return Vector3d(-a.z(), 0.0, +a.x()).normalized();
    }
    else {
        return Vector3d(0.0, +a.z(), -a.y()).normalized();
    }
}

double  Mapping1D::Jacobian(double xi) const {
    return tangent(v1, vc, v2, xi).norm();
}

Mapping2D::Mapping2D(const ShortList2D& vs) :
    SB(ShortList1D({vs[0], vs[1]})),
    SH(ShortList1D({0.5*(vs[0] + vs[2]), 0.5*(vs[1] + vs[3])})),
    ST(ShortList1D({vs[2], vs[3]})) {

}

Mapping2D::Mapping2D(const LargeList2D& vs) :
        SB(LargeList1D({vs[0], vs[1], vs[2]})),
        SH(LargeList1D({vs[3], vs[4], vs[5]})),
        ST(LargeList1D({vs[6], vs[7], vs[8]})) {

}

Vector3d Mapping2D::operator()(double xi, double eta) const {
    return Mapping1D::get(SB(xi), SH(xi), ST(xi), eta);
}

Vector3d Mapping2D::normal(double xi, double eta) const {
    Mapping1D SL({SB.v1, SH.v1, ST.v1});
    Mapping1D SV({SB.vc, SH.vc, ST.vc});
    Mapping1D SR({SB.v2, SH.v2, ST.v2});

    Vector3d W_x = Mapping1D::tangent(SL(eta), SV(eta), SR(eta), xi);
    Vector3d W_e = Mapping1D::tangent(SB(xi), SH(xi), ST(xi), eta);

    return W_x.cross(W_e).normalized();
}

double Mapping2D::Jacobian(double xi, double eta) const {
    Mapping1D SL({SB.v1, SH.v1, ST.v1});
    Mapping1D SV({SB.vc, SH.vc, ST.vc});
    Mapping1D SR({SB.v2, SH.v2, ST.v2});

    Vector3d W_x = Mapping1D::tangent(SL(eta), SV(eta), SR(eta), xi);
    Vector3d W_e = Mapping1D::tangent(SB(xi), SH(xi), ST(xi), eta);

    return W_x.cross(W_e).norm();
}

Mapping3D::Mapping3D(const ShortList3D& vs)
    : vs(vs) { }

Vector3d Mapping3D::operator()(double xi, double eta, double chi) {
    Vector3d res = {0.0, 0.0, 0.0};
    res += ((1 - xi) * (1 - eta) * (1 - chi)) * vs[0];
    res += ((1 + xi) * (1 - eta) * (1 - chi)) * vs[1];
    res += ((1 - xi) * (1 + eta) * (1 - chi)) * vs[2];
    res += ((1 + xi) * (1 + eta) * (1 - chi)) * vs[3];
    res += ((1 - xi) * (1 - eta) * (1 + chi)) * vs[4];
    res += ((1 + xi) * (1 - eta) * (1 + chi)) * vs[5];
    res += ((1 - xi) * (1 + eta) * (1 + chi)) * vs[6];
    res += ((1 + xi) * (1 + eta) * (1 + chi)) * vs[7];
    res /= 8.0;
    return res;
}

double Mapping3D::Jacobian(double xi, double eta, double chi) {
    // Производные по направлениям (xi, eta, chi)
    // Формулы актуальны в силу линейности отображения по каждой из переменных
    Vector3d N1 = {0.0, 0.0, 0.0};
    N1 += ((1 - eta) * (1 - chi)) * (vs[1] - vs[0]);
    N1 += ((1 + eta) * (1 - chi)) * (vs[3] - vs[2]);
    N1 += ((1 + eta) * (1 + chi)) * (vs[7] - vs[6]);
    N1 += ((1 - eta) * (1 + chi)) * (vs[5] - vs[4]);

    Vector3d N2 = {0.0, 0.0, 0.0};
    N2 += ((1 - xi) * (1 - chi)) * (vs[2] - vs[0]);
    N2 += ((1 + xi) * (1 - chi)) * (vs[3] - vs[1]);
    N2 += ((1 - xi) * (1 + chi)) * (vs[6] - vs[4]);
    N2 += ((1 + xi) * (1 + chi)) * (vs[7] - vs[5]);

    Vector3d N3 = {0.0, 0.0, 0.0};
    N3 += ((1 - xi) * (1 - eta)) * (vs[4] - vs[0]);
    N3 += ((1 + xi) * (1 - eta)) * (vs[5] - vs[1]);
    N3 += ((1 - xi) * (1 + eta)) * (vs[6] - vs[2]);
    N3 += ((1 + xi) * (1 + eta)) * (vs[7] - vs[3]);

    return N1.dot(N2.cross(N3)) / (8.0*8.0*8.0);
}

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a + (b - a) * i / (n - 1.0);
    }
    return res;
}

// Одномерный интеграл
template <class V>
V simpson(const std::vector<double>& t, const std::function<V(double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(t.size()) - 1; ++i) {
        res += (t[i + 1] - t[i]) * (f(t[i]) + 4.0 * f(0.5 * (t[i] + t[i + 1])) + f(t[i + 1]));
    }
    return res / 6.0;
}

// Двумерный интеграл
template <class V>
V simpson(
        const std::vector<double>& xi,
        const std::vector<double>& eta,
        const std::function<V(double, double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(xi.size()) - 1; ++i) {
        double xi1 = xi[i];
        double xi2 = xi[i + 1];
        double xic = 0.5*(xi1 + xi2);
        for (int j = 0; j < int(eta.size()) - 1; ++j) {
            double eta1 = eta[j];
            double eta2 = eta[j + 1];
            double etac = 0.5 * (eta1 + eta2);

            res += (xi2 - xi1) * (eta2 - eta1) * (
                    f(xi1, eta1) + 4 * f(xi1, etac) + f(xi1, eta2) +
                    4 * f(xic, eta1) + 16 * f(xic, etac) + 4 * f(xic, eta2) +
                    f(xi2, eta1) + 4 * f(xi2, etac) + f(xi2, eta2)
            );
        }
    }
    return res / 36.0;
}

// Трехмерный интеграл
template <class V>
V simpson(
        const std::vector<double>& xi,
        const std::vector<double>& eta,
        const std::vector<double>& chi,
        const std::function<V(double, double, double)>& f, V zero) {
    V res = zero;
    for (int i = 0; i < int(xi.size()) - 1; ++i) {
        double xi1 = xi[i];
        double xi2 = xi[i + 1];
        double xic = 0.5 * (xi1 + xi2);
        for (int j = 0; j < int(eta.size()) - 1; ++j) {
            double eta1 = eta[j];
            double eta2 = eta[j + 1];
            double etac = 0.5 * (eta1 + eta2);

            for (int k = 0; k < int(chi.size()) - 1; ++k) {
                double chi1 = chi[k];
                double chi2 = chi[k + 1];
                double chic = 0.5 * (chi1 + chi2);


                res += (xi2 - xi1) * (eta2 - eta1) * (chi2 - chi1) * (
                        64.0 * f(xic, etac, chic) +
                        16.0 * (
                                f(xic, etac, chi1) + f(xic, etac, chi2) +
                                f(xic, eta1, chic) + f(xic, eta2, chic) +
                                f(xi1, etac, chic) + f(xi2, etac, chic)
                        ) +
                        4.0 * (
                                f(xic, eta1, chi1) + f(xic, eta1, chi2) +
                                f(xic, eta2, chi1) + f(xic, eta2, chi2) +
                                f(xi1, etac, chi1) + f(xi1, etac, chi2) +
                                f(xi2, etac, chi1) + f(xi2, etac, chi2) +
                                f(xi1, eta1, chic) + f(xi1, eta2, chic) +
                                f(xi2, eta1, chic) + f(xi2, eta2, chic)
                        ) +
                        f(xi1, eta1, chi1) + f(xi1, eta1, chi2) +
                        f(xi1, eta2, chi1) + f(xi1, eta2, chi2) +
                        f(xi2, eta1, chi1) + f(xi2, eta1, chi2) +
                        f(xi2, eta2, chi1) + f(xi2, eta2, chi2)

                );
            }
        }
    }
    return res / (6.0 * 6.0 * 6.0);
}

void checkout() {

    // Основной тестовый двумерный набор
    LargeList2D ML_2D = {
            Vector3d(0.1, 1.0, 0.0), Vector3d(0.7, 1.1, 0.0), Vector3d(1.3, 1.1, 0.0),
            Vector3d(0.2, 1.5, 0.0), Vector3d(0.8, 1.7, 0.0), Vector3d(1.4, 1.8, 0.0),
            Vector3d(0.0, 2.0, 0.0), Vector3d(0.5, 2.4, 0.0), Vector3d(1.2, 2.6, 0.0)
    };

    // Основной тестовый трехмерный набор
    ShortList3D ML_3D = {
            Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, 0.2, 0.1),
            Vector3d(0.0, 1.0, 0.3), Vector3d(1.0, 1.2, 0.2),
            Vector3d(0.0, 0.1, 1.0), Vector3d(1.5, 0.5, 1.0),
            Vector3d(0.4, 1.0, 1.0), Vector3d(1.6, 1.5, 1.3),
    };

    // Характерная длина, площадь, объем
    double L1 = (ML_2D[0] - ML_2D[8]).norm();
    double L2 = L1 * L1;
    double L3 = L1 * L1 * L1;

    // Центр ячейки
    Vector3d C = ML_2D[4];

    std::cout << std::setprecision(6);

    {
        std::cout << "Длина и нормаль простого отрезка\n";
        ShortList1D vs = {ML_2D[0], ML_2D[2]};
        Mapping1D map1D(vs);

        // int normal(t) * dl(t) = int tangent(t) normal(t) dt = length() * normal()
        Vector3d res1 = length(vs) * normal(vs, C);

        std::function<Vector3d(double)> func =
                [&map1D](double xi) -> Vector3d {
                    return map1D.Jacobian(xi) * map1D.normal(xi);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res2 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Длина и нормаль криволинейного отрезка\n";

        LargeList1D vs = {ML_2D[0], ML_2D[1], ML_2D[2]};
        Mapping1D map1D(vs);

        // int normal(t) * dl(t) = int tangent(t) normal(t) dt = length() * normal()
        Vector3d res1 = length(vs) * normal(vs, C);

        std::function<Vector3d(double)> func =
                [&map1D](double xi) -> Vector3d {
                    return map1D.Jacobian(xi) * map1D.normal(xi);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Площадь простой двумерной ячейки\n";

        ShortList2D vs = {ML_2D[0], ML_2D[2], ML_2D[6], ML_2D[8]};
        Mapping2D map2D(vs);

        double res1 = area(vs);

        std::function<double(double, double)> func =
                [&map2D](double xi, double eta) -> double {
                    return map2D.Jacobian(xi, eta);
                };

        double res2 = simpson(
                linspace(-1.0, 1.0, 51),
                linspace(-1.0, 1.0, 51),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Площадь и нормаль простой изогнутой двумерной ячейки\n";

        ShortList2D vs = {ML_2D[0], ML_2D[2], ML_2D[6], ML_2D[8]};
        vs[2].z() += 0.1;
        vs[3].z() += 0.2;
        Vector3d C2 = {C.x(), C.y(), C.z() - 0.3};

        Mapping2D map2D(vs);

        Vector3d res1 = area(vs) * normal(vs, C2);

        std::function<Vector3d(double, double)> func =
                [&map2D](double xi, double eta) -> Vector3d {
                    return map2D.Jacobian(xi, eta) * map2D.normal(xi, eta);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 51),
                linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0});

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L2 << "\n\n";
    }

    {
        std::cout << "Площадь криволинейной двумерной ячейки\n";

        LargeList2D& vs = ML_2D;
        Mapping2D map2D(vs);

        double res1 = area(vs);

        std::function<double(double, double)> func =
                [&map2D](double xi, double eta) -> double {
                    return map2D.Jacobian(xi, eta);
                };

        double res2 = simpson(
                linspace(-1.0, 1.0, 101),
                linspace(-1.0, 1.0, 101),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем обычной двумерной ячейки в осесимметричной постновке\n";

        ShortList2D vs = {ML_2D[0], ML_2D[2], ML_2D[6], ML_2D[8]};
        Mapping2D map2D(vs);

        double res1 = volume_as(vs);

        std::function<double(double, double)> func =
                [&map2D](double xi, double eta) -> double {
                    return map2D(xi, eta).y() * map2D.Jacobian(xi, eta);
                };

        double res2 = simpson(
                linspace(-1.0, 1.0, 51),
                linspace(-1.0, 1.0, 51),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем криволинейной ячейки в осесимметричной постновке\n";

        LargeList2D& vs = ML_2D;
        Mapping2D map2D(vs);

        double res1 = volume_as(vs);


        std::function<double(double, double)> func =
                [&map2D](double xi, double eta) -> double {
                    return map2D(xi, eta).y() * map2D.Jacobian(xi, eta);
                };

        double res2 = simpson(
                linspace(-1.0, 1.0, 501),
                linspace(-1.0, 1.0, 501),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L2 << "\n\n";
    }

    {
        std::cout << "Объем трехмерной ячейки\n";

        ShortList3D& vs = ML_3D;
        Mapping3D map3D(vs);

        double res1 = volume(vs);

        std::function<double(double, double, double)> func =
                [&map3D](double xi, double eta, double chi) -> double {
                    return map3D.Jacobian(xi, eta, chi);
                };

        double res2 = simpson(
                linspace(-1.0, 1.0, 41),
                linspace(-1.0, 1.0, 41),
                linspace(-1.0, 1.0, 41),
                func, 0.0);

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << std::abs(res1 - res2) / L3 << "\n\n";
    }

    {
        std::cout << "Барицентр простой двумерной ячейки\n";

        ShortList2D vs = {ML_2D[0], ML_2D[2], ML_2D[6], ML_2D[8]};
        Mapping2D map2D(vs);

        double S = area(vs);

        Vector3d res1 = centroid(vs, S);

        std::function<Vector3d(double, double)> func =
                [&map2D](double xi, double eta) -> Vector3d {
                    return map2D(xi, eta) * map2D.Jacobian(xi, eta);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 51),
                linspace(-1.0, 1.0, 51),
                func, {0.0, 0.0, 0.0}) / S;

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Барицентр криволинейной двумерной ячейки\n";

        LargeList2D& vs = ML_2D;
        Mapping2D map2D(vs);

        double S = area(vs);

        Vector3d res1 = centroid(vs, S);

        std::function<Vector3d(double, double)> func =
                [&map2D](double xi, double eta) -> Vector3d {
                    return map2D(xi, eta) * map2D.Jacobian(xi, eta);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 501),
                linspace(-1.0, 1.0, 501),
                func, {0.0, 0.0, 0.0}) / S;

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }

    {
        std::cout << "Барицентр трехмерной ячейки\n";

        ShortList3D& vs = ML_3D;
        Mapping3D map3D(vs);

        double V = volume(vs);
        Vector3d res1 = centroid(vs, V);

        std::function<Vector3d(double, double, double)> func =
                [&map3D](double xi, double eta, double chi) -> Vector3d {
                    return map3D(xi, eta, chi) * map3D.Jacobian(xi, eta, chi);
                };

        Vector3d res2 = simpson(
                linspace(-1.0, 1.0, 11),
                linspace(-1.0, 1.0, 11),
                linspace(-1.0, 1.0, 11),
                func, {0.0, 0.0, 0.0}) / V;

        std::cout << "    res (formula): " << res1 << "\n";
        std::cout << "    res (numeric): " << res2 << "\n";
        std::cout << "    error: " << (res1 - res2).norm() / L1 << "\n\n";
    }
}

} // geom
} // zephyr