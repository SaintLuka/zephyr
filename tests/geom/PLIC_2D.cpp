// PLIC реконструкция границы на двумерных декартовых сетках

#include <iomanip>
#include <atomic>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/sections.h>

#include <zephyr/geom/geom.h>
#include <zephyr/math/funcs.h>

using namespace zephyr;
using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Rectangle;

Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

// Нормаль для теста с плоскостью
static Vector3d some_n = (Vector3d{0.24, 0.13, 0.0}).normalized();

// Характеристическая функция (функция-индикатор)
using InFunction = std::function<bool(const Vector3d &)>;

// Пространственная функция
using SpFunction = std::function<double(const Vector3d &)>;

// 1. Гладкая замкнутая кривая
InFunction smooth_func = [](const Vector3d &v) -> bool {
    const double R1 = 0.3;
    const double e1 = 0.2;
    const double R2 = 0.7;
    const double e2 = 0.3;
    const double k = 5.0;

    double phi = std::atan2(v.y(), v.x());
    double r = v.norm();
    double r1 = R1 * (1.0 + e1 * std::sin(k * phi));
    double r2 = R2 * (1.0 + e2 * std::sin(k * phi));
    return r1 < r && r < r2;
};

// 2. Некоторая угловатая функция
InFunction angle_func = [](const Vector3d &v) -> bool {
    static bool initialized = false;
    static Polygon poly(10);
    if (!initialized) {
        double R = 0.8;
        double r = 0.382 * R;
        for (int i = 0; i < 10; i += 2) {
            Vector3d v1 = to_cartesian(R, M_PI*(i + 0)/5.0);
            Vector3d v2 = to_cartesian(r, M_PI*(i + 1)/5.0);
            poly.set(i, v1);
            poly.set(i + 1, v2);
        }
        initialized = true;
    }
    return poly.inside(v);
};

// 3. Функция с диффузной границей
SpFunction diffuse_func = [](const Vector3d &v) -> double {
    const double R1 = 0.3;
    const double e1 = 0.2;
    const double R2 = 0.7;
    const double e2 = 0.3;
    const double k = 5.0;

    double phi = std::atan2(v.y(), v.x());
    double r = v.norm();
    double r1 = R1 * (1.0 + e1 * std::sin(k * phi));
    double r2 = R2 * (1.0 + e2 * std::sin(k * phi));
    double h = 0.08;
    return 0.25 * (1 + math::sign_p(r - r1, h)) * (1 + math::sign_p(r2 - r, h));
};

// 4. Зашумленная функция (случайный шум)
InFunction noise_func = [](const Vector3d &v) -> bool {
    return false;
};

// Объемная доля
static Storable<double> a;

// Центральные разности
static Storable<Vector3d> n1;
static Storable<double> p1;
static Storable<double> e1;

// Формула Youngs
static Storable<Vector3d> n2;
static Storable<double> p2;
static Storable<double> e2;

// Схема ELVIRA
static Storable<Vector3d> n3;
static Storable<double> p3;
static Storable<double> e3;

// Моя формула
static Storable<Vector3d> n4;
static Storable<double> p4;
static Storable<double> e4;

inline bool mixed(EuCell& cell) {
    return 0.0 < cell[a] && cell[a] < 1.0;
}

class Stencil {
    std::array<std::array<double, 3>, 3> arr;

public:
    Stencil(EuCell& cell) {
        for (int i: {-1, 0, 1}) {
            for (int j: {-1, 0, 1}) {
                C(i, j) = cell.neib(i, j)[a];
            }
        }
    }

    double& C(int i, int j)       { return arr[i + 1][j + 1]; }
    double  C(int i, int j) const { return arr[i + 1][j + 1]; }

    double& operator()(int i, int j)       { return C(i, j); }
    double  operator()(int i, int j) const { return C(i, j); }

    Vector3d Youngs(double hx, double hy) const {
        Vector3d norm = Vector3d::Zero();
        norm.x() = ((C(-1, +1) + 2 * C(-1, 0) + C(-1, -1)) - (C(+1, +1) + 2 * C(+1, 0) + C(+1, -1))) / hx;
        norm.y() = ((C(+1, -1) + 2 * C(0, -1) + C(-1, -1)) - (C(+1, +1) + 2 * C(0, +1) + C(-1, +1))) / hy;
        return norm.normalized();
    }

    Vector3d ELVIRA(double hx, double hy) const {
        Vector3d YN = Youngs(hx, hy);

        // Функции столбцов y(x)
        double y1 = (C(-1, -1) + C(-1, 0) + C(-1, 1)) * hy;
        double y2 = (C( 0, -1) + C( 0, 0) + C( 0, 1)) * hy;
        double y3 = (C(+1, -1) + C(+1, 0) + C(+1, 1)) * hy;

        double sgn_ny = math::sign(YN.y());

        double n_xl = (y2 - y1) / hx;
        double n_xc = (y3 - y1) / (2 * hx);
        double n_xr = (y3 - y2) / hx;

        // Функции строк x(y)
        double x1 = (C(-1, -1) + C(0, -1) + C(+1, -1)) * hx;
        double x2 = (C(-1,  0) + C(0,  0) + C(+1,  0)) * hx;
        double x3 = (C(-1, +1) + C(0, +1) + C(+1, +1)) * hx;

        double sgn_nx = math::sign(YN.x());

        double n_yb = (x2 - x1) / hy;
        double n_yc = (x3 - x1) / (2 * hy);
        double n_yt = (x3 - x2) / hy;

        // Нормали для тестов
        std::array ns = {
            Vector3d{-n_xl, sgn_ny, 0.0}.normalized(),
            Vector3d{-n_xc, sgn_ny, 0.0}.normalized(),
            Vector3d{-n_xr, sgn_ny, 0.0}.normalized(),
            Vector3d{sgn_nx, -n_yb, 0.0}.normalized(),
            Vector3d{sgn_nx, -n_yc, 0.0}.normalized(),
            Vector3d{sgn_nx, -n_yt, 0.0}.normalized(),
        };

        std::array<double, 6> errors;

        for (int i = 0; i < 6; ++i) {
            Vector3d n = ns[i];
            double p = quad_find_section(C(0, 0), n, hx, hy);

            errors[i] = 0.0;            
            errors[i] += std::pow(C(-1,  0) - quad_volume_fraction(p + hx * n.x(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1,  0) - quad_volume_fraction(p - hx * n.x(), n, hx, hy), 2);
            errors[i] += std::pow(C( 0, -1) - quad_volume_fraction(p + hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C( 0, +1) - quad_volume_fraction(p - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(-1, -1) - quad_volume_fraction(p + hx * n.x() + hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1, +1) - quad_volume_fraction(p - hx * n.x() - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(-1, +1) - quad_volume_fraction(p + hx * n.x() - hy * n.y(), n, hx, hy), 2);
            errors[i] += std::pow(C(+1, -1) - quad_volume_fraction(p - hx * n.x() + hy * n.y(), n, hx, hy), 2);
        }

        int err_min = std::ranges::min_element(errors) - errors.begin();
        return ns[err_min];
    }
};

void make_interface(EuMesh& mesh) {
    mesh.for_each([](EuCell& cell) {
        if (!mixed(cell)) {
            for (auto n: {n1, n2, n3, n4}) cell[n] = Vector3d::Zero();
            for (auto p: {p1, p2, p3, p4}) cell[p] = 0.0;
            for (auto e: {e1, e2, e3, e4}) cell[e] = 0.0;
            return;
        }

        // Размеры ячейки
        double hx = cell.hx();
        double hy = cell.hy();

        // Простая производная
        Vector3d grad = Vector3d::Zero();
        for (auto face: cell.faces()) {
            double a_f = 0.5 * (cell[a] + face.neib(a));
            Vector3d S = face.area() * face.normal();
            grad += a_f * S;
        }
        cell[n1] = -grad.normalized();

        Stencil C(cell);
        cell[n2] = C.Youngs(hx, hy);
        cell[n3] = C.ELVIRA(hx, hy);

        // Моя формула
        grad = Vector3d::Zero();
        for (auto face: cell.faces()) {
            double a_f = face_fraction(cell[a], face.neib(a));
            Vector3d S = face.area() * face.normal();
            grad += a_f * S;
        }
        cell[n4] = -grad.normalized();

        // Пересчитываем позиции прямых
        cell[p1] = quad_find_section(cell[a], cell[n1], hx, hy);
        cell[p2] = quad_find_section(cell[a], cell[n2], hx, hy);
        cell[p3] = quad_find_section(cell[a], cell[n3], hx, hy);
        cell[p4] = quad_find_section(cell[a], cell[n4], hx, hy);
    });
}

void calc_errors(EuMesh& mesh, InFunction func, int nx) {
    std::atomic err1{0.0};
    std::atomic err2{0.0};
    std::atomic err3{0.0};
    std::atomic err4{0.0};

    mesh.for_each([func, nx, &err1, &err2, &err3, &err4](EuCell& cell) {
        // Нулевые погрешности
        cell[e1] = cell[e2] = cell[e3] = cell[e4] = 0.0;

        if (!mixed(cell)) {
            return;
        }
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                return;
            }
        }

        double hx = cell.hx();
        double hy = cell.hy();

        std::array counter = {0, 0, 0, 0};
        for (int i = 0; i < nx; ++i) {
            double x = hx * ((i + 0.5) / nx - 0.5);
            for (int j = 0; j < nx; ++j) {
                double y = hy * ((j + 0.5) / nx - 0.5);

                // Точка относительно центра ячейки
                Vector3d r = {x, y, 0.0};

                bool inside = func(cell.center() + r);
                if (inside != (r.dot(cell[n1]) < cell[p1])) { ++counter[0]; }
                if (inside != (r.dot(cell[n2]) < cell[p2])) { ++counter[1]; }
                if (inside != (r.dot(cell[n3]) < cell[p3])) { ++counter[2]; }
                if (inside != (r.dot(cell[n4]) < cell[p4])) { ++counter[3]; }
            }
        }

        // Интегральная метрика L1, как у Aulisa
        double xi = cell.volume() / (nx * nx);
        cell[e1] = xi * counter[0];
        cell[e2] = xi * counter[1];
        cell[e3] = xi * counter[2];
        cell[e4] = xi * counter[3];

        err1 += cell[e1];
        err2 += cell[e2];
        err3 += cell[e3];
        err4 += cell[e4];
    });

    double hx{NAN}, hy{NAN};
    for (auto cell: mesh) {
        if (mixed(cell)) {
            hx = cell.hx();
            hy = cell.hy();
            break;
        }
    }

    std::cout << "hx, hy               Central      Youngs      ELVIRA       CSIR\n";

    std::cout << std::setprecision(5) << std::fixed;
    std::cout << hx << ", " << hy << ":";

    std::cout << std::setprecision(2) << std::scientific;

    std::cout << std::setw(12) << err1;
    std::cout << std::setw(12) << err2;
    std::cout << std::setw(12) << err3;
    std::cout << std::setw(12) << err4 << "\n";
}

EuMesh body(EuMesh& mesh, Storable<double> p, Storable<Vector3d> n) {
    EuMesh clipped(2, false);
    for (auto cell: mesh) {
        if (cell[a] <= 0.0 || (cell[a] < 0.5 && cell[n].isZero())) {
            continue;
        }
        if (cell[a] >= 1.0 || (cell[a] > 0.5 && cell[n].isZero())) {
            clipped.push_back(cell.polygon());
            continue;
        }

        Vector3d point = cell.center() + cell[p] * cell[n];
        auto poly = cell.polygon().clip(point, cell[n]);
        clipped.push_back(poly);
    }
    return clipped;
}

void save_mesh(EuMesh& mesh) {
    Variables vars = {"level", "flag"};
    vars.append("a", a);
    vars.append("n1", n1);
    vars.append("n2", n2);
    vars.append("n3", n3);
    vars.append("n4", n4);
    vars.append("p1", p1);
    vars.append("p2", p2);
    vars.append("p3", p3);
    vars.append("p4", p4);
    vars.append("e1", e1);
    vars.append("e2", e2);
    vars.append("e3", e3);
    vars.append("e4", e4);
    vars.append<bool>("mixed", mixed);
    vars.append<double>("delta", [](EuCell& cell) -> double {
        return std::min(std::abs(cell[a]), std::abs(1.0 - cell[a]));
    });

    VtuFile::save("output/mesh.vtu", mesh, vars);

    auto body_central = body(mesh, p1, n1);
    VtuFile::save("output/body(central).vtu", body_central, {}, false, true);

    auto body_youngs = body(mesh, p2, n2);
    VtuFile::save("output/body(youngs).vtu", body_youngs, {}, false, true);

    auto body_elvira = body(mesh, p3, n3);
    VtuFile::save("output/body(elvira).vtu", body_elvira, {}, false, true);

    auto body_csir = body(mesh, p4, n4);
    VtuFile::save("output/body(csir).vtu", body_csir, {}, false, true);
}

// Адаптировать, если ячейка или сосед смешанные
void set_flag(EuCell& cell) {
    cell.set_flag(-1);
    if (mixed(cell)) {
        cell.set_flag(1);
        return;
    }
    for (auto face: cell.faces()) {
        auto neib = face.neib();
        if (mixed(neib)) {
            cell.set_flag(1);
            return;
        }
    }
}

// Выставить значения объемной доли и провести адаптацию
void initialize(EuMesh& mesh, std::function<void(EuCell&)> set_alpha) {
    mesh.set_distributor(Distributor::initializer(set_alpha));
    mesh.for_each(set_alpha);
    if (mesh.adaptive()) {
        for (int i = 0; i <= mesh.max_level(); ++i) {
            mesh.for_each(set_flag);
            mesh.make_shuba(1);
            mesh.refine();
        }
    }
}

// Для плоскости
void show_plain(EuMesh& mesh) {
    initialize(mesh, [](EuCell& cell) {
        double p = -cell.center().dot(some_n);
        cell[a] = quad_volume_fraction(p, some_n, cell.hx(), cell.hy());
    });

    make_interface(mesh);

    auto under_plain = [](const Vector3d& v) -> bool {
        return v.dot(some_n) < 0.0;
    };
    calc_errors(mesh, under_plain, 200);

    save_mesh(mesh);
}

// Для обычной характеристической функции
void show_classic(EuMesh& mesh, InFunction func) {
    initialize(mesh, [func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 100000);
    });

    make_interface(mesh);

    calc_errors(mesh, func, 200);

    save_mesh(mesh);
}

void show_diffuse(EuMesh& mesh) {
    mesh.set_max_level(0);

    mesh.for_each([](EuCell& cell) {
        double alpha = cell.integrate_low(diffuse_func, 10) / cell.volume();
        if (alpha < 1.0e-6) alpha = 0.0;
        if (alpha > 1.0 - 1.0e-6) alpha = 1.0;
        cell[a] = alpha;
    });

    make_interface(mesh);

    save_mesh(mesh);
}

void show_noise(EuMesh& mesh, InFunction func) {
    mesh.set_max_level(0);

    mesh.for_each([func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 1000);
    });

    // Добавить шум
    auto ampl = mesh.add<double>("ampl");
    mesh.for_each([ampl](EuCell& cell) {
        double avg_a = 0.0;
        for (int i = - 2; i <= 2; ++i) {
            for (int j = - 2; j <= 2; ++j) {
                avg_a += cell.neib(i, j)[a];
            }
        }
        avg_a /= 25;
        cell[ampl] = 0.5 - std::abs(avg_a - 0.5);
    });

    mesh.for_each([ampl](EuCell& cell) {
        double r = 2.0 * rand() / double(RAND_MAX) - 1.0;
        double alpha = cell[a] + 0.08 * r * cell[ampl];
        cell[a] = math::between(alpha, 0.0, 1.0);
    });

    make_interface(mesh);

    save_mesh(mesh);
}

int main() {
    utils::mpi::handler mpi_init;
    utils::threads::on();

    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    //gen.set_sizes(30, 50);
    gen.set_size(50);

    EuMesh mesh(gen);
    mesh.set_max_level(1);

    a = mesh.add<double>("a");
    std::tie(p1, p2, p3, p4) = mesh.add<double  >("p1", "p2", "p3", "p4");
    std::tie(n1, n2, n3, n4) = mesh.add<Vector3d>("n1", "n2", "n3", "n4");
    std::tie(e1, e2, e3, e4) = mesh.add<double  >("e1", "e2", "e3", "e4");

    int test = 1;

    switch (test) {
        case 0: show_plain(mesh); return 0;
        case 1: show_classic(mesh, smooth_func); return 0;
        case 2: show_classic(mesh, angle_func); return 0;
        case 3: show_diffuse(mesh); return 0;
        case 4: show_noise(mesh, smooth_func); return 0;
        default: return 0;
    }
}

/*
Classic interface
hx       hy          Central      Youngs      ELVIRA       CSIR
0.04000, 0.04000:    1.04e-02    3.17e-03    3.38e-03    2.64e-03
0.02000, 0.02000:    5.31e-03    9.18e-04    7.65e-04    5.88e-04
0.01000, 0.01000:    2.66e-03    3.45e-04    1.63e-04    1.39e-04
0.00500, 0.00500:    1.33e-03    1.56e-04    4.02e-05    3.52e-05
0.00250, 0.00250:    6.67e-04    7.63e-05    1.01e-05    8.88e-06
0.00125, 0.00125:    3.33e-04    3.71e-05    2.61e-06    2.26e-06
0.00062, 0.00062:    1.67e-04    1.85e-05    6.74e-07    6.03e-07
*/