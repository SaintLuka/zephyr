// PLIC реконструкция интерфейса на трёхмерных декартовых сетках
#include "PLIC.h"

// Нормаль для теста с плоскостью
static Vector3d some_n = Vector3d{0.24512, 0.48123, 0.91253}.normalized();

// 0. Область под плоскостью
auto plain_func = [](const Vector3d& v) -> bool {
    return v.dot(some_n) < 0.0;
};

// 1. Сфера
InFunction sphere_func = [](const Vector3d &v) -> bool {
    const double R = 0.8;
    return (v + 0.5 * Vector3d::UnitZ()).norm() < R;
};

// 2. Тор
InFunction torus_func = [](const Vector3d &v) -> bool {
    const double R = 0.55;
    const double r = 0.28;
    double phi = std::atan2(v.y(), v.x());
    Vector3d c = to_cartesian(R, phi);
    return (v - c).norm() < r;
};

// 3. Сфера с диффузной границей
SpFunction diffuse_func = [](const Vector3d &v) -> double {
    const double R = 0.75;
    const double h = 0.2;
    double r = (v + 0.5 * Vector3d::UnitZ()).norm();
    return 0.5 * (1 + math::sign_p(R - r, h));
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

// Моя формула (2D)
static Storable<Vector3d> n4;
static Storable<double> p4;
static Storable<double> e4;

// Моя формула (3D)
static Storable<Vector3d> n5;
static Storable<double> p5;
static Storable<double> e5;

// Смешанная ячейка?
inline bool mixed(EuCell& cell) {
    return 0.0 < cell[a] && cell[a] < 1.0;
}

inline bool boundary(EuCell& cell) {
    for (auto face: cell.faces()) {
        if (face.is_boundary()) {
            return true;
        }
    }
    return false;
}

void make_interface(EuMesh& mesh) {
    mesh.for_each([](EuCell& cell) {
        if (!mixed(cell)) {
            for (auto n: {n1, n2, n3, n4, n5}) cell[n] = Vector3d::Zero();
            for (auto p: {p1, p2, p3, p4, p5}) cell[p] = 0.0;
            for (auto e: {e1, e2, e3, e4, e5}) cell[e] = 0.0;
            return;
        }

        // Размеры ячейки
        double hx = cell.hx();
        double hy = cell.hy();
        double hz = cell.hz();

        // Простая производная
        Vector3d grad = Vector3d::Zero();
        for (auto face: cell.faces()) {
            double a_f = 0.5 * (cell[a] + face.neib(a));
            grad += a_f * face.area_n();
        }
        cell[n1] = -grad.normalized();

        Stencil3D C(cell, a);
        cell[n2] = C.Youngs(hx, hy, hz);
        cell[n3] = C.ELVIRA(hx, hy, hz);

        // Моя формула (2D)
        grad = Vector3d::Zero();
        for (auto face: cell.faces()) {
            double a_f = face_fraction(cell[a], face.neib(a));
            grad += a_f * face.area_n();
        }
        cell[n4] = -grad.normalized();

        // Моя формула (3D)
        // Собрать доли в соседних ячейках
        std::array<double, Side3D::count()> a_neib;
        for (auto side: Side3D::items()) {
            a_neib[side] = cell.face(side).neib(a);
        }

        // Объемные доли на гранях
        auto a_f = face_fractions(cell[a], a_neib);

        // И обычный Гаусс
        grad = Vector3d::Zero();
        for (auto face: cell.faces()) {
            grad += a_f[face.side()] * face.area_n();
        }
        cell[n5] = -grad.normalized();

        cell[p1] = cube_find_section(cell[a], cell[n1], hx, hy, hz);
        cell[p2] = cube_find_section(cell[a], cell[n2], hx, hy, hz);
        cell[p3] = cube_find_section(cell[a], cell[n3], hx, hy, hz);
        cell[p4] = cube_find_section(cell[a], cell[n4], hx, hy, hz);
        cell[p5] = cube_find_section(cell[a], cell[n5], hx, hy, hz);
    });
}

void calc_errors(EuMesh& mesh, InFunction func, int nx) {
    // Число несовпадающих точек
    std::atomic<size_t> err1_l1{0}, err1_li{0};
    std::atomic<size_t> err2_l1{0}, err2_li{0};
    std::atomic<size_t> err3_l1{0}, err3_li{0};
    std::atomic<size_t> err4_l1{0}, err4_li{0};
    std::atomic<size_t> err5_l1{0}, err5_li{0};

    mesh.for_each([&, nx](EuCell& cell) {
        // Нулевые погрешности
        cell[e1] = cell[e2] = cell[e3] = cell[e4] = cell[e5] = 0.0;

        if (!mixed(cell)) return;
        if (boundary(cell)) return;

        double hx = cell.hx();
        double hy = cell.hy();
        double hz = cell.hz();

        std::array counter = {0, 0, 0, 0, 0};
        for (int i = 0; i < nx; ++i) {
            double x = hx * ((i + 0.5) / nx - 0.5);
            for (int j = 0; j < nx; ++j) {
                double y = hy * ((j + 0.5) / nx - 0.5);
                for (int k = 0; k < nx; ++k) {
                    double z = hz * ((k + 0.5) / nx - 0.5);

                    // Точка относительно центра ячейки
                    Vector3d r = {x, y, z};

                    bool inside = func(cell.center() + r);
                    if (inside != (r.dot(cell[n1]) < cell[p1])) { ++counter[0]; }
                    if (inside != (r.dot(cell[n2]) < cell[p2])) { ++counter[1]; }
                    if (inside != (r.dot(cell[n3]) < cell[p3])) { ++counter[2]; }
                    if (inside != (r.dot(cell[n4]) < cell[p4])) { ++counter[3]; }
                    if (inside != (r.dot(cell[n5]) < cell[p5])) { ++counter[4]; }
                }
            }
        }

        // Интегральная метрика L1, как у Aulisa
        double xi = cell.volume() / (nx * nx * nx);
        cell[e1] = xi * counter[0];
        cell[e2] = xi * counter[1];
        cell[e3] = xi * counter[2];
        cell[e4] = xi * counter[3];
        cell[e5] = xi * counter[4];

        err1_l1 += counter[0];
        err2_l1 += counter[1];
        err3_l1 += counter[2];
        err4_l1 += counter[3];
        err5_l1 += counter[4];

        update_max(err1_li, counter[0]);
        update_max(err2_li, counter[1]);
        update_max(err3_li, counter[2]);
        update_max(err4_li, counter[3]);
        update_max(err5_li, counter[4]);
    });

    double hx{NAN}, hy{NAN}, hz{NAN};
    for (auto cell: mesh) {
        if (mixed(cell)) {
            hx = cell.hx();
            hy = cell.hy();
            hz = cell.hz();
            break;
        }
    }

    std::cout << "hx, hy, hz                    Central      Youngs     CSIR  2D     ELVIRA     CSIR  3D\n";

    std::cout << std::setprecision(5) << std::fixed;
    std::cout << hx << ", " << hy << ", " << hz << ":";
    std::cout << std::setprecision(2) << std::scientific;

    double li_coeff = 1.0 / (nx * nx * nx);
    std::cout << std::setw(12) << li_coeff * err1_li;
    std::cout << std::setw(12) << li_coeff * err2_li;
    std::cout << std::setw(12) << li_coeff * err4_li;
    std::cout << std::setw(12) << li_coeff * err3_li;
    std::cout << std::setw(12) << li_coeff * err5_li << "\t (L_inf norm)\n";

    std::cout << std::setprecision(5) << std::fixed;
    std::cout << hx << ", " << hy << ", " << hz << ":";
    std::cout << std::setprecision(2) << std::scientific;

    double l1_coeff = hx * hy * hz / (nx * nx * nx);
    std::cout << std::setw(12) << l1_coeff * err1_l1;
    std::cout << std::setw(12) << l1_coeff * err2_l1;
    std::cout << std::setw(12) << l1_coeff * err4_l1;
    std::cout << std::setw(12) << l1_coeff * err3_l1;
    std::cout << std::setw(12) << l1_coeff * err5_l1 << "\t (L_1 norm)\n";

    std::cout.flush();
}

EuMesh body(EuMesh& mesh, Storable<double> p, Storable<Vector3d> n) {
    EuMesh clipped(3, false);

    for (auto& cell: mesh) {
        if (boundary(cell)) continue;

        if (cell[a] <= 0.0 || (cell[a] < 0.5 && cell[n].isZero())) {
            continue;
        }
        if (cell[a] >= 1.0 || (cell[a] > 0.5 && cell[n].isZero())) {
            clipped.push_back(cell.polyhedron());
            continue;
        }

        Vector3d point = cell.center() + cell[p] * cell[n];
        auto poly = cell.polyhedron();
        auto clip = poly.clip(point, cell[n]);
        if (clip.checkout() != 0) {
            std::cout << "WTF? " << clip.checkout() << " " << cell[a] << " " << cell[p] << " "
                      << cell[n].transpose() << "\n";
            continue;
        }
        if (!clip.empty()) {
            clipped.push_back(clip);
        }
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
    vars.append("n5", n5);
    vars.append("p1", p1);
    vars.append("p2", p2);
    vars.append("p3", p3);
    vars.append("p4", p4);
    vars.append("p5", p5);
    vars.append("e1", e1);
    vars.append("e2", e2);
    vars.append("e3", e3);
    vars.append("e4", e4);
    vars.append("e5", e5);
    vars.append<bool>("mixed", mixed);
    vars.append<double>("delta", [](EuCell& cell) -> double {
        return std::min(std::abs(cell[a]), std::abs(1.0 - cell[a]));
    });

    VtuFile::save("output/mesh.vtu", mesh, vars);

    auto body_central = body(mesh, p1, n1);
    VtuFile::save("output/body(central).vtu", body_central, {}, true);

    auto body_youngs = body(mesh, p2, n2);
    VtuFile::save("output/body(youngs).vtu", body_youngs, {}, true);

    auto body_elvira = body(mesh, p3, n3);
    VtuFile::save("output/body(elvira).vtu", body_elvira, {}, true);

    auto body_csir_2D = body(mesh, p4, n4);
    VtuFile::save("output/body(csir_2D).vtu", body_csir_2D, {}, true);

    auto body_csir_3D = body(mesh, p5, n5);
    VtuFile::save("output/body(csir_3D).vtu", body_csir_3D, {}, true);
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
    std::cout << "  setup init\t(" << mesh.n_cells() << " cells)\n";
    mesh.for_each(set_alpha);
    if (mesh.adaptive()) {
        for (int i = 0; i <= mesh.max_level(); ++i) {
            std::cout << "  setup " << i << " / " << mesh.max_level() << "\t(" << mesh.n_cells() << " cells)\n";
            mesh.for_each(set_flag);
            mesh.make_shuba(2);
            mesh.refine();
        }
    }
}

// Для плоскости
void show_plain(EuMesh& mesh) {
    initialize(mesh, [](EuCell& cell) {
        double p = -cell.center().dot(some_n);
        cell[a] = cube_volume_fraction(p, some_n, cell.hx(), cell.hy(), cell.hz());
    });

    make_interface(mesh);

    calc_errors(mesh, plain_func, 200);

    save_mesh(mesh);
}

// Для обычной характеристической функции
void show_classic(EuMesh& mesh, InFunction func) {
    std::cout << "Initialize\n";
    initialize(mesh, [func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 1000000);
    });

    std::cout << "Make interfaces\n";
    make_interface(mesh);

    std::cout << "Compute errors\n";
    calc_errors(mesh, func, 200);

    std::cout << "Save mesh\n";
    save_mesh(mesh);
}

void show_diffuse(EuMesh& mesh) {
    mesh.set_max_level(0);

    mesh.for_each([](EuCell& cell) {
        double alpha = cell.integrate_low(diffuse_func, 10) / cell.volume();
        if (alpha < 1.0e-8) alpha = 0.0;
        if (alpha > 1.0 - 1.0e-8) alpha = 1.0;
        cell[a] = alpha;
    });

    make_interface(mesh);

    save_mesh(mesh);
}

// Исследование сходимости
void convergence(EuMesh& mesh, InFunction func, int nx) {
    auto set_alpha = [func, nx](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, nx * nx * nx);
    };

    mesh.set_distributor(Distributor::initializer(set_alpha));
    mesh.for_each(set_alpha);
    make_interface(mesh);
    calc_errors(mesh, func, nx);

    if (mesh.adaptive()) {
        for (int i = 0; i < mesh.max_level(); ++i) {
            mesh.for_each(set_flag);
            mesh.make_shuba(2);
            mesh.refine();
            make_interface(mesh);
            calc_errors(mesh, func, nx);
        }
    }
    save_mesh(mesh);
}

int main() {
    utils::mpi::handler mpi_init;
    utils::threads::on();

    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -0.5, 0.5);
    gen.set_nx(22);
    //gen.set_sizes(21, 17, 11);

    EuMesh mesh(gen);
    mesh.set_max_level(0);

    a = mesh.add<double>("a");
    std::tie(p1, p2, p3, p4, p5) = mesh.add<double  >("p1", "p2", "p3", "p4", "p5");
    std::tie(n1, n2, n3, n4, n5) = mesh.add<Vector3d>("n1", "n2", "n3", "n4", "n5");
    std::tie(e1, e2, e3, e4, e5) = mesh.add<double  >("e1", "e2", "e3", "e4", "e5");

    int test = 6;

    switch (test) {
        case 0: show_plain(mesh); break;
        case 1: show_classic(mesh, sphere_func); break;
        case 2: show_classic(mesh, torus_func); break;
        case 3: show_diffuse(mesh); break;
        case 4: convergence(mesh, sphere_func, 200); break;
        case 5: convergence(mesh, torus_func, 200); break;
        default: break;
    }
    std::cout << "\a";
    return 0;
}