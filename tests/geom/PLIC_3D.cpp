// PLIC реконструкция интерфейса на трёхмерных декартовых сетках

#include <boost/type_traits/type_with_alignment.hpp>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/sections.h>

#include <zephyr/geom/geom.h>

#include <zephyr/math/funcs.h>

using namespace zephyr;
using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Cuboid;

// Радиус и угол в декартовы координаты
Vector3d to_cartesian(double r, double phi) {
    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

// Нормаль для теста с плоскостью
static Vector3d some_n = (Vector3d{0.24, 0.13, 0.852}).normalized();

// Характеристическая функция (функция-индикатор)
using InFunction = std::function<bool(const Vector3d &)>;

// Пространственная функция
using SpFunction = std::function<double(const Vector3d &)>;


// 0. Нижняя полуплоскость
InFunction plain_func = [](const Vector3d& v) -> bool {
    return v.dot(some_n) < 0.0;
};

// 1. Гладкая замкнутая кривая (тор)
InFunction smooth_func = [](const Vector3d &v) -> bool {
    const double R = 0.6;
    const double r = 0.3;
    double phi = std::atan2(v.y(), v.x());
    Vector3d c = to_cartesian(R, phi);
    return (v - c).norm() < r;
};

// 2. Некоторая угловатая функция
InFunction angle_func = [](const Vector3d &v) -> bool {
    return false;
};


static Storable<double> a;
static Storable<double> a2;
// Простые производные
static Storable<Vector3d> n1;
static Storable<double> p1;
static Storable<double> e1;
// Формула из 2D
static Storable<Vector3d> n2;
static Storable<double> p2;
static Storable<double> e2;
// Моя формула
static Storable<Vector3d> n3;
static Storable<double> p3;
static Storable<double> e3;

inline bool mixed(EuCell& cell) {
    return 0.0 < cell[a] && cell[a] < 1.0;
}

inline double get_ae(double p, const Vector3d& n) {
    if (std::abs(n.z()) > 1.0e-15) {
        double z = (p - 0.5 * n.x() - 0.5 * n.y()) / n.z();
        return math::between(0.5 + z, 0.0, 1.0);
    }
    return n.x() + n.y() > 0.0 ? 0.0 : 1.0;
}

struct out_t {
    double p;
    Vector3d n;
    double ae;
};

out_t edge_fraction(double a0, double a1, double a2) {
    constexpr double eps_a = 1.0e-6;
    if (a1 <= eps_a || a2 <= eps_a) {
        return {.p=NAN, .n={NAN, NAN, NAN}, .ae=0.0};
    }
    if (a1 >= 1.0 - eps_a || a2 >= 1.0 - eps_a) {
        return {.p=NAN, .n={NAN, NAN, NAN}, .ae=1.0};
    }
    if (a0 <= eps_a || a0 >= 1.0 - eps_a) {
        // return 0 или 1
        double ae = 0.5 + 0.5 * math::sign_p(a1 + a2 - 1.0, eps_a);
        return {.p=NAN, .n={NAN, NAN, NAN}, .ae=ae};
    }

    // Остались случаи:
    //  eps_a < a0 < 1 - eps_a,
    //  eps_a < a1 < 1 - eps_a,
    //  eps_a < a2 < 1 - eps_a

    // Это как первая итерация в геометрическом методе, работает отлично!
    // Для простого случая сразу дает ответ.
    Vector3d n = {a0 - a1, a0 - a2, 1}; n.normalize();
    double p = (a0 - 0.5) * n.z();
    double ae = get_ae(p, n);

    constexpr int max_iters = 20;
    constexpr double eps_r = 1.0e-12;
    constexpr double eps_p = 1.0e-12;
    constexpr double eps_n = 1.0e-12;
    constexpr double eps_ae = 1.0e-12;

    for (int counter = 0; counter < max_iters; ++counter) {
        double A0 = cube_volume_fraction(p, n);
        double A1 = cube_volume_fraction(p - n.x(), n);
        double A2 = cube_volume_fraction(p - n.y(), n);

        Vector4d F = { a0 - A0, a1 - A1, a2 - A2, 0.0 };

        // Условие выхода по максимальной невязке
        double err_1 = F.cwiseAbs().maxCoeff();
        if (err_1 < eps_r) {
            break;
        }

        vf_grad_t D0 = cube_volume_fraction_grad(p, n);
        vf_grad_t D1 = cube_volume_fraction_grad(p - n.x(), n);
        vf_grad_t D2 = cube_volume_fraction_grad(p - n.y(), n);

        Matrix4d M;
        M << D0.dp, D0.dn.x(),         D0.dn.y(),         D0.dn.z(),
             D1.dp, D1.dn.x() - D1.dp, D1.dn.y(),         D1.dn.z(),
             D2.dp, D2.dn.x(),         D2.dn.y() - D2.dp, D2.dn.z(),
             0.0,   n.x(),             n.y(),             n.z();

        Vector4d delta = M.inverse() * F;

        double   delta_p = delta[0];
        Vector3d delta_n = {delta[1], delta[2], delta[3]};

        double lambda = 1.0;
        if (n.z() + delta_n.z() < 0.0) {
            // Запрет на выход из диапазона n.z() > 0
            lambda = 0.95 * std::abs(n.z() / delta_n.z());
        }
        delta_p *= lambda;
        delta_n *= lambda;

        n += delta_n; n.normalize();
        p += delta_p;

        double new_ae = get_ae(p, n);

        if (delta_n.norm() < eps_n) {
            break;
        }
        if (std::abs(delta_p) < eps_p) {
            break;
        }
        if (std::abs(ae - new_ae) < eps_ae) {
            ae = new_ae;
            break;
        }
        ae = new_ae;
    }
    return {.p=p, .n=n, .ae=ae};
}

double e2f(double ax1, double ax2, double ay1, double ay2) {
    auto [min_x, max_x] = math::sorted(ax1, ax2);
    auto [min_y, max_y] = math::sorted(ay1, ay2);
    return 0.5 * (min_x * max_y + max_x * min_y + max_x * max_y - min_x * min_y);
}

Vector3d get_normal3D(EuCell& cell, bool verbose = false) {
    double A = cell[a];

    double AL = cell.face(Side3D::L).neib(a);
    double AR = cell.face(Side3D::R).neib(a);
    double AB = cell.face(Side3D::B).neib(a);
    double AT = cell.face(Side3D::T).neib(a);
    double AX = cell.face(Side3D::X).neib(a);
    double AF = cell.face(Side3D::F).neib(a);

    double a_lb = edge_fraction(A, AL, AB).ae;
    double a_lt = edge_fraction(A, AL, AT).ae;;
    double a_lx = edge_fraction(A, AL, AX).ae;;
    double a_lf = edge_fraction(A, AL, AF).ae;;

    double a_rb = edge_fraction(A, AR, AB).ae;;
    double a_rt = edge_fraction(A, AR, AT).ae;;
    double a_rx = edge_fraction(A, AR, AX).ae;;
    double a_rf = edge_fraction(A, AR, AF).ae;;

    double a_bx = edge_fraction(A, AB, AX).ae;;
    double a_bf = edge_fraction(A, AB, AF).ae;;
    double a_bl = edge_fraction(A, AB, AL).ae;;
    double a_br = edge_fraction(A, AB, AR).ae;;

    double a_tx = edge_fraction(A, AT, AX).ae;;
    double a_tf = edge_fraction(A, AT, AF).ae;;
    double a_tl = edge_fraction(A, AT, AL).ae;;
    double a_tr = edge_fraction(A, AT, AR).ae;;

    double a_xl = edge_fraction(A, AX, AL).ae;;
    double a_xr = edge_fraction(A, AX, AR).ae;;
    double a_xb = edge_fraction(A, AX, AB).ae;;
    double a_xt = edge_fraction(A, AX, AT).ae;;

    double a_fl = edge_fraction(A, AF, AL).ae;;
    double a_fr = edge_fraction(A, AF, AR).ae;;
    double a_fb = edge_fraction(A, AF, AB).ae;;
    double a_ft = edge_fraction(A, AF, AT).ae;;

    double a_l = e2f(a_lb, a_lt, a_lx, a_lf);
    double a_r = e2f(a_rb, a_rt, a_rx, a_rf);
    double a_b = e2f(a_bl, a_br, a_bx, a_bf);
    double a_t = e2f(a_tl, a_tr, a_tx, a_tf);
    double a_x = e2f(a_xl, a_xr, a_xb, a_xt);
    double a_f = e2f(a_fl, a_fr, a_fb, a_ft);

    Vector3d n = {a_l - a_r, a_b - a_t, a_x - a_f};
    n.normalize();

    if (verbose) {
        double max_a = std::max({A, AL, AR, AB, AT, AX, AF});
        double min_a = std::min({A, AL, AR, AB, AT, AX, AF});

        std::cout << "A " << A << "\n";
        std::cout << "AL " << AL << "\n";
        std::cout << "AR " << AR << "\n";
        std::cout << "AB " << AB << "\n";
        std::cout << "AT " << AT << "\n";
        std::cout << "AX " << AX << "\n";
        std::cout << "AF " << AF << "\n";

        double h = cell.linear_size();
        double p = cube_find_section(A, some_n) * h;
        double P2 = cube_find_section(A, n) * h;

        EuMesh nice1(3, false);
        Polyhedron poly1 = cell.polyhedron().clip(cell.center() + p * some_n, some_n);
        nice1.push_back(poly1);
        VtuFile::save("output/that_cell1.vtu", nice1, {}, false, true);

        EuMesh nice2(3, false);
        Polyhedron poly2 = cell.polyhedron().clip(cell.center() + P2 * n, n);
        nice2.push_back(poly2);
        VtuFile::save("output/that_cell2.vtu", nice2, {}, false, true);

        std::cout << "a_lb: " << edge_fraction(A, AL, AB).ae << "\n";
        std::cout << "a_lt: " << edge_fraction(A, AL, AT).ae << "\n";
        std::cout << "a_lx: " << edge_fraction(A, AL, AX).ae << "\n";
        std::cout << "a_lf: " << edge_fraction(A, AL, AF).ae << "\n";

        std::cout << "a_rb: " << edge_fraction(A, AR, AB).ae << "\n";
        std::cout << "a_rt: " << edge_fraction(A, AR, AT).ae << "\n";
        std::cout << "a_rx: " << edge_fraction(A, AR, AX).ae << "\n";
        std::cout << "a_rf: " << edge_fraction(A, AR, AF).ae << "\n";

        std::cout << "a_bx: " << edge_fraction(A, AB, AX).ae << "\n";
        std::cout << "a_bf: " << edge_fraction(A, AB, AF).ae << "\n";
        std::cout << "a_bl: " << edge_fraction(A, AB, AL).ae << "\n";
        std::cout << "a_br: " << edge_fraction(A, AB, AR).ae << "\n";

        std::cout << "a_tx: " << edge_fraction(A, AT, AX).ae << "\n";
        std::cout << "a_tf: " << edge_fraction(A, AT, AF).ae << "\n";
        std::cout << "a_tl: " << edge_fraction(A, AT, AL).ae << "\n";
        std::cout << "a_tr: " << edge_fraction(A, AT, AR).ae << "\n";

        std::cout << "a_xl: " << edge_fraction(A, AX, AL).ae << "\n";
        std::cout << "a_xr: " << edge_fraction(A, AX, AR).ae << "\n";
        std::cout << "a_xb: " << edge_fraction(A, AX, AB).ae << "\n";
        std::cout << "a_xt: " << edge_fraction(A, AX, AT).ae << "\n";

        std::cout << "a_fl: " << edge_fraction(A, AF, AL).ae << "\n";
        std::cout << "a_fr: " << edge_fraction(A, AF, AR).ae << "\n";
        std::cout << "a_fb: " << edge_fraction(A, AF, AB).ae << "\n";
        std::cout << "a_ft: " << edge_fraction(A, AF, AT).ae << "\n\n";

        std::cout << "a_l: " << e2f(a_lb, a_lt, a_lx, a_lf) << "\n";
        std::cout << "a_r: " << e2f(a_rb, a_rt, a_rx, a_rf) << "\n";
        std::cout << "a_b: " << e2f(a_bl, a_br, a_bx, a_bf) << "\n";
        std::cout << "a_t: " << e2f(a_tl, a_tr, a_tx, a_tf) << "\n";
        std::cout << "a_x: " << e2f(a_xl, a_xr, a_xb, a_xt) << "\n";
        std::cout << "a_f: " << e2f(a_fl, a_fr, a_fb, a_ft) << "\n";

        std::cout << "some_n: " << some_n.transpose() << "\n";
        std::cout << "mine n: " << n.transpose() << "\n";
        //std::cout << "POLY1\n" << poly1 << "\n";
        //std::cout << "POLY2\n" << poly2 << "\n";
    }

    return n;
}

void make_interface(EuMesh& mesh) {
    for (EuCell& cell: mesh) {
        if (!mixed(cell)) {
            cell[n1] = Vector3d::Zero();
            cell[n2] = Vector3d::Zero();
            cell[n3] = Vector3d::Zero();
            cell[p1] = 0.0;
            cell[p2] = 0.0;
            cell[p3] = 0.0;
            continue;
        }

        Vector3d grad1 = Vector3d::Zero();
        Vector3d grad2 = Vector3d::Zero();
        for (auto& face: cell.faces()) {
            double a_f1 = 0.5 * (cell[a] + face.neib(a));
            Vector3d S = face.area() * face.normal();
            grad1 -= a_f1 * S;

            double a_f2 = face_fraction(cell[a], face.neib(a));
            grad2 -= a_f2 * S;
        }
        cell[n1] = grad1.normalized();
        cell[n2] = grad2.normalized();
        cell[n3] = get_normal3D(cell);

        Vector3d c = cell.center();

        //if (std::abs(c.x()) < 0.85 && std::abs(c.y()) < 0.85 && std::abs(c.z()) < 0.35 && cell[n3].dot(some_n) < 0.9999) {
        //    Vector3d kek = get_normal3D(cell, true);
        //    throw std::runtime_error("enough");
        //}

        double h = cell.linear_size();
        cell[p1] = cube_find_section(cell[a], cell[n1]) * h;
        cell[p2] = cube_find_section(cell[a], cell[n2]) * h;
        cell[p3] = cube_find_section(cell[a], cell[n3]) * h;
    }
}

EuMesh body(EuMesh& mesh, Storable<double> p, Storable<Vector3d> n) {
    EuMesh clipped(3, false);

    for (auto& cell: mesh) {
        if (cell[a] <= 0.0 || (cell[a] < 0.5 && cell[n].isZero())) {
            continue;
        }
        if (cell[a] >= 1.0 || (cell[a] > 0.5 && cell[n].isZero())) {
            clipped.push_back(cell.polyhedron());
            continue;
        }

        EuMesh clip_i(3, false);

        Vector3d point = cell.center() + cell[p] * cell[n];
        auto poly1 = cell.polyhedron();
        auto poly = poly1.clip(point, cell[n]);
        if (poly.checkout() != 0) {
            std::cout << "WTF? " << poly.checkout() << " " << cell[a] << " " << cell[p] << " " << cell[n].transpose() << "\n";
            continue;
        }
        clipped.push_back(poly);
    }
    return clipped;
}

std::tuple<double, double> calc_errors(EuMesh& mesh, InFunction func, int nx) {
    double err1 = 0.0;
    double err2 = 0.0;

    mesh.for_each([func, nx, &err1, &err2](EuCell& cell) {
        if (!mixed(cell)) {
            cell[e1] = 0.0;
            cell[e2] = 0.0;
            return;
        }

        Box box = cell.bbox();

        double pos1 = cell[p1];
        double pos2 = cell[p2];
        Vector3d norm1 = cell[n1];
        Vector3d norm2 = cell[n2];

        int counter1 = 0;
        int counter2 = 0;
        for (int i = 0; i < nx; ++i) {
            double x = box.vmin.x() + (box.vmax.x() - box.vmin.x()) * (i + 0.5) / nx;
            for (int j = 0; j < nx; ++j) {
                double y = box.vmin.y() + (box.vmax.y() - box.vmin.y()) * (j + 0.5) / nx;
                Vector3d r = {x, y, 0.0};
                bool in0 = func(r);
                bool in1 = (r - cell.center()).dot(norm1) < pos1;
                bool in2 = (r - cell.center()).dot(norm2) < pos2;
                if (in0 != in1) { ++counter1; }
                if (in0 != in2) { ++counter2; }
            }
        }
        if (counter1 == nx * nx) {
            // Выяснить, что здесь такое вообще
            cell[e1] = 0.0;
            cell[e2] = 0.0;
        }
        else {
            cell[e1] = double(counter1) / (nx * nx);// * cell.volume();
            cell[e2] = double(counter2) / (nx * nx);// * cell.volume();

            err1 += cell[e1];
            err2 += cell[e2];
        }
    });

    int mix_count = mesh.sum([](EuCell& cell) -> int {
        return mixed(cell) ? 1 : 0;
    }, 0);
    err1 /= mix_count;
    err2 /= mix_count;

    return {err1, err2};
}

void save_mesh(EuMesh& mesh) {
    Variables vars;
    vars.append("a", a);
    vars.append("a2", a2);
    vars.append("p1", p1);
    vars.append("p2", p2);
    vars.append("p3", p3);
    vars.append("n1", n1);
    vars.append("n2", n2);
    vars.append("n3", n3);
    vars.append("e1", e1);
    vars.append("e2", e2);
    vars.append("e3", e3);

    VtuFile::save("output/mesh.vtu", mesh, vars);

    auto body1 = body(mesh, p1, n1);
    VtuFile::save("output/body1.vtu", body1, {}, false, true);

    auto body2 = body(mesh, p2, n2);
    VtuFile::save("output/body2.vtu", body2, {}, false, true);

    auto body3 = body(mesh, p3, n3);
    VtuFile::save("output/body3.vtu", body3, {}, false, true);
}

// Для плоскости
void show_plain(EuMesh& mesh) {
    mesh.for_each([](EuCell& cell) {
        double h = cell.linear_size();
        double p = -cell.center().dot(some_n) / h;
        cell[a] = cube_volume_fraction(p, some_n);
    });

    make_interface(mesh);

    //auto [error1, error2] = calc_errors(mesh, func, 10);
    //std::cout << "Error1: " << error1 << "\n";
    //std::cout << "Error2: " << error2 << "\n";

    save_mesh(mesh);
}

// Для обычной характеристической функции
void show_classic(EuMesh& mesh, InFunction func) {
    mesh.for_each([func](EuCell& cell) {
        cell[a] = cell.volume_fraction(func, 100000);
    });

    make_interface(mesh);

    //auto [error1, error2] = calc_errors(mesh, func, 10);
    //std::cout << "Error1: " << error1 << "\n";
    //std::cout << "Error2: " << error2 << "\n";

    save_mesh(mesh);
}

// В норме L_inf: ц.р. не сходятся, мой метод с 1-ым порядком
// В норме L_1: ц.р. имеют 1-ый порядок, мой метод 2-ой
// Средняя от L_inf: ц.р. имеют 1-ый порядок, мой метод 2-ой
// Во всех случаях обычные разности дают погрешности на порядок больше
int main() {
    utils::threads::on();

    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -0.5, 0.5);
    gen.set_nx(40);

    EuMesh mesh(gen);

    a = mesh.add<double>("a");
    a2 = mesh.add<double>("a2");
    p1 = mesh.add<double>("p1");
    n1 = mesh.add<Vector3d>("n1");
    e1 = mesh.add<double>("e1");
    p2 = mesh.add<double>("p2");
    n2 = mesh.add<Vector3d>("n2");
    e2 = mesh.add<double>("e2");
    p3 = mesh.add<double>("p3");
    n3 = mesh.add<Vector3d>("n3");
    e3 = mesh.add<double>("e3");

    int test = 1;

    switch (test) {
        case 0: show_plain(mesh); return 0;
        case 1: show_classic(mesh, smooth_func); return 0;
        //case 2: show_classic(mesh, angle_func); return 0;
        //case 3: show_diffuse(mesh); return 0;
        //case 4: show_noise(mesh, smooth_func); return 0;
        default: return 0;
    }
}