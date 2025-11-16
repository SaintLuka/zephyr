#include <iostream>
#include <map>

#include <zephyr/math/funcs.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/io/vtr_file.h>
#include <zephyr/utils/pyplot.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr;
using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::math;
using namespace zephyr::io;
using namespace zephyr::utils;

using generator::Cuboid;

// Записать три исходных кубика
void save_origin() {
    auto poly1 = Polyhedron::Cube();
    auto poly2 = Polyhedron::Cube();
    auto poly3 = Polyhedron::Cube();

    poly2.move(Vector3d::UnitX());
    poly3.move(Vector3d::UnitY());

    EuMesh orig(3, false);
    orig.push_back(poly1);
    orig.push_back(poly2);
    orig.push_back(poly3);

    VtuFile::save("out/orig.vtu", orig, {});
}

// Записать последовательность сечений
void save_history(const std::vector<double>& ps, const std::vector<Vector3d>& ns) {
    auto poly1 = Polyhedron::Cube();
    auto poly2 = Polyhedron::Cube();
    auto poly3 = Polyhedron::Cube();

    poly2.move(Vector3d::UnitX());
    poly3.move(Vector3d::UnitY());

    PvdFile pvd("clipped", "out");
    pvd.polyhedral = true;

    for (int i = 0; i < ps.size(); ++i) {
        Vector3d p = ps[i] * ns[i];

        auto clip1 = poly1.clip(p, ns[i]);
        auto clip2 = poly2.clip(p, ns[i]);
        auto clip3 = poly3.clip(p, ns[i]);

        EuMesh clipped(3, false);
        clipped.push_back(clip1);
        clipped.push_back(clip2);
        clipped.push_back(clip3);
        pvd.save(clipped, i);
    }
}

// Типы сечений трёх кубов
struct cases_t {
    int c1{-13}, c2{-13}, c3{-13};

    cases_t() = default;

    // Типы сечений для заданной прямой
    cases_t(double p, const Vector3d& n) {
        c1 = cube_section_case(p, n);
        c2 = cube_section_case(p - n.x(), n);
        c3 = cube_section_case(p - n.y(), n);

        // Простейший случай (есть чистая ячейка)
        if (std::abs(c1) == 0 || std::abs(c2) == 0 || std::abs(c3) == 0) {
            c1 = c2 = c3 = -13;
        }

        // Инвертируем относительно основной
        if (c1 < 0) { c1 *= -1; c2 *= -1; c3 *= -1; }

        // c2, c3 эквивалентны, упорядочим, чтобы снизить размерность
        if (c2 > c3) std::swap(c2, c3);
    }

    // Проверка на равенство
    bool operator!=(const cases_t &other) const {
        return c1 != other.c1 || c2 != other.c2 || c3 != other.c3;
    }

    // Сравнение в лексикографическом порядке
    bool operator<(const cases_t &other) const {
        if (c1 < other.c1) { return true; }
        if (c1 > other.c1) { return false; }

        if (c2 < other.c2) { return true; }
        if (c2 > other.c2) { return false; }

        return c3 < other.c3;
    }
};

// Получить сечение ребра для заданной плоскости
inline double get_ae(double p, const Vector3d& n) {
    if (std::abs(n.z()) > 1.0e-15) {
        double z = (p - 0.5 * n.x() - 0.5 * n.y()) / n.z();
        return between(0.5 + z, 0.0, 1.0);
    }
    return n.x() + n.y() > 0.0 ? 0.0 : 1.0;
}

// Ошибка в определении объемных долей
inline double get_error(double p, const Vector3d& n, double a0, double a1, double a2) {
    Vector3d errs = {cube_volume_fraction(p, n) - a0,
                     cube_volume_fraction(p - n.x(), n) - a1,
                     cube_volume_fraction(p - n.y(), n) - a2};
    return errs.cwiseAbs().maxCoeff();
}

// Последовательность плоскостей во время итераций
struct res_t {
    double   p = NAN;              // Расстояние до плоскости
    Vector3d n = {NAN, NAN, NAN};  // Нормаль плоскости
    double   ae = NAN;             // Искомое сечение ребра

    int count    = -1; // Число итераций
    int out_case = -1; // Тип завершения

    // При записи данных с итераций
    std::vector<double>   ps_history;
    std::vector<Vector3d> ns_history;

    res_t() = default;

    // Когда хотим просто установить ae и не париться
    explicit res_t(double ae_in) {
        p = NAN;
        n = {NAN, NAN, NAN};
        ae = ae_in;
        count = 0;
        ps_history.push_back(p);
        ns_history.push_back(n);
    }

    // Когда хотим просто установить плоскость и не париться
    res_t(double p_in, const Vector3d& n_in) {
        append(p_in, n_in);
    }

    // Добавить плоскость
    void append(double p_in, const Vector3d& n_in) {
        ps_history.push_back(p_in);
        ns_history.push_back(n_in);
        p = p_in;
        n = n_in;
        ae = get_ae(p, n);
    }

    void plot(double a0, double a1, double a2) {
        auto steps = np::arange<double>(ps_history.size());
        std::vector<double> errs(ps_history.size());
        std::vector<double> a_es(ps_history.size());

        for (size_t i = 0; i < ps_history.size(); ++i) {
            errs[i] = get_error(ps_history[i], ns_history[i], a0, a1, a2);
            a_es[i] = get_ae(ps_history[i], ns_history[i]);
        }

        std::cout << "a0, a1, a2: " << a0 << ", " << a1 << ", " << a2 << "\n";
        std::cout << "  p, n: " << ps_history.back() << ", " << ns_history.back().transpose() << "\n";
        std::cout << "  ae: " << a_es.back() << "\n";
        std::cout << "  count: " << count << "\n";
        std::cout << "  out_case: " << out_case << "\n";
        std::cout << "  error: " << errs.back() << "\n";

        // Отклонение от последнего значения
        for (size_t i = 0; i < ps_history.size(); ++i) {
            //a_es[i] = std::abs(a_es[i] - a_es.back());
        }

        pyplot plt;

        plt.figure({.figsize={9.0, 4.0}, .dpi=170});

        plt.subplot(1, 2, 1);
        plt.title("Погрешность от итерации");
        plt.grid(true);
        plt.semilogy(steps, errs, "k");

        plt.subplot(1, 2, 2);
        plt.title("Погрешность $\\alpha_e$");
        plt.grid(true);
        plt.plot(steps, a_es, "k");

        plt.tight_layout();
        plt.show();
    }
};

// Опции для итерационных процедур
struct IterOpts {
    int max_iters = 50;

    bool series  = false;
    bool verbose = false;
};

// Почти как предыдущее, но решается система трёх + одно уравнений.
// В данном случае ищется нормаль и параметр p.
// Считает очень быстро, и сходится как будто бы везде при dumper = 0.95.
// Понятно, что среднее число итераций < 10, но есть тяжелая сходимость
// на границах куба.
// Наверное, лучший из всех вариантов.
// Всего 7 итераций не считая тех областей, где сходимости нет.
// Причем можно сразу ставить дампер равный 1.0.
// Эхх, вот если бы я только знал области предельные, чтобы просто сразу их исключать.
res_t edge_fraction(double a0, double a1, double a2, const IterOpts& opts) {
    constexpr double eps_a = 1.0e-6;
    if (a1 <= eps_a || a2 <= eps_a) {
        return res_t(0.0);
    }
    if (a1 >= 1.0 - eps_a || a2 >= 1.0 - eps_a) {
        return res_t(1.0);
    }
    if (a0 <= eps_a || a0 >= 1.0 - eps_a) {
        // return 0 или 1
        double ae = 0.5 + 0.5 * sign_p(a1 + a2 - 1.0, eps_a);
        return res_t(ae);
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

    constexpr double eps_r = 1.0e-12;
    constexpr double eps_n = 1.0e-12;
    constexpr double eps_ae = 1.0e-12;

    res_t series;
    int counter = 0;
    int out_case = 10;
    for (; counter < opts.max_iters; ++counter) {
        if (opts.series) { series.append(p, n); }

        double A0 = cube_volume_fraction(p, n);
        double A1 = cube_volume_fraction(p - n.x(), n);
        double A2 = cube_volume_fraction(p - n.y(), n);

        Vector4d F = { a0 - A0, a1 - A1, a2 - A2, 0.0 };

        // Условие выхода по максимальной невязке
        double err_1 = F.cwiseAbs().maxCoeff();
        if (err_1 < eps_r) {
            if (opts.verbose) { std::cout << "break 1 (residual)\n"; }
            out_case = 1;
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
            if (opts.verbose) { std::cout << "break 2 (delta normal)\n"; }
            out_case = 2;
            break;
        }

        if (std::abs(ae - new_ae) < eps_ae) {
            if (opts.verbose) { std::cout << "break 3 (delta ae)\n"; }
            out_case = 3;
            break;
        }
        ae = new_ae;
    }
    series.append(p, n);
    series.count = counter;
    series.out_case = out_case;
    return series;
}

inline std::vector<double> chebyshev(int n) {
    std::vector res(n, 0.0);
    for (int i = 0; i < n; ++i) {
        res[i] = 0.5 + 0.5 * (-std::cos((2 * i + 1) * M_PI / (2 * n)) / std::cos(M_PI / (2 * n)));
    }
    return res;
}

// Построить функцию edge_fraction(a0, a1, a2);
void ae_full(const std::string& filename) {
    threads::on();

    int N = 101;
    std::vector xs = chebyshev(N);
    std::vector ys = chebyshev(N);
    std::vector zs = chebyshev(N);

#define table(x) std::vector(N, std::vector(N, std::vector(N, x)))

    auto a_e = table(0.0);
    auto iters = table(0);
    auto error = table(0.0);
    auto cases = table(cases_t{});
    auto type = table(0);
    auto out_case = table(0);
    auto dist = table(0.0);
    auto norm = table(Vector3d{});

    Stopwatch elapsed(true);
    threads::parallel_for(0, N, [&](int i) {
        double a1 = xs[i];
        for (int j = 0; j < N; ++j) {
            double a2 = ys[j];
            for (int k = 0; k < N; ++k) {
                double a0 = zs[k];

                auto res = edge_fraction(a0, a1, a2, {.max_iters=50});

                a_e [i][j][k] = res.ae;

                dist[i][j][k] = res.p;
                norm[i][j][k] = res.n;
                iters[i][j][k] = res.count;
                error[i][j][k] = get_error(res.p, res.n, a0, a1, a2);
                cases[i][j][k] = cases_t(res.p, res.n);
                out_case[i][j][k] = res.out_case;
            }
        }
    });

    // Внешние границы, a_i in {0, 1}
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            a_e[0  ][i][j] = 0.0;
            a_e[N-1][i][j] = 1.0;
            a_e[i][0  ][j] = 0.0;
            a_e[i][N-1][j] = 1.0;
            a_e[i][j][0]   = heav(xs[i] + ys[j] - 1.0);
            a_e[i][j][N-1] = heav(xs[i] + ys[j] - 1.0);
        }
    }
    // Плоскость a_e(a0, a1, 1 - a1) = 0.5
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            a_e[i][N - i - 1][j] = 0.5;
        }
    }

    // Считаем число случаев
    std::map<cases_t, int> all_cases;
    for (int i = 0; i < xs.size(); ++i) {
        for (int j = 0; j < ys.size(); ++j) {
            for (int k = 0; k < zs.size(); ++k) {
                all_cases[cases[i][j][k]] = 0;
            }
        }
    }

    int counter = 0;
    for (auto& c: all_cases) {
        c.second = counter;
        ++counter;
    }

    // Пронумеровать случаи единственным индексом
    for (int i = 0; i < xs.size(); ++i) {
        for (int j = 0; j < ys.size(); ++j) {
            for (int k = 0; k < zs.size(); ++k) {
                type[i][j][k] = all_cases[cases[i][j][k]];
            }
        }
    }
    elapsed.stop();

    std::cout << "Cases number: " << all_cases.size() << "\n";
    std::cout << "Compute ae (elapsed): " << elapsed.extended_time() << "\n";

    VtrFile file(filename, xs, ys, zs);

    file.point_data("a_edge", a_e);
    file.point_data("n_iters", iters);
    file.point_data("type", type);

    file.save();
}

int main() {
    save_origin();

    double a0 = 0.95;
    double a1 = 0.125;
    double a2 = 0.479167;

    res_t res = edge_fraction(a0, a1, a2, {.max_iters=150, .series=true, .verbose=true});
    save_history(res.ps_history, res.ns_history);
    res.plot(a0, a1, a2);

    ae_full("out/edge_fraction.vtr");

    return 0;
}
