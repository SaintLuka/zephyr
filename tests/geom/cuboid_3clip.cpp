#include <iostream>
#include <map>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/utils/numpy.h>
#include <zephyr/utils/pyplot.h>
#include <zephyr/math/calc/roots.h>
#include <zephyr/math/funcs.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr;
using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::math;
using namespace zephyr::io;

using generator::Cuboid;
using utils::Stopwatch;

// Записать три исходных кубика
void save_origin();

Vector3d plot_target(double a0,
    const std::vector<double>& as,
    const std::vector<Vector3d>& cs,
    const std::vector<double>& ws);

// Типы сечений трёх кубов
struct cases_t {
    int c1{-13}, c2{-13}, c3{-13};

    cases_t() = default;

    cases_t(int c1_in, int c2_in, int c3_in);

    // Проверка на равенство
    bool operator!=(const cases_t &other) const;

    // Сравнение в лексикографическом порядке
    bool operator<(const cases_t &other) const;
};

struct Info {
    double   p;       // Точка плоскости
    Vector3d n;       // Нормаль плоскости
    double   ae;      // Искомое сечение ребра

    cases_t  cases;   // Типы сечений кубов

    int      count;   // Число итераций
    Vector3d errs;    // Погрешности каждой
    double   max_err; // Максимальная из трёх

    static double get_ae(double p, const Vector3d& n) {
        if (std::abs(n.z()) > 1.0e-15) {
            double z = (p - 0.5 * n.x() - 0.5 * n.y()) / n.z();
            return between(0.5 + z, 0.0, 1.0);
        }
        return n.x() + n.y() > 0.0 ? 0.0 : 1.0;
    }

    static cases_t get_cases(double p, const Vector3d& n) {
        return cases_t{cube_section_case(n, p),
                       cube_section_case(n, p - n.x()),
                       cube_section_case(n, p - n.y())};
    }

    static Vector3d get_errors(double p, const Vector3d& n, double a0, double a1, double a2) {
        return {cube_volume_fraction(n, p) - a0,
                cube_volume_fraction(n, p - n.x()) - a1,
                cube_volume_fraction(n, p - n.y()) - a2};
    }

    static double get_error(double p, const Vector3d& n, double a0, double a1, double a2) {
        return get_errors(p, n, a0, a1, a2).cwiseAbs().maxCoeff();
    }

    Info(double p, const Vector3d& n) : p(p), n(n) {
        ae = get_ae(p, n);
        cases = get_cases(p, n);

        count = -1;
        errs = {NAN, NAN, NAN};
        max_err = NAN;
    }

    void set_errors(double a0, double a1, double a2) {
        errs = get_errors(p, n, a0, a1, a2);
        max_err = errs.cwiseAbs().maxCoeff();
    }
};

// Последовательность плоскостей во время итераций
struct Series {
    std::vector<double>   ps;
    std::vector<Vector3d> ns;

    // Число итераций метода, не обязательно равно ps.size() или ns.size().
    // Если все итерации не хранятся, тогда ps и ns содержат по одному элементу.
    int count = 0;

    // Добавить плоскость
    void append(double p, const Vector3d& n) {
        ps.push_back(p);
        ns.push_back(n);
    }

    // Записать в pvd последовательность сечений
    void save(double a0, double a1, double a2);

    // Вывести данные о сходимости
    void plot(double a0, double a1, double a2);
};

// Опции для итерационных процедур
struct IterOpts {
    int max_iters = 150;
    double dumper = 0.5;

    bool series  = false;
    bool verbose = false;
};

// Геометрический способ
Series ae_geom(double a0, double a1, double a2, const IterOpts& opts) {
    Vector3d n = Vector3d::UnitZ();
    double p = 0.5 * (min(a0, a1, a2) + max(a0, a1, a2) - 1.0);

    Series series;
    for (; series.count < opts.max_iters; ++series.count) {
        // Находим точки для каждой в отдельности
        double p0 = cube_find_section(n, a0);
        double p1 = cube_find_section(n, a1);
        double p2 = cube_find_section(n, a2);

        // Три точки для разных плоскостей
        Vector3d P0 = p0 * n;
        Vector3d P1 = Vector3d::UnitX() + p1 * n;
        Vector3d P2 = Vector3d::UnitY() + p2 * n;

        // Типа усредненная точка плоскости (здесь для статистики)
        p = 0.5 * (min(p0, p1 + n.x(), p2 + n.y()) + max(p0, p1 + n.x(), p2 + n.y()));

        if (opts.series) { series.append(p, n); }

        // Строим плоскость по трём точкам
        Vector3d n2 = (P1 - P0).cross(P2 - P0).normalized();

        if (n2.z() * n.z() < 0.0) {
            double t = n.z() / (n.z() - n2.z());
            Vector3d n_t = (1.0 - t) * n + t * n2;
            n_t.z() = 0.0;
            n2 = n_t.normalized();
        }
        double step = opts.dumper;
        Vector3d n3 = step * n2 + (1.0 - step) * n;
        n = n3.normalized();

        double err_1 = std::abs(cube_volume_fraction(n, p0) - a0);
        double err_2 = std::abs(cube_volume_fraction(n, p1) - a1);
        double err_3 = std::abs(cube_volume_fraction(n, p2) - a2);

        double err = max(err_1, err_2, err_3);

        if (err < 1.0e-14) {
            break;
        }
    }

    series.append(p, n);

    return series;
}

// Первый вариант метода Ньютона. Итерации по нормали и по параметру alpha_e.
// Кое-как работает, была просто идея сделать сразу ограничения на alpha_e,
// такая идея не выгорела. На плохих условиях сходимость неизвестно куда.
Series ae_newt_1(double a0, double a1, double a2, const IterOpts& opts) {
    Vector3d n = Vector3d::UnitZ();
    double ae = 0.5 * (min(a0, a1, a2) + max(a0, a1, a2));
    double p = ae - 0.5;

    Series series;
    for (; series.count < opts.max_iters; ++series.count) {
        double p0 = +0.5 * n.x() + 0.5 * n.y() + (ae - 0.5) * n.z();
        double p1 = -0.5 * n.x() + 0.5 * n.y() + (ae - 0.5) * n.z();
        double p2 = +0.5 * n.x() - 0.5 * n.y() + (ae - 0.5) * n.z();

        p = p0;
        if (opts.series) { series.append(p, n); }

        // Значения
        double A0 = cube_volume_fraction(n, p0);
        double A1 = cube_volume_fraction(n, p1);
        double A2 = cube_volume_fraction(n, p2);

        // Градиенты
        vf_grad_t D0 = cube_volume_fraction_grad(n, p0);
        vf_grad_t D1 = cube_volume_fraction_grad(n, p1);
        vf_grad_t D2 = cube_volume_fraction_grad(n, p2);

        Vector4d F = { a0 - A0, a1 - A1, a2 - A2, 0.0 };

        //std::cout << "  A_i: " << A0 << ", " << A1 << ", " << A2 << "\n";
        //std::cout << "  RHS: " << F.transpose() << "\n";

        double err = F.cwiseAbs().maxCoeff();
        if (err < 1.0e-14) {
            break;
        }

        Matrix4d M;

        M(0, 0) = D0.dn.x() + 0.5 * D0.dp;
        M(0, 1) = D0.dn.y() + 0.5 * D0.dp;
        M(0, 2) = D0.dn.z() + (ae - 0.5) * D0.dp;
        M(0, 3) = n.z() * D0.dp;

        M(1, 0) = D1.dn.x() - 0.5 * D1.dp;
        M(1, 1) = D1.dn.y() + 0.5 * D1.dp;
        M(1, 2) = D1.dn.z() + (ae - 0.5) * D1.dp;
        M(1, 3) = n.z() * D1.dp;

        M(2, 0) = D2.dn.x() + 0.5 * D2.dp;
        M(2, 1) = D2.dn.y() - 0.5 * D2.dp;
        M(2, 2) = D2.dn.z() + (ae - 0.5) * D2.dp;
        M(2, 3) = n.z() * D2.dp;

        M(3, 0) = n.x();
        M(3, 1) = n.y();
        M(3, 2) = n.z();
        M(3, 3) = 0.0;

        Vector4d delta = M.inverse() * F;

        //std::cout << "delta: " << M << std::endl;
        //std::cout << "Det: " << M.determinant() << "\n";

        if constexpr (false) {
            // Сделать ae в интервале [0, 1], не работает
            double C = 1.0;
            if (ae + delta[3] > 1.0) { C = std::abs(1.0 - ae) / delta[3]; }
            if (ae + delta[3] < 0.0) { C = std::abs(ae) / delta[3]; }
            delta *= C;
        }

        if (delta.hasNaN()) {
            break;
        }

        n.x() += delta[0];
        n.y() += delta[1];
        n.z() += delta[2];
        n.normalize();

        ae += delta[3];

        //std::cout << "curr_ae: " << ae << "\n";

        if (std::abs(delta[3]) < 1.0e-14) {
            break;
        }

        if (n.hasNaN() || std::isnan(ae)) {
            break;
        }
    }

    if (n.hasNaN()) {
        std::cout << a0 << " " << a1 << " " << a2 << std::endl;
        throw std::runtime_error("Nan detected");
    }
    if (std::isnan(ae)) {
        n.z() = 0.0;
        n.normalize();
    }

    series.append(p, n);
    return series;
}

// Не совсем Ньютон, итерации по двум уравнениям, при этом p0 не дифференцируется
// В отличие от предыдущего легко сделать дампер, тогда даже при отсутствии решения
// сходится к чему-то и фиксируется
// С дампером 0.95 сходится практически везде за < 10 итераций.
// Опять же есть области с очень плохой сходимостью.
// При уменьшении дампера до 0.5 практически всё сходится, но среднее число
// итераций даже в простых случаях, очевидно, возрастает.
Series ae_newt_2(double a0, double a1, double a2, const IterOpts& opts) {
    Vector3d n = Vector3d::UnitZ();
    double p = 0.5 * (min(a0, a1, a2) + max(a0, a1, a2) - 1.0);

    // Это как первая итерация в геометрическом методе, работает отлично!
    // Для простого случая сразу дает ответ.
    if constexpr (true) {
        Vector3d ex = Vector3d::UnitX();
        Vector3d ey = Vector3d::UnitY();
        Vector3d ez = Vector3d::UnitZ();

        // Три точки для разных плоскостей
        Vector3d P0 = (a0 - 0.5) * ez;
        Vector3d P1 = ex + (a1 - 0.5) * ez;
        Vector3d P2 = ey + (a2 - 0.5) * ez;

        // Строим плоскость по трём точкам
        n = (P1 - P0).cross(P2 - P0).normalized();

        p = P0.dot(n);
    }

    Series series;
    for (; series.count < opts.max_iters; ++series.count) {
        double p0 = cube_find_section(n, a0);

        p = p0;
        if (opts.series) {
            series.append(p, n);
        }

        double A1 = cube_volume_fraction(n, p0 - n.x());
        double A2 = cube_volume_fraction(n, p0 - n.y());

        vf_grad_t D1 = cube_volume_fraction_grad(n, p0 - n.x());
        vf_grad_t D2 = cube_volume_fraction_grad(n, p0 - n.y());

        Vector3d F = { a1 - A1, a2 - A2, 0.0 };

        double err = F.cwiseAbs().maxCoeff();
        if (err < 1.0e-14) {
            break;
        }

        Matrix3d M;

        M(0, 0) = D1.dn.x() - D1.dp;
        M(0, 1) = D1.dn.y();
        M(0, 2) = D1.dn.z();

        M(1, 0) = D2.dn.x();
        M(1, 1) = D2.dn.y() - D2.dp;
        M(1, 2) = D2.dn.z();

        M(2, 0) = n.x();
        M(2, 1) = n.y();
        M(2, 2) = n.z();

        Vector3d delta = M.inverse() * F;

        Vector3d n_new = (n + delta).normalized();

        // Дампер
        if (n_new.z() * n.z() < 0.0) {
            double t = n.z() / (n.z() - n_new.z());
            Vector3d n_temp = (1.0 - t) * n + t * n_new;
            n_temp.z() = 0.0;
            n_new = n_temp.normalized();
        }
        double step = opts.dumper;
        n_new = (step * n_new + (1.0 - step) * n).normalized();

        err = std::abs(n.dot(n_new) - 1.0);

        // Может в этом смысле сойтись к слабому решению
        if (err < 1.0e-4) { break; }

        n = n_new;
    }

    series.append(p, n);
    return series;
}

// Почти как предыдущее, но решается система трёх + одно уравнений.
// В данном случае ищется нормаль и параметр p.
// Считает очень быстро, и сходится как будто бы везде при dumper = 0.95.
// Понятно, что среднее число итераций < 10, но есть тяжелая сходимость
// на границах куба.
// Наверное, лучший из всех вариантов.
// Всего 7 итераций не считая тех областей, где сходимости нет.
// Причем можно сразу ставить дампер равный 1.0.
// Эхх, вот если бы я только знал области предельные, чтобы просто сразу их исключать.
Series ae_newt_3(double a0, double a1, double a2, const IterOpts& opts) {
    a0 = between(a0, 0.0, 1.0);
    a1 = between(a1, 0.0, 1.0);
    a2 = between(a2, 0.0, 1.0);

    Vector3d n = Vector3d::UnitZ();
    double p = 0.5 * (min(a0, a1, a2) + max(a0, a1, a2) - 1.0);

    // Это как первая итерация в геометрическом методе, работает отлично!
    // Для простого случая сразу дает ответ.
    if constexpr (true) {
        Vector3d ex = Vector3d::UnitX();
        Vector3d ey = Vector3d::UnitY();
        Vector3d ez = Vector3d::UnitZ();

        // Три точки для разных плоскостей
        Vector3d P0 = (a0 - 0.5) * ez;
        Vector3d P1 = ex + (a1 - 0.5) * ez;
        Vector3d P2 = ey + (a2 - 0.5) * ez;

        // Строим плоскость по трём точкам
        n = (P1 - P0).cross(P2 - P0).normalized();

        p = P0.dot(n);
    }


    const double eps = 1.0e-15;

    int less_count = 0;
    int greater_count = 0;

    Series series;
    for (; series.count < opts.max_iters; ++series.count) {
        if (opts.series) { series.append(p, n); }

        double A0 = cube_volume_fraction(n, p);
        double A1 = cube_volume_fraction(n, p - n.x());
        double A2 = cube_volume_fraction(n, p - n.y());

        Vector4d F = { a0 - A0, a1 - A1, a2 - A2, 0.0 };

        // Максимум невязки
        double err_1 = F.cwiseAbs().maxCoeff();
        if (err_1 < eps) { if (opts.verbose) { std::cout << "break 1\n"; } break; }

        vf_grad_t D0 = cube_volume_fraction_grad(n, p);
        vf_grad_t D1 = cube_volume_fraction_grad(n, p - n.x());
        vf_grad_t D2 = cube_volume_fraction_grad(n, p - n.y());

        Matrix4d M;
        M << D0.dn.x(),         D0.dn.y(),         D0.dn.z(), D0.dp,
             D1.dn.x() - D1.dp, D1.dn.y(),         D1.dn.z(), D1.dp,
             D2.dn.x(),         D2.dn.y() - D2.dp, D2.dn.z(), D2.dp,
             n.x(),             n.y(),             n.z(),     0.0;

        Vector4d delta = M.inverse() * F;

        // Максимум приращения
        double err_2 = delta.cwiseAbs().maxCoeff();
        if (err_2 < eps) { if (opts.verbose) { std::cout << "break 2\n"; } break; }

        double   delta_p = delta[3];
        Vector3d delta_n = Vector3d{delta[0], delta[1], delta[2]};

        Vector3d n_new = (n + delta_n).normalized();

        // Дампер
        if (n_new.z() * n.z() < 0.0) {
            double t = n.z() / (n.z() - n_new.z());
            Vector3d n_temp = (1.0 - t) * n + t * n_new;
            n_temp.z() = 0.0;
            n_new = n_temp.normalized();
        }
        double step = opts.dumper;
        n_new = (step * n_new + (1.0 - step) * n).normalized();

        // Может в этом смысле сойтись к слабому решению
        //double err_3 = std::abs(n.dot(n_new) - 1.0);
        //if (err_3 < eps) { if (opts.verbose) { std::cout << "break 3\n"; } break; }

        Vector3d x = delta_n.normalized();
        Vector3d y = (n_new - n).normalized();
        double C = x.dot(y);
        p += between(C, 0.0, 1.0) * delta_p;


        double new_ae = Info::get_ae(p, n);

        // Выход по достижению точности по alpha_e
        //double err_4 = std::abs(new_ae - ae);
        //if (err_4 < eps) {  if (opts.verbose) { std::cout << "break 4\n"; } break; }

        // Выход если застопорились вне диапазона
        if (new_ae <= eps) {
            ++less_count;
            if (less_count > 2) {
                break;
            }
        }
        else {
            less_count = 0;
        }

        // Выход если застопорились вне диапазона
        if (new_ae >= 1.0 - eps) {
            ++greater_count;
            if (greater_count > 2) {
                break;
            }
        }
        else {
            greater_count = 0;
        }

        n = n_new;
    }

    series.append(p, n);
    return series;
}

// Решаем как задачу оптимизации
Series ae_optimize(double a0, double a1, double a2, const IterOpts& opts);


// Построить функцию alpha_e(a0, a1, a2);
void ae_full(const std::string& filename) {
    Cuboid gen = Cuboid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    gen.set_sizes(10, 10, 10);

    // Создать сетку
    EuMesh mesh(gen);
    auto ae = mesh.add<double>("ae");
    auto count = mesh.add<int>("count");
    auto err = mesh.add<double>("err");
    auto cs1 = mesh.add<int>("cs1");
    auto cs2 = mesh.add<int>("cs2");
    auto cs3 = mesh.add<int>("cs3");
    auto type = mesh.add<int>("type");

    mesh.set_max_level(5);

    std::map<cases_t, int> all_cases;

    auto setup_values = [&]() {
        mesh.for_each([&](EuCell& cell) {
            double a0 = cell.x();
            double a1 = cell.y();
            double a2 = cell.z();

            for (auto face:cell.faces()) {
                if (face.is_boundary()) {
                    a0 = face.x();
                    a1 = face.y();
                    a2 = face.z();
                    break;
                }
            }

            auto res = ae_newt_3(a0, a1, a2, {.max_iters=150, .dumper=0.95});

            double   p = res.ps.back();
            Vector3d n = res.ns.back();

            Info out(p, n);
            out.count = res.count;
            out.set_errors(a0, a1, a2);

            cell[ae]    = out.ae;
            cell[count] = out.count;
            cell[err]   = out.max_err;
            cell[cs1]   = out.cases.c1;
            cell[cs2]   = out.cases.c2;
            cell[cs3]   = out.cases.c3;

            all_cases[out.cases] = 0;
        });
    };

    auto setup_flags = [&]() {
        mesh.for_each([&](EuCell& cell) {
            cell.set_flag(-1);
            for (auto face: cell.faces()) {
                auto neib = face.neib();
                if (std::abs(cell[ae] - neib[ae]) > 0.08) {
                    cell.set_flag(1);
                    break;
                }
            }
        });
    };

    Stopwatch elapsed(true);
    for (int i = 0; i < mesh.max_level() + 3; ++i) {
        std::cout << "Level " << i << "\n";
        setup_values();
        setup_flags();
        mesh.refine();
    }
    std::cout << "Level " << mesh.max_level() << "\n";
    setup_values();
    elapsed.stop();

    std::cout << "Cases number: " << all_cases.size() << "\n";
    std::cout << "Compute ae (elapsed): " << elapsed.extended_time() << "\n";

    int counter = 0;
    for (auto& c: all_cases) {
        c.second = counter;
        ++counter;
    }

    // Пронумеровать случаи единственным индексом
    mesh.for_each([&](EuCell& cell) {
        cases_t cases(cell[cs1], cell[cs2], cell[cs3]);
        cell[type] = all_cases[cases];
    });

    Variables vars = {"center", "level"};
    vars.append("ae", ae);
    vars.append("count", count);
    vars.append("err",   err);
    vars.append("case1", cs1);
    vars.append("case2", cs2);
    vars.append("case3", cs3);
    vars.append("type", type);

    VtuFile::save(filename, mesh, vars);
}

int main() {
    //save_origin();

    double a0 = 0.0026;
    double a1 = 0.25;
    double a2 = 0.915;

    Series res = ae_optimize(a0, a1, a2, {.max_iters=150, .dumper=0.9, .series=true, .verbose=true});
    res.save(a0, a1, a2);
    res.plot(a0, a1, a2);

    //calc_alpha_edge("out/ae.vtu");

    return 0;
}

struct out_t {
    double   val  = 0.0;
    Vector3d grad = Vector3d::Zero();
};


// Функция для минимизации и её градиент
out_t target_func(const Vector3d& n, double a0,
    const std::vector<double>& as,
    const std::vector<Vector3d>& cs,
    const std::vector<double>& ws) {

    if (cs.size() != as.size() || ws.size() != as.size()) {
        throw std::invalid_argument("Array sizes mismatch");
    }

    double h = 1.0e-5;
    Vector3d ex = Vector3d::UnitX();
    Vector3d ey = Vector3d::UnitY();
    Vector3d ez = Vector3d::UnitZ();

    // Как будто бы изи
    out_t out{.val = 0.0, .grad = Vector3d::Zero()};
    double p0 = cube_find_section(n, a0);
    for (int i = 0; i < int(as.size()); ++i) {
        double a_i  = cube_volume_fraction     (n, p0 - cs[i].dot(n));
        auto grad_i = cube_volume_fraction_grad(n, p0 - cs[i].dot(n));

        out.val += ws[i] * std::pow(a_i - as[i], 2);
        out.grad += 2.0 * ws[i] * (a_i - as[i]) * (grad_i.dn - grad_i.dp * cs[i]);
    }

    debug_code {
        Vector3d ograd = Vector3d::Zero();
        for (int i = 0; i < int(as.size()); ++i) {
            double a_p, a_m;

            a_p  = cube_volume_fraction(n + h * ex, p0 - cs[i].dot(n + h * ex));
            a_m  = cube_volume_fraction(n - h * ex, p0 - cs[i].dot(n - h * ex));
            ograd.x() += ws[i] * (std::pow(a_p - as[i], 2) - std::pow(a_m - as[i], 2)) / (2 * h);

            a_p  = cube_volume_fraction(n + h * ey, p0 - cs[i].dot(n + h * ey));
            a_m  = cube_volume_fraction(n - h * ey, p0 - cs[i].dot(n - h * ey));
            ograd.y() += ws[i] * (std::pow(a_p - as[i], 2) - std::pow(a_m - as[i], 2)) / (2 * h);

            a_p  = cube_volume_fraction(n + h * ez, p0 - cs[i].dot(n + h * ez));
            a_m  = cube_volume_fraction(n - h * ez, p0 - cs[i].dot(n - h * ez));
            ograd.z() += ws[i] * (std::pow(a_p - as[i], 2) - std::pow(a_m - as[i], 2)) / (2 * h);
        }

        if ((out.grad - ograd).norm() > 1.0e-5) {
            std::cout << out.grad.transpose() << "  " << ograd.transpose() << std::endl;
            throw std::invalid_argument("Output gradient");
        }
    }

    return out;
}

// Проекция вектора на касательную плоскость к сфере
Vector3d sphere_tangent(const Vector3d& p, const Vector3d& vec) {
    return vec - vec.dot(p) * p;
}

// Движение по геодезической
Vector3d exp_map(const Vector3d& p, const Vector3d& dir, double step) {
    double norm = dir.norm();
    if (norm < 1.0e-14) { return p; }

    return p * std::cos(step * norm) + (dir / norm) * std::sin(step * norm);
}

Series ae_optimize(double a0, double a1, double a2, const IterOpts& opts) {
    std::vector as = {a1, a2};
    std::vector ws = {1.0, 1.0};
    std::vector<Vector3d> cs = {Vector3d::UnitX(), Vector3d::UnitY() };

    plot_target(a0, as, cs, ws);

    // Так, ну а здесь пытаемся честно оптимизировать

    auto func = [a0, as, cs, ws](Vector3d& n) -> out_t {
        return target_func(n, a0, as, cs, ws);
    };

    // Делаем итерации, градиентный спуск

    // Начальное приближение
    Vector3d n = Vector3d::UnitZ();
    double p = 0.5 * (min(a0, a1, a2) + max(a0, a1, a2) - 1.0);

    // Это как первая итерация в геометрическом методе, работает отлично!
    // Для простого случая сразу дает ответ.
    if constexpr (true) {
        Vector3d ex = Vector3d::UnitX();
        Vector3d ey = Vector3d::UnitY();
        Vector3d ez = Vector3d::UnitZ();

        // Три точки для разных плоскостей
        Vector3d P0 = (a0 - 0.5) * ez;
        Vector3d P1 = ex + (a1 - 0.5) * ez;
        Vector3d P2 = ey + (a2 - 0.5) * ez;

        // Строим плоскость по трём точкам
        n = (P1 - P0).cross(P2 - P0).normalized();

        p = P0.dot(n);
    }

    int max_iter = 1000;
    double ntol = 1.0e-12; // Погрешность по разнице нормалей
    double stol = 1.0e-12; // Погрешность по размеру шага
    double ftol = 1.0e-12; // Погрешность по значению функции
    double gtol = 1.0e-12; // Погрешность по норме градиента

    double init_step = 0.25 * M_PI;

    double n_err = NAN;
    double f_err = NAN;

    Vector3d curr_n = n;
    Vector3d prev_n = n;
    double prev_f = std::numeric_limits<double>::max();

    Series series;
    for (; series.count < opts.max_iters; ++series.count) {
        if (opts.series) {
            series.append(cube_find_section(curr_n, a0), curr_n);
        }

        out_t res = func(curr_n);

        double curr_f = res.val;

        // Направление градиентного спуска
        Vector3d desc_grad = -sphere_tangent(curr_n, res.grad);


        if (opts.verbose) {
            std::cout << "  " << series.count << ")\tn: " << curr_n.transpose() << ";\tf: " << curr_f << ";\tgrad: " << desc_grad.transpose() << "\n";
        }

        // Backtracking line search

        // Условие Армихо: f(x + a * dir) <= f(x) + c1 * a * grad f * dir
        double dir_deriv = -desc_grad.dot(desc_grad);

        const double rho = 0.5; // Коэффициент уменьшения шага
        const double c1 = 1.0e-4;

        int ls_iter = 0;
        Vector3d new_n;
        out_t new_res;

        double alpha = init_step;
        while (true) {

            new_n = exp_map(curr_n, desc_grad, alpha);
            new_res = func(new_n);

            if (new_res.val <= curr_f + c1 * alpha * dir_deriv) {
                break;
            }

            alpha *= rho;
            ++ls_iter;
            if (ls_iter > 20) {
                alpha = 0.0;
                break;
            }
        }

        if (alpha > 1.0e-12) {
            curr_n = new_n;
        }
        else {
            std::cout << "Line search failed\n";
            break;
        }

        // Градиент достиг нулевого значения (экстремум)
        double g_err = res.grad.norm();
        if (g_err < gtol) {
            std::cout << " grad tol: " << g_err << std::endl;
        //    break;
        }

        // Функция перестала изменяться
        f_err = std::abs(curr_f - prev_f);
        if (f_err < ftol) {
            std::cout << " func tol: " << f_err << std::endl;
            break;
        }

        // Нормаль перестала изменяться
        n_err = std::abs(curr_n.dot(prev_n) - 1.0);
        if (n_err < ntol) {
            std::cout << "  norm tol: " << n_err << std::endl;
            break;
        }

        prev_n = curr_n;
        prev_f = curr_f;
    }

    series.append(cube_find_section(curr_n, a0), curr_n);
    return series;
}

cases_t::cases_t(int c1_in, int c2_in, int c3_in)
    : c1(c1_in), c2(c2_in), c3(c3_in) {

    // Простейший случай (есть чистая ячейка)
    if (c1 == 0 || c2 == 0 || c3 == 0) {
        c1 = c2 = c3 = -13;
    }

    // Инвертируем относительно основной
    if (c1 < 0) { c1 *= -1; c2 *= -1; c3 *= -1; }

    // c2, c3 эквивалентны, упорядочим, чтобы снизить размерность
    if (c2 > c3) std::swap(c2, c3);
}

bool cases_t::operator!=(const cases_t &other) const {
    return c1 != other.c1 || c2 != other.c2 || c3 != other.c3;
}

bool cases_t::operator<(const cases_t &other) const {
    if (c1 < other.c1) { return true; }
    if (c1 > other.c1) { return false; }

    if (c2 < other.c2) { return true; }
    if (c2 > other.c2) { return false; }

    return c3 < other.c3;
}

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

    VtuFile::save("out/orig.vtu", orig, {}, false, true);
}

Vector3d plot_target(double a0,
    const std::vector<double>& as,
    const std::vector<Vector3d>& cs,
    const std::vector<double>& ws) {

    // Создать случайные точки на сфере
    int M = 20;
    int N = 6 * M * M;

    auto row = np::linspace(-1.0 + 1.0 / M, 1.0 - 1.0 / M, M);

    auto xs = np::zeros(N);
    auto ys = np::zeros(N);
    auto zs = np::zeros(N);

    // Нормализуем
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            std::array<Vector3d, 6> ns = {
                Vector3d{-1.0, row[i], row[j]},
                Vector3d{+1.0, row[i], row[j]},
                Vector3d{row[i], -1.0, row[j]},
                Vector3d{row[i], +1.0, row[j]},
                Vector3d{row[i], row[j], -1.0},
                Vector3d{row[i], row[j], +1.0}
            };
            for (auto& n: ns) n.normalize();

            for (int k = 0; k < 6; ++k) {
                int m = 6 * (i * M + j) + k;
                xs[m] = ns[k].x();
                ys[m] = ns[k].y();
                zs[m] = ns[k].z();
            }
        }
    }

    // Значения функции
    auto fs = np::zeros(N);

    for (int i = 0; i < N; ++i) {
        Vector3d n = {xs[i], ys[i], zs[i]};
        fs[i] = target_func(n, a0, as, cs, ws).val;
    }


    utils::pyplot plt;

    plt.figure();
    plt.scatter3D(xs, ys, zs, {.c=fs, .colorbar=true});
    plt.tight_layout();
    plt.show();

    // Наилучший результат в наборе
    size_t idx = std::min_element(fs.begin(), fs.end()) - fs.begin();
    return Vector3d{xs[idx], ys[idx], std::abs(zs[idx])};
}

void Series::save(double a0, double a1, double a2) {
    auto poly1 = Polyhedron::Cube();
    auto poly2 = Polyhedron::Cube();
    auto poly3 = Polyhedron::Cube();

    poly2.move(Vector3d::UnitX());
    poly3.move(Vector3d::UnitY());

    PvdFile pvd("clipped", "out");
    pvd.hex_only = false;
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

void Series::plot(double a0, double a1, double a2) {
    auto steps = np::arange<double>(ps.size());
    std::vector<double> errs(ps.size());
    std::vector<double> a_es(ps.size());

    for (size_t i = 0; i < ps.size(); ++i) {
        errs[i] = Info::get_error(ps[i], ns[i], a0, a1, a2);
        a_es[i] = Info::get_ae(ps[i], ns[i]);
    }

    std::cout << "a0, a1, a2: " << a0 << ", " << a1 << ", " << a2 << "\n";
    std::cout << "  p, n: " << ps.back() << ", " << ns.back().transpose() << "\n";
    std::cout << "  ae: " << a_es.back() << "\n";

    // Отклонение от последнего значения
    for (size_t i = 0; i < ps.size(); ++i) {
        a_es[i] = std::abs(a_es[i] - a_es.back());
    }

    utils::pyplot plt;

    plt.figure({.figsize={9.0, 4.0}, .dpi=170});

    plt.subplot(1, 2, 1);
    plt.title("Погрешность от итерации");
    plt.grid(true);
    plt.semilogy(steps, errs, "k");

    plt.subplot(1, 2, 2);
    plt.title("Погрешность $\\alpha_e$");
    plt.grid(true);
    plt.semilogy(steps, a_es, "k");

    plt.tight_layout();
    plt.show();
}