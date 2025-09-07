// Проверка производных на AMR-шаблонах

#include <iomanip>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

#include <zephyr/utils/matplotlib.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;

namespace plt = zephyr::utils::matplotlib;


struct _U_ {
    // Метка целевой ячейки
    int target;

    // Значение функции
    double u;

    // Точный градиент
    Vector2d grad;

    // Гаусс
    Vector2d grad_g;
    Vector2d err_g;

    // Оригинал МНК
    Vector2d grad_o;
    Vector2d err_o;

    // Новый МНК
    Vector2d grad_n;
    Vector2d err_n;

    double err_gauss() const { return err_g.norm(); }
    double err_orig()  const { return err_o.norm(); }
    double err_new()   const { return err_n.norm(); }
};

_U_ U;

// Переменные для сохранения
double get_interesting(EuCell& cell)  { return double(cell(U).target); }
double get_u(EuCell& cell) { return cell(U).u; }
double get_ux(EuCell& cell) { return cell(U).grad.x(); }
double get_uy(EuCell& cell) { return cell(U).grad.y(); }
double get_du_dx_o(EuCell& cell) { return cell(U).grad_o.x(); }
double get_du_dy_o(EuCell& cell) { return cell(U).grad_o.y(); }
double get_err_x_o(EuCell& cell) { return cell(U).err_o.x(); }
double get_err_y_o(EuCell& cell) { return cell(U).err_o.y(); }
double get_du_dx_n(EuCell& cell) { return cell(U).grad_n.x(); }
double get_du_dy_n(EuCell& cell) { return cell(U).grad_n.y(); }
double get_err_x_n(EuCell& cell) { return cell(U).err_n.x(); }
double get_err_y_n(EuCell& cell) { return cell(U).err_n.y(); }
double get_du_dx_g(EuCell& cell) { return cell(U).grad_g.x(); }
double get_du_dy_g(EuCell& cell) { return cell(U).grad_g.y(); }
double get_err_x_g(EuCell& cell) { return cell(U).err_g.x(); }
double get_err_y_g(EuCell& cell) { return cell(U).err_g.y(); }


// Расчет градиента методом Гаусса
void gauss(EuCell &cell) {
    if (!cell(U).target) {
        cell(U).grad_g = Vector2d::Zero();
        cell(U).err_g  = Vector2d::Zero();
        return;
    }

    double zc = cell(U).u;

    Vector2d F = Vector2d::Zero();
    for (auto &face: cell.faces()) {
        auto neib = face.neib();

        double zn = neib(U).u;

        Vector3d S = face.normal() * face.area();

        double d1 = std::abs((face.center() - cell.center()).dot(face.normal()));
        double d2 = std::abs((face.center() - neib.center()).dot(face.normal()));
        double a1 = d2 / (d1 + d2);
        double a2 = d1 / (d1 + d2);

        double zf = a1 * zc + a2 * zn;

        F.x() += zf * S.x();
        F.y() += zf * S.y();
    }
    F /= cell.volume();

    cell(U).grad_g = F;
    cell(U).err_g = (cell(U).grad_g - cell(U).grad).cwiseAbs();
}

// Расчет градиента старым МНК
void LSM_old(EuCell &cell) {
    if (!cell(U).target) {
        cell(U).grad_o = Vector2d::Zero();
        cell(U).err_o  = Vector2d::Zero();
        return;
    }

    double zc = cell(U).u;

    Vector3d F = Vector3d::Zero();
    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        auto neib = face.neib();

        double zn = neib(U).u;

        Vector3d dr = neib.center() - cell.center();

        double weight = 1.0 / std::pow(dr.norm(), 2);

        A += weight * dr * dr.transpose();

        F += (weight * (zn - zc)) * dr;
    }

    Vector3d grad = A.colPivHouseholderQr().solve(F);
    cell(U).grad_o = {grad.x(), grad.y()};
    cell(U).err_o = (cell(U).grad_o - cell(U).grad).cwiseAbs();
}

// Расчет градиента новым МНК
void LSM_new(EuCell &cell) {
    if (!cell(U).target) {
        cell(U).grad_n = Vector2d::Zero();
        cell(U).err_n  = Vector2d::Zero();
        return;
    }

    double zc = cell(U).u;

    Vector3d F = Vector3d::Zero();
    Matrix3d A = Matrix3d::Zero();
    for (auto &face: cell.faces()) {
        auto neib = face.neib();

        double zn = neib(U).u;

        Vector3d dr = neib.center() - cell.center();

        Vector3d drn = dr.dot(face.normal()) * face.normal();

        double weight = face.area() / std::pow(drn.norm(), 3);

        A += weight * drn * dr.transpose();

        F += (weight * (zn - zc)) * drn;
    }

    Vector3d grad = A.colPivHouseholderQr().solve(F);
    cell(U).grad_n = {grad.x(), grad.y()};
    cell(U).err_n = (cell(U).grad_n - cell(U).grad).cwiseAbs();
}

// Количество двумерных шаблонов
constexpr int n_templates() { return 11; }

// Все двумерные шаблоны
// num -- номер шаблона
// H -- линейный размер целевой ячейки
EuMesh get_template(int num, double H) {
    if (0 <= num && num <= 5) {
        Rectangle rect(-1.5 * H, 1.5 * H, -1.5 * H, 1.5 * H);
        rect.set_nx(3);

        EuMesh mesh(rect);
        mesh.add<_U_>("U");
        mesh.set_max_level(1);

        if (num == 0) {
            mesh.set_max_level(0);
        } else if (num == 1) {
            Vector3d c1 = {-H, 0.0, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        } else if (num == 2) {
            Vector3d c1 = {-H, 0.0, 0.0};
            Vector3d c2 = {0.0, +H, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H ||
                    (cell.center() - c2).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        } else if (num == 3) {
            Vector3d c1 = {-H, 0.0, 0.0};
            Vector3d c2 = {+H, 0.0, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H ||
                    (cell.center() - c2).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        } else if (num == 4) {
            Vector3d c1 = {-H, 0.0, 0.0};
            Vector3d c2 = {0.0, +H, 0.0};
            Vector3d c3 = {+H, 0.0, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H ||
                    (cell.center() - c2).norm() < 0.1 * H ||
                    (cell.center() - c3).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        } else {
            Vector3d c1 = {-H, 0.0, 0.0};
            Vector3d c2 = {0.0, +H, 0.0};
            Vector3d c3 = {+H, 0.0, 0.0};
            Vector3d c4 = {0.0, -H, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H ||
                    (cell.center() - c2).norm() < 0.1 * H ||
                    (cell.center() - c3).norm() < 0.1 * H ||
                    (cell.center() - c4).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        }

        mesh.refine();
        for (auto cell: mesh) {
            if (cell.center().norm() < 0.1 * H) {
                cell(U).target = 1;
            }
        }

        return mesh;
    } else {
        Rectangle rect(-2.5 * H, 1.5 * H,
                       -1.5 * H, 2.5 * H);
        rect.set_nx(2);

        EuMesh mesh(rect);
        mesh.add<_U_>("U");
        mesh.set_max_level(2);

        if (num == 6) {
            Vector3d c1 = {0.5 * H, -0.5 * H, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
        } else {
            Vector3d c1 = {0.5 * H, -0.5 * H, 0.0};
            Vector3d c2 = {0.5 * H, +1.5 * H, 0.0};
            for (auto cell: mesh) {
                if ((cell.center() - c1).norm() < 0.1 * H ||
                    (cell.center() - c2).norm() < 0.1 * H) {
                    cell.set_flag(1);
                }
            }
            mesh.refine();

            if (num == 8) {
                Vector3d x1 = {H, 0.0, 0.0};
                for (auto cell: mesh) {
                    if ((cell.center() - x1).norm() < 0.1 * H) {
                        cell.set_flag(1);
                    }
                }
            } else if (num == 9) {
                Vector3d x1 = {0.0, H, 0.0};
                for (auto cell: mesh) {
                    if ((cell.center() - x1).norm() < 0.1 * H) {
                        cell.set_flag(1);
                    }
                }
            } else if (num == 10) {
                Vector3d x1 = {H, 0.0, 0.0};
                Vector3d x2 = {0.0, H, 0.0};
                for (auto cell: mesh) {
                    if ((cell.center() - x1).norm() < 0.1 * H ||
                        (cell.center() - x2).norm() < 0.1 * H) {
                        cell.set_flag(1);
                    }
                }
            }
        }

        mesh.refine();
        for (auto cell: mesh) {
            if (cell.center().norm() < 0.1 * H) {
                cell(U).target = 1;
            }
        }

        return mesh;
    }
}

// Линейная функция
double test_linear_func(double x, double y) {
    return 1.4 * x - 0.8 * y + 3.0;
}

// Линейная функция
Vector2d test_linear_grad(double x, double y) {
    return {1.4, -0.8};
}

// Произвольная функция однородная по X
double test_x_func(double x, double y) {
    return 5.0 * std::sin(0.2 * (y + 2.0)) / std::cos(0.4);
}

// Произвольная функция однородная по X
Vector2d test_x_grad(double x, double y) {
    return {0.0, std::cos(0.2 * (y + 2.0)) / std::cos(0.4)};
}

// Произвольная функция однородная по Y
double test_y_func(double x, double y) {
    return 5.0 * std::sin(0.2 * (x + 2.0)) / std::cos(0.4);
}

// Произвольная функция однородная по Y
Vector2d test_y_grad(double x, double y) {
    return {std::cos(0.2 * (x + 2.0)) / std::cos(0.4), 0.0};
}

// Квадратичная функция однородная по X
double test_x_sq_func(double x, double y) {
    return 0.4 * y * y - 0.9 * y + 1.0;
}

// Квадратичная функция однородная по X
Vector2d test_x_sq_grad(double x, double y) {
    return {0.0, 0.8 * y - 0.9};
}

// Квадратичная функция однородная по Y
double test_y_sq_func(double x, double y) {
    return 0.4 * x * x - 0.9 * x + 1.0;
}

// Квадратичная функция однородная по Y
Vector2d test_y_sq_grad(double x, double y) {
    return {0.8 * x - 0.9, 0.0};
}

// Сложная полиномиальная функция
double test_arb_func(double x, double y) {
    double f(0.5), f_x(0.9), f_y(-1.3), f_xx(0.03), f_xy(-0.04), f_yy(-0.02);
    double f_xxx(0.07), f_xxy(0.13), f_xyy(0.11), f_yyy(-0.08);
    double f_xxxx(0.15), f_xxxy(-0.2), f_xxyy(0.25), f_xyyy(0.2), f_yyyy(-0.12);

    double x2(x * x), x3(x * x * x), x4(x * x * x * x);
    double y2(y * y), y3(y * y * y), y4(y * y * y * y);

    double res = f + f_x * x + f_y * y;
    res += (f_xx * x2 + 2 * f_xy * x * y + f_yy * y2) / 2.0;
    res += (f_xxx * x3 + 3 * f_xxy * x2 * y + 3 * f_xyy * x * y2 + f_yyy * y3) / 6.0;
    res += (f_xxxx * x4 + 4 * f_xxxy * x3 * y + 6 * f_xxyy * x2 * y2 + 4 * f_xyyy * x * y3 + f_yyyy * y4) / 24.0;

    return res;
}

// Сложная полиномиальная функция
Vector2d test_arb_grad(double x, double y) {
    double f(0.5), f_x(0.9), f_y(-1.3), f_xx(0.03), f_xy(-0.04), f_yy(-0.02);
    double f_xxx(0.07), f_xxy(0.13), f_xyy(0.11), f_yyy(-0.08);
    double f_xxxx(0.15), f_xxxy(-0.2), f_xxyy(0.25), f_xyyy(0.2), f_yyyy(-0.12);

    double x2(x * x), x3(x * x * x), x4(x * x * x * x);
    double y2(y * y), y3(y * y * y), y4(y * y * y * y);

    double df_dx = f_x;
    df_dx += (2 * f_xx * x + 2 * f_xy * y) / 2.0;
    df_dx += (3 * f_xxx * x2 + 6 * f_xxy * x * y + 3 * f_xyy * y2) / 6.0;
    df_dx += (4 * f_xxxx * x3 + 12 * f_xxxy * x2 * y + 12 * f_xxyy * x * y2 + 4 * f_xyyy * y3) / 24.0;

    double df_dy = f_y;
    df_dy += (2 * f_xy * x + 2 * f_yy * y) / 2.0;
    df_dy += (3 * f_xxy * x2 + 6 * f_xyy * x * y + 3 * f_yyy * y2) / 6.0;
    df_dy += (4 * f_xxxy * x3 + 12 * f_xxyy * x2 * y + 12 * f_xyyy * x * y2 + 4 * f_yyyy * y3) / 24.0;

    return {df_dx, df_dy};
}

// Задать на сетке все данные
void set_data(EuMesh& mesh,
              std::function<double(double, double)> func,
              std::function<Vector2d(double, double)> grad) {

    for (auto cell: mesh) {
        cell(U).u    = func(cell.x(), cell.y());
        cell(U).grad = grad(cell.x(), cell.y());
    }

    mesh.for_each(gauss);
    mesh.for_each(LSM_old);
    mesh.for_each(LSM_new);
}

int main() {
    PvdFile pvd("mesh", "output");
    pvd.variables += {"target", get_interesting};
    pvd.variables += {"u", get_u };
    pvd.variables += {"du/dx", get_ux };
    pvd.variables += {"du/dx_o", get_du_dx_o };
    pvd.variables += {"du/dx_n", get_du_dx_n };
    pvd.variables += {"du/dx_g", get_du_dx_g };
    pvd.variables += {"err_x_o", get_err_x_o };
    pvd.variables += {"err_x_n", get_err_x_n };
    pvd.variables += {"err_x_g", get_err_x_g };

    pvd.variables += {"du/dy", get_uy };
    pvd.variables += {"du/dy_o", get_du_dy_o };
    pvd.variables += {"du/dy_n", get_du_dy_n };
    pvd.variables += {"du/dy_g", get_du_dy_g };
    pvd.variables += {"err_y_o", get_err_y_o };
    pvd.variables += {"err_y_n", get_err_y_n };
    pvd.variables += {"err_y_g", get_err_y_g };

    // Записать все шаблоны в файл
    for (int k = 0; k < n_templates(); ++k) {
        EuMesh mesh = get_template(k, 1.0);
        set_data(mesh, test_linear_func, test_linear_grad);
        pvd.save(mesh, k + 1);
    }

    // Тестирование сеточной сходимости
    for (int k = 0; k < n_templates(); ++k) {
        std::vector<double> hs;
        std::vector<double> lsm_old;
        std::vector<double> lsm_new;
        std::vector<double> gauss;

        std::vector<double> acc_1st;
        std::vector<double> acc_2nd;

        for (int j = 0; j < 10; ++j) {
            double xi = 0.3;
            double h = std::pow(xi, j);
            hs.push_back(h);

            EuMesh mesh = get_template(k, h);

            //set_data(mesh, test_linear_func, test_linear_grad);
            //set_data(mesh, test_x_func, test_x_grad);
            set_data(mesh, test_y_func, test_y_grad);
            //set_data(mesh, test_x_sq_func, test_x_sq_grad);
            //set_data(mesh, test_y_sq_func, test_y_sq_grad);
            //set_data(mesh, test_arb_func, test_arb_grad);

            for (auto cell: mesh) {
                if (cell(U).target) {
                    lsm_old.push_back(cell(U).err_orig());
                    lsm_new.push_back(cell(U).err_new());
                    gauss.push_back(cell(U).err_gauss());
                    break;
                }
            }

            if (j < 1) {
                double E0min = std::min({lsm_old[0], lsm_new[0], gauss[0]});
                double E0max = std::max({lsm_old[0], lsm_new[0], gauss[0]});
                acc_1st.push_back(E0max);
                acc_2nd.push_back(E0min);
            }
            else {
                acc_1st.push_back(acc_1st[0] * std::pow(xi, j));
                acc_2nd.push_back(acc_2nd[0] * std::pow(xi, 2 * j));
            }
        }

        plt::figure_size(5.4, 3.6, 250);
        plt::figure(k + 1);

        plt::xlabel("Размер ячейки");
        plt::ylabel("Погрешность");

        plt::grid(true);

        plt::loglog(hs, acc_1st, {{"color", "black"}, {"linestyle", "dotted"}});
        plt::loglog(hs, acc_2nd, {{"color", "black"}, {"linestyle", "dotted"}});

        plt::loglog(hs, gauss, {{"color", "blue"}, {"label", "Гаусс"}, {"marker", "o"}});
        plt::loglog(hs, lsm_old, {{"color", "green"}, {"label", "МНК"}, {"marker", "o"}});
        plt::loglog(hs, lsm_new, {{"color", "orange"}, {"label", "сМНК"}, {"marker", "o"}});

        plt::legend();

        plt::tight_layout();
    }

    plt::show();

    return 0;
}