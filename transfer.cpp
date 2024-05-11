/// @file Решение задачи переноса с CRP решателем.

#include "fast.h"

#include <zephyr/math/solver/transfer.h>

/// @brief Наследуем собственный решатель от Transfer, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Transfer {
public:

    /// @brief Задаем скорость переноса
    Vector3d velocity(const Vector3d& c) const override {
        return {0.7, 0.5, 0.0};

        double rx = c.x() - 0.5;
        double ry = c.y() - 0.5;
        double phi = std::atan2(ry, rx);
        double v = 0.5 * std::sqrt(rx * rx + ry * ry);
        return {v * std::sin(phi), -v * std::cos(phi), 0.0};

         //0.40 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
    }
};

// Получим тип данных
Solver::State U = Solver::datatype();
Solver solver;
// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  { return cell(U).u1; }

double get_lvl(AmrStorage::Item& cell)  { return cell.level; }

double normal_x(AmrStorage::Item& cell) { return cell(U).n.x(); }

double normal_y(AmrStorage::Item& cell) { return cell(U).n.y(); }

double point_x(AmrStorage::Item& cell) { return cell(U).p.x(); }

double point_y(AmrStorage::Item& cell) { return cell(U).p.y(); }

double get_over(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    if (u < 0.0) {
        return u;
    }
    else if (u <= 1.0) {
        return 0.0 / 0.0;
    }
    else {
        return u - 1.0;
    }
}
double get_fr(AmrStorage::Item& cell){
    return  cell(U).flow_u[0];
}
double get_ft(AmrStorage::Item& cell){
    return  cell(U).flow_u[1];
}
double get_fl(AmrStorage::Item& cell){
    return  cell(U).flow_u[2];
}
double get_fb(AmrStorage::Item& cell){
    return  cell(U).flow_u[3];
}
double get_u2(AmrStorage::Item& cell){
    return  cell(U).u2;
}
double get_u3(AmrStorage::Item& cell){
    return  (cell(U).flow_u[3] + cell(U).flow_u[2]  + cell(U).flow_u[1] + cell(U).flow_u[0] +  cell(U).u1)/5 ;
}
double get_close(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    return std::abs(u < 0.5 ? u : 1.0 - u);
}

double get_lambda_alpha(AmrStorage::Item& cell) {
    return cell(U).lambda_alpha;
}
double get_lambda_beta(AmrStorage::Item& cell) {
    return cell(U).lambda_beta;
}
double get_sp(AmrStorage::Item& cell){
    if (cell(U).lambda_beta*cell(U).lambda_alpha >0) return 1; // 1 - usual point
    return 0; // sonic point
}

// Начальное условие в виде полосы
void setup_initial_1(EuMesh& mesh, bool exact = false) {
    // Полоса задана двумя линиями
    auto inside = [](const Vector3d& v) -> bool {
        const Vector3d p1 = {0.1, 0.0, 0.0};
        const Vector3d p2 = {0.3, 0.0, 0.0};
        const Vector3d n1 = {-1.0, 0.0, 0.0};
        const Vector3d n2 = {+1.0, 0.0, 0.0};

        return (v - p1).dot(n1) < 0.0 &&
               (v - p2).dot(n2) < 0.0;
    };

    for (auto cell: mesh) {
        if (!exact) {
            cell(U).u1 = inside(cell.center());
        }
        else {
            double vol_frac = cell.approx_vol_fraction(inside);
            if (0.0 < vol_frac && vol_frac < 1.0) {
                vol_frac = cell.volume_fraction(inside, 10000);
            }
            cell(U).u1 = vol_frac;
        }
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде квадрата
void setup_initial_2(EuMesh& mesh, bool exact = false) {
    auto inside = [](const Vector3d &v) -> bool {
        const double x_min = 0.1;
        const double x_max = 0.4;
        const double y_min = 0.1;
        const double y_max = 0.4;

        return x_min < v.x() && v.x() < x_max &&
               y_min < v.y() && v.y() < y_max;
    };
    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);
            if (!exact) {
                cell(U).u1 = inside(cell.center());
                //cell(U).flow_u[1] = cell(U).flow_u[3] = 0.5 * cell(U).u1;
                //cell(U).flow_u[0] = cell(U).flow_u[2] = 0.7 * cell(U).u1;
                for (int k = 0; k < 4; k++) {
                    cell(U).flow_u[k] = cell(U).u1;
                    cell(U).flow_u_tmp[k] = cell(U).flow_u[k];
                }
            } else {
                double vol_frac = cell.approx_vol_fraction(inside);
                if (0.0 < vol_frac && vol_frac < 1.0) {
                    vol_frac = cell.volume_fraction(inside, 10000);
                }
                cell(U).u1 = vol_frac;
            }
            cell(U).u2 = 0.0;
        }
    }

    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);
            if (solver.velocity(cell.vs<+1, 0>()).x() > 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  cell(U).u1/solver.velocity(cell.vs<+1, 0>()).x();
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  mesh(i+1,j).data(U).u1/solver.velocity(cell.vs<+1, 0>()).x();

            if (solver.velocity(cell.vs<-1, 0>()).x() > 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =  mesh(i-1,j).data(U).u1/solver.velocity(cell.vs<-1, 0>()).x();
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =   cell(U).u1/solver.velocity(cell.vs<+1, 0>()).x();

            if (solver.velocity(cell.vs<0, 1>()).y() > 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  cell(U).u1/solver.velocity(cell.vs<0, 1>()).y();
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  mesh(i,j+1).data(U).u1/solver.velocity(cell.vs<0, 1>()).y();

            if (solver.velocity(cell.vs<0, -1>()).y() > 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  mesh(i,j-1).data(U).u1/solver.velocity(cell.vs<0, -1>()).y();
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  cell(U).u1/solver.velocity(cell.vs<0, 1>()).y();
        }
    }
}

// Начальное условие в виде круга
void setup_initial_3(EuMesh& mesh, bool exact = true) {
    const double R = 0.15;
    const Vector3d c = {0.25, 0.25, 0.0};

    for (auto cell: mesh) {
        auto poly = cell.polygon();
        cell(U).u1 = poly.disk_clip_area(c, R) / cell.volume();
        for (int k = 0; k < 4; k++) {
            cell(U).flow_u[k] = cell(U).u1;
            cell(U).flow_u_tmp[k] = cell(U).flow_u[k];
        }

        cell(U).u2 = 0.0;
    }

    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);
            if (solver.velocity(cell.vs<+1, 0>()).x() > 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  cell(U).u1;
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  mesh(i+1,j).data(U).u1;

            if (solver.velocity(cell.vs<-1, 0>()).x() > 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =  mesh(i-1,j).data(U).u1;
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =   cell(U).u1;

            if (solver.velocity(cell.vs<0, 1>()).y() > 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  cell(U).u1;
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  mesh(i,j+1).data(U).u1;

            if (solver.velocity(cell.vs<0, -1>()).y() > 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  mesh(i,j-1).data(U).u1;
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  cell(U).u1;
        }
    }


}

// Объем тела
double volume(EuMesh& cells) {
    double sum = 0.0;
    for (auto cell: cells) {
        if (cell.volume() >= 0)
            sum += cell.volume() * cell(U).u1;
    }
    return sum;
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");
    PvdFile pvd_scheme("scheme", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"u2",  get_u2};
    pvd.variables += {"u3",  get_u3};
    pvd.variables += {"fr",  get_fr};
    pvd.variables += {"ft",  get_ft};
    pvd.variables += {"fl",  get_fl};
    pvd.variables += {"fb",  get_fb};
    //pvd.variables += {"lvl", get_lvl};
    //pvd.variables += {"n.x", normal_x};
    //pvd.variables += {"n.y", normal_y};
    //pvd.variables += {"p.x", point_x};
    //pvd.variables += {"p.y", point_y};
    //pvd.variables += {"over", get_over};
    //pvd.variables += {"close", get_close};
    //pvd.variables += {"l_a", get_lambda_alpha};
    //pvd.variables += {"l_b", get_lambda_beta};
    //pvd.variables += {"sp", get_sp};

    pvd_scheme.variables += {"u",  get_u};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 1.0, false);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

    // Создать решатель
    solver.set_CFL(0.5);
    solver.set_version(4);
    solver.dir_splitting(false);
    solver.set_mnt(true);

    // Создать сетку
    EuMesh mesh(U, &rect);

    // Настраиваем адаптацию
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Адаптация под начальные данные
    int n_init_loops = mesh.is_adaptive() ? mesh.max_level() + 2 : 0;
    for (int k = n_init_loops; k >= 0; --k) {
        setup_initial_3(mesh, 0);
        solver.set_flags(mesh);
        mesh.refine();
    }
    solver.update_interface(mesh);
    solver.prep_ver4(mesh);
    double init_volume = volume(mesh);
    std::cout << "Начальный объем: " << init_volume << "\n";

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time < 0.8 && n_step < 50) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh);
            std::cout << "\tLoss: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh.locals(), curr_time);

            AmrStorage body = solver.body(mesh);
            pvd_body.save(body, curr_time);

            AmrStorage scheme = solver.scheme(mesh);
            pvd_scheme.save(scheme, curr_time);

            next_write += 0.0;
        }

        // Шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        //solver.set_flags(mesh);

        // Адаптировать сетку
        //mesh.refine();

        n_step += 1;
        curr_time += solver.dt();
    }

    return 0;
}