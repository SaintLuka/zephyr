#include <cmath>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;
using generator::Rectangle;
using zephyr::utils::Stopwatch;
using zephyr::utils::threads;

double left_x(const EuCell& cell) {
    return 0.0;
}

double right_x(const EuCell& cell) {
  	double y = cell.y();
    return std::sin(5) * std::cos(M_PI * std::pow(y, 2) / 2);
}

double left_y(const EuCell& cell) {
  	double x = cell.x();
    return std::sin(5 * x);
}

double right_y(const EuCell& cell) {
    double x = cell.x();
  	return -std::sin(5 * x);
}

double func(const EuCell& cell) {
  	double x = cell.x();
    double y = cell.y();
    return std::sin(5 * x) * (
        (25 + std::pow(M_PI * y, 2)) * std::cos(M_PI * std::pow(y, 2) / 2) +
        M_PI * std::sin(M_PI * std::pow(y, 2) / 2)
    );
}

template <typename T = double>
using operator_t = std::function<T(const EuCell&, Storable<T>)>;

operator_t<> OpL = [](const EuCell& cell, Storable<double> u) {
	double res = 0.0;
    for (auto face: cell.faces()) {
    	if (face.is_boundary()) {
          	double coef = face.area() /
                        (face.center() - cell.center()).norm() /
                        cell.volume();
    		res += cell[u] * coef;
    	} else {
          	double coef = face.area() /
                        (face.neib_center() - cell.center()).norm() /
                        cell.volume();
    		res += (cell[u] - face.neib(u)) * coef;
    	}
    }
    return res;
};

template <typename T = double>
using rhs_t = std::function<T(const EuCell&)>;

rhs_t<> RHS = [](const EuCell& cell) {
    double res = func(cell);
    for (auto face: cell.faces()) {
    	if (face.is_boundary()) {
          	double coef = face.area() /
                        (face.center() - cell.center()).norm() /
                        cell.volume();
    		if (face.normal().x() < 0) {
               	res += left_x(cell) * coef;
            } else if (face.normal().x() > 0) {
              	res += right_x(cell) * coef;
            } else if (face.normal().y() < 0) {
              	res += left_y(cell) * coef;
            } else if (face.normal().y() > 0) {
            	res += right_y(cell) * coef;
            }
    	}
	}
    return res;
};

double aim_func(const EuCell& cell) {
  	double x = cell.x(), y = cell.y();
    return std::sin(5 * x) * std::cos(M_PI * std::pow(y, 2) / 2);
}

template <typename T = double>
struct options_t {
   int max_iters = 20000;
   double rtol = 1e-6;
   bool use_initial = false;
   std::function<void(EuMesh&, Storable<T>)> callback = nullptr;
};

template <typename T = double>
void solver(EuMesh& mesh, Storable<T> u, const operator_t<T>& OpL, const rhs_t<T>& RHS, const options_t<T>& opts) {
    //Использовать начальное приближение?
	if (!opts.use_initial) {
    	mesh.for_each([u](EuCell cell) {
            cell[u] = 0.0;
    	});
	}

    //Переменные для хранения на сетке. Невязка и правая часть (без невязки никак не обойтись(?))
    Storable<T> r_cur = mesh.add<T>("r_cur");
    Storable<T> rhs = mesh.add<T>("rhs");

    //Инициализация правой части
    mesh.for_each([rhs, RHS](EuCell cell) {
      	cell[rhs] = RHS(cell);
    });

    int counter = 0;
    T tau, max_err;

    //Начало цикла
    do {
		T vdot_r_cur = 0, vdot_lr_cur = 0;

        //Обновление невязки и подсчет скалярных произведений (r_cur, r_cur) и (L[r_cur], r_cur)
        mesh.for_each([u, r_cur, rhs, OpL, &vdot_r_cur, &vdot_lr_cur](EuCell cell) {
            cell[r_cur] = OpL(cell, u) - cell[rhs];
            vdot_r_cur += cell[r_cur] * cell[r_cur];
            vdot_lr_cur += OpL(cell, r_cur) * cell[r_cur];
        });

        tau = vdot_r_cur / vdot_lr_cur; //Обновление значения тау

        //Обновление текущего приближения u
      	max_err = mesh.max([u, r_cur, tau](EuCell cell){
            cell[u] -= tau * cell[r_cur];
            return cell[r_cur];
      	});
        counter++;
        if (counter % 1000 == 0) {
            std::cout << "Iteration: " << counter << std::endl;
        }
    } while ((max_err > opts.rtol) && (counter < opts.max_iters));
    if (opts.callback) {
      	opts.callback(mesh, u);
    }
}

int main() {
  	//Многопоточность
	threads::on();

    //Создание сетки
    Rectangle gen(0.0, 1.0, 0.0, sqrt(2), true);  // [x_min, x_max, y_min, y_max]
    gen.set_nx(200);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});
    EuMesh mesh(gen);

    //Переменные для хранения на сетке
    Storable<double> u_cur = mesh.add<double>("u_cur");
    Storable<double> u_orig = mesh.add<double>("u_orig");

    //Файл для записи
    PvdFile pvd("poisson", "output");

    //Переменные для хранения на файле
    pvd.variables.append("u_approx", u_cur);
    pvd.variables.append("u_orig", u_orig);

    //Точное значение
    mesh.for_each(
        [u_orig](EuCell cell) {
			cell[u_orig] = aim_func(cell);
        }
    );

    //Замер скорости
    Stopwatch elapsed;

    //Гиперпарметры
    options_t<> opts;

    //Решатель
    elapsed.resume();
    solver<>(mesh, u_cur, OpL, RHS, opts);
    elapsed.stop();

    std::cout << "Elapsed time: " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)" << std::endl;

    pvd.save(mesh, 1);
    return 0;
}