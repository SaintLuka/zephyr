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
   double rtol = 0.0;
   bool use_initial = false;
   std::function<void(EuMesh&, Storable<T>)> callback = nullptr;
};

// Метод сопряженных градиентов
template <typename T = double>
void cg_solver(EuMesh& mesh, Storable<T> u, const operator_t<T>& OpL, const rhs_t<T>& RHS, const options_t<T>& opts) {
  	//Использовать начальное приближение?
	if (!opts.use_initial) {
    	mesh.for_each([u](EuCell cell) {
            cell[u] = 0.0;
    	});
	}

    //Переменные для хранения на сетке. Невязка и правая часть (без невязки никак не обойтись(?))
    Storable<T> r_cur = mesh.add<T>("r_cur");
    Storable<T> p_cur = mesh.add<T>("p_cur");

    // Инициализация начальных значений для r_cur и p_cur
    T max_rhs = mesh.max([r_cur, p_cur, u, OpL, RHS](EuCell cell) {
      	cell[r_cur] = OpL(cell, u) - RHS(cell);
        cell[p_cur] = -cell[r_cur];
        return RHS(cell);
    });

    int counter = 0;
    T tau, max_r;
    T vdot_r_cur = 0, vdot_r_prev = 1, vdot_lp_cur;
    T beta;

    vdot_r_cur = mesh.sum([r_cur](EuCell cell) {
  		return cell[r_cur] * cell[r_cur];
	}, 0.0);

    vdot_lp_cur = mesh.sum([p_cur, OpL](EuCell cell) {
        return OpL(cell, p_cur) * cell[p_cur];
    }, 0.0);

    do {
        tau = vdot_r_cur / vdot_lp_cur;
        vdot_r_prev = vdot_r_cur;
        vdot_r_cur = 0;
        max_r = mesh.max([u, r_cur, p_cur, OpL, tau, &vdot_r_cur](EuCell cell) {
          	cell[u] += tau * cell[p_cur];
          	cell[r_cur] += tau * OpL(cell, p_cur);
            vdot_r_cur += cell[r_cur] * cell[r_cur];
            return cell[r_cur];
        });

        beta = vdot_r_cur / vdot_r_prev;

        mesh.for_each([r_cur, p_cur, beta](EuCell cell) {
          	cell[p_cur] = -cell[r_cur] + beta * cell[p_cur];
        });

        vdot_lp_cur = mesh.sum([p_cur, OpL](EuCell cell) {
          	return OpL(cell, p_cur) * cell[p_cur];
        }, 0.0);

        counter++;
        if (counter % 1000 == 0) {
            std::cout << "Iteration: " << counter << "; Error: " << max_r / max_rhs << std::endl;
        }
    } while ((max_r / max_rhs > opts.rtol) && (counter < opts.max_iters));
    if (opts.callback) {
      	opts.callback(mesh, u);
    }
}

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
    T max_rhs = mesh.max([rhs, RHS](EuCell cell) {
      	cell[rhs] = RHS(cell);
        return cell[rhs];
    });

    int counter = 0;
    T tau, max_r;
    T vdot_r_cur = 0, vdot_lr_cur = 1;

    //Начало цикла
    do {
        tau = vdot_r_cur / vdot_lr_cur; //Обновление значения тау

        vdot_r_cur = 0;
        vdot_lr_cur = 0;

        // Обновление текущего приближения u, подсчет ошибки, обновление значения невязки
    	max_r = mesh.max([u, r_cur, rhs, OpL, tau, &vdot_r_cur, &vdot_lr_cur](EuCell cell){
			cell[u] -= tau * cell[r_cur];

    		// OpL использует соседние ячейки, что может быть некорректно в данном месте
            // Может возникнуть ситуация, когда для текущей ячейки соседние еще не обновили значение u
			cell[r_cur] = OpL(cell, u) - cell[rhs];
			vdot_r_cur += cell[r_cur] * cell[r_cur];
			vdot_lr_cur += OpL(cell, r_cur) * cell[r_cur]; //Аналогичная ситуация только для значения r_cur
			return cell[r_cur];
		});
        // Ниже в закомментирован код для данного метода из теории
//        // Обновление значения u
//      	mesh.for_each([u, r_cur, tau](EuCell cell){
//            cell[u] -= tau * cell[r_cur];
//      	});
//        // Подсчет ошибки и обновление невязки
//    	max_err = mesh.max([u, r_cur, rhs, OpL](EuCell cell){
//			cell[r_cur] = OpL(cell, u) - cell[rhs];
//			return cell[r_cur];
//		});
//        // Подсчет скалярных произведений для вычисления тау
//    	mesh.for_each([r_cur, OpL, &vdot_r_cur, &vdot_lr_cur](EuCell cell){
//			vdot_r_cur += cell[r_cur] * cell[r_cur];
//			vdot_lr_cur += OpL(cell, r_cur) * cell[r_cur];
//		});
        counter++;
        if (counter % 1000 == 0) {
            std::cout << "Iteration: " << counter << "; Error: " << max_r / max_rhs << std::endl;
        }
    } while ((max_r / max_rhs > opts.rtol) && (counter < opts.max_iters));
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