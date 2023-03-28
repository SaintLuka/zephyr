#include "fast.h"

#include <zephyr/math/cfd/cir.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>

using namespace zephyr::phys;
using namespace zephyr::math;

struct _U_ {
    double rho1, rho2;
    Vector3d v1, v2;
    double p1, p2;
    double e1, e2;

    PState get_state1() const {
        return {rho1, v1, p1, e1};
    }

    void set_state2(const PState& z) {
        rho2 = z.density;
        v2 = z.velocity;
        p2 = z.pressure;
        e2 = z.energy;
    }

    void swap() {
        std::swap(rho1, rho2);
        std::swap(v1, v2);
        std::swap(p1, p2);
        std::swap(e1, e2);
    }
};

// Для быстрого доступа по типу
_U_ U;

/// Переменные для сохранения
double get_rho(Storage::iterator cell) { return cell(U).rho1; }
double get_u(Storage::iterator cell)   { return cell(U).v1.x(); }
double get_v(Storage::iterator cell)   { return cell(U).v1.y(); }
double get_V(Storage::iterator cell)   { return cell(U).v1.norm(); }
double get_p(Storage::iterator cell)   { return cell(U).p1 - 1.0_bar; }

/// @brief Уровень звукового давления в дБ
double get_spl(Storage::iterator cell) {
    return 20.0 * std::log(1.0 + std::abs(cell(U).p1 - 1.0_bar) / 20.0e-6_Pa) / std::log(10.0);
}

inline double sqr(double x) { return x*x; }

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"density", get_rho};
    pvd.variables += {"velocity.x", get_u};
    pvd.variables += {"velocity.y", get_v};
    pvd.variables += {"|velocity|", get_V};
    pvd.variables += {"pressure", get_p};
    pvd.variables += {"SPL", get_spl};

    // Тестовая задача
    IdealGas eos("Air");
    double T0 = 20.0_C;
    double P0 = 1.0_bar;
    double R0 = eos.density_pt(P0, T0);

    // Создаем одномерную сетку
    Rectangle rect(0.0, 4.0, -2.0, 4.0);
    rect.set_nx(100);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Создать сетку
    Mesh mesh(U, &rect);

    // Заполняем начальные данные
    for (auto cell: mesh.cells()) {
        cell(U).rho1 = R0;
        cell(U).v1   = Vector3d(0.0, 0.0, 0.0);
        cell(U).p1   = P0;
        cell(U).e1   = eos.energy_rp(cell(U).rho1, cell(U).p1);
    }

    // Число Куранта
    double CFL = 0.5;

    NumFlux::Ptr nf = CIR1::create();

    double time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;

    double max_time = 0.05;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\t\tВремя: " << std::setw(6) << std::setprecision(3) << time << "\n";
            pvd.save(mesh, time);
            next_write += max_time / 200;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto cell: mesh.cells()) {
            // скорость звука
            double c = eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
            for (auto &face: cell.faces()) {
                // Нормальная составляющая скорости
                double vn = cell(U).v1.dot(face.normal());

                // Максимальное по модулю СЗ
                double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

                // Условие КФЛ
                dt = std::min(dt, cell.volume() / face.area() / lambda);
            }
        }
        dt *= CFL;

        // Расчет по схеме CIR
        for (auto cell: mesh.cells()) {
            // Примитивный вектор в ячейке
            PState zc = cell(U).get_state1();

            // Консервативный вектор в ячейке
            QState qc(zc);

            // Переменная для потока
            Flux flux;
            for (auto& face: cell.faces()) {
                // Внешняя нормаль
                auto &normal = face.normal();

                // Примитивный вектор соседа
                PState zn(zc);

                if (face.is_boundary()) {
                    const double d = 0.5;

                    // Граничные условия типа "стенка"
                    double vn = zc.velocity.dot(face.normal());
                    zn.velocity -= 2.0 * vn * face.normal();

                    // Камертон
                    if (face.normal().x() < 0.0 && fabs(face.y()) < d) {
                        // Амплитуда колебаний
                        double A = 0.0;
                        if (fabs(face.y()) < d) {
                            A = sqr(std::cos(M_PI_2 * face.y() / d));
                        }

                        // Непосредственно колебание
                        double T = std::sin(2.0 * M_PI * 440.0 * time);

                        // Затухание
                        double R = std::exp(-100.0 * time);

                        zn.velocity.x() += 1.0e-5 * R * T * A;
                    }
                }
                else {
                    zn = face.neib()(U).get_state1();
                }

                // Значение на грани со стороны ячейки
                PState zm = zc.in_local(normal);

                // Значение на грани со стороны соседа
                PState zp = zn.in_local(normal);

                // Численный поток на грани
                auto loc_flux = nf->flux(zm, zp, eos);
                loc_flux.to_global(normal);

                // Суммируем поток
                flux.vec() += loc_flux.vec() * face.area();
            }

            // Новое значение в ячейке (консервативные переменные)
            QState Qc = qc.vec() - dt * flux.vec() / cell.volume();

            // Новое значение примитивных переменных
            PState Zc(Qc, eos);

            cell(U).set_state2(Zc);
        }

        // Обновляем слои
        for (auto cell: mesh.cells()) {
            cell(U).swap();
        }

        n_step += 1;
        time += dt;
    }

    return 0;
}
