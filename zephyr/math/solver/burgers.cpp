#include <zephyr/math/solver/burgers.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/rotate.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using namespace geom;

using utils::threads;
using utils::mpi;

BurgersSolver::BurgersSolver() {
    m_CFL = 0.5;
    m_limiter = Limiter("MC");
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

BurgersSolver::Parts BurgersSolver::add_types(EuMesh& mesh) {
    part.init = mesh.add<Vector3d>("init");
    part.d_dx = mesh.add<Vector3d>("d_dx");
    part.d_dy = mesh.add<Vector3d>("d_dy");
    part.d_dz = mesh.add<Vector3d>("d_dz");
    part.half = mesh.add<Vector3d>("half");
    part.next = mesh.add<Vector3d>("next");
    return part;
}

void BurgersSolver::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

void BurgersSolver::set_accuracy(int acc) {
    m_acc = std::min(std::max(1, acc), 2);  // 1 или 2
}

void BurgersSolver::set_limiter(const std::string& lim) {
    m_limiter = Limiter(lim);
}

double BurgersSolver::CFL() const {
    return m_CFL;
}

double BurgersSolver::dt() const {
    return m_dt;
}

void BurgersSolver::set_max_dt(double dt) {
    m_max_dt = dt;
}

Vector3d boundary_value(const Vector3d &zc, const Vector3d &normal, Boundary flag) {
    return zc;
}

void BurgersSolver::update(EuMesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    if (m_acc == 1) {
        mesh.sync(part.init);
        fluxes(mesh);
    }
    else {
        mesh.sync(part.init);
        compute_grad(mesh);

        mesh.sync(part.d_dx, part.d_dy, part.d_dz);
        fluxes_stage1(mesh);

        mesh.sync(part.half);
        fluxes_stage2(mesh);
    }

    // Обновляем слои
    swap(mesh);
}

void BurgersSolver::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        return cell.incircle_diameter() / cell[part.init].cwiseAbs().maxCoeff();
    });

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void BurgersSolver::compute_grad(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        auto grad = gradient::LSM<Vector3d>(cell, part.init, boundary_value);
        grad = gradient::limiting<Vector3d>(cell, m_limiter, grad, part.init, boundary_value);

        cell[part.d_dx] = grad.x;
        cell[part.d_dy] = grad.y;
        cell[part.d_dz] = grad.z;
    });
}

Vector3d diff_flux(const Vector3d& q) {
    return {0.5 * q.x() * q.x(), q.x() * q.y(), q.x() * q.z()};
}

Vector3d quasi_upwind(const Vector3d& zL, const Vector3d& zR) {
    double s = 0.5 * (zL.x() + zR.x());
    return s >= 0.0 ? diff_flux(zL) : diff_flux(zR);
}

Vector3d cir_like(const Vector3d& zL, const Vector3d& zR) {
    Vector3d Fl = diff_flux(zL);
    Vector3d Fr = diff_flux(zR);
    double s = 0.5 * (zL.x() + zR.x());
    return 0.5 * (Fl + Fr) + 0.5 * std::abs(s) * (zL - zR);
}

Vector3d calc_flux(const Vector3d& zL, const Vector3d& zR) {
    return quasi_upwind(zL, zR);
    //return cir_like(zL, zR);
}

void BurgersSolver::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Консервативный вектор в ячейке
        Vector3d z_c = cell[part.init];

        // Переменная для потока
        Vector3d flux = Vector3d::Zero();
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            Vector3d z_n;
            if (!face.is_boundary()) {
                z_n = face.neib(part.init);
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            Vector3d zm = z_c; Rotate::to_local(zm, normal);

            // Значение на грани со стороны соседа
            Vector3d zp = z_n; Rotate::to_local(zp, normal);

            // Численный поток на грани
            Vector3d loc_flux = calc_flux(zm, zp);
            Rotate::to_global(loc_flux, normal);

            // Суммируем поток
            flux += loc_flux * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        Vector3d S = {
            z_c.x() * (cell[part.d_dy].y() + cell[part.d_dz].z()),
            z_c.y() * (cell[part.d_dx].x() + cell[part.d_dz].z()),
            z_c.z() * (cell[part.d_dx].z() + cell[part.d_dy].y())
        };
        z_c -= (m_dt / cell.volume()) * flux;// + m_dt * S;

        // Новое значение примитивных переменных
        cell[part.next] = z_c;
    });
}

void BurgersSolver::fluxes_stage1(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор в ячейке
        Vector3d z_c = cell[part.init];

        // Переменная для потока
        Vector3d flux = Vector3d::Zero();
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа
            Vector3d z_n;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n = neib[part.init];
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n = boundary_value(z_c, normal, face.flag());
            }

            auto face_extra = FaceExtra::Direct(
                    z_c, cell(part.d_dx), cell(part.d_dy), cell(part.d_dz),
                    z_n, neib(part.d_dx), neib(part.d_dy), neib(part.d_dz),
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            Vector3d zm = face_extra.m(z_c);

            // Переводим в локальную систему координат
            Rotate::to_local(zm, normal);

            // Численный поток на грани
            Vector3d loc_flux = diff_flux(zm);
            Rotate::to_global(loc_flux, normal);

            // Суммируем поток
            flux += loc_flux * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        z_c -= (0.5 * m_dt / cell.volume()) * flux;

        // Значение примитивных переменных на полушаге
        cell[part.half] = z_c;
    });
}

void BurgersSolver::fluxes_stage2(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Центр ячейки
        Vector3d cell_c = cell.center();

        // Примитивный вектор на полуслое
        Vector3d z_c = cell[part.init];

        // Примитивный вектор на полуслое
        Vector3d z_ch = cell[part.half];

        // Переменная для потока (суммирование по промежуточным)
        Vector3d flux = Vector3d::Zero();
        for (auto &face: cell.faces()) {
            // Внешняя нормаль и центр грани
            auto  normal = face.normal();
            auto &face_c = face.center();

            // Возвращает саму ячейку, если соседа не существует
            auto neib = face.neib();

            // Примитивный вектор соседа (на предыдущем и на полушаге)
            Vector3d z_n, z_nh;
            Vector3d neib_c;
            if (!face.is_boundary()) {
                neib_c = neib.center();
                z_n  = neib[part.init];
                z_nh = neib(part.half);
            }
            else {
                neib_c = face.symm_point(cell_c);
                z_n  = boundary_value(z_c,  normal, face.flag());
                z_nh = boundary_value(z_ch, normal, face.flag());
            }

            // Параметры интерполяции с предыдущего (!) слоя
            auto face_extra = FaceExtra::Direct(
                    z_c, cell(part.d_dx), cell(part.d_dy), cell(part.d_dz),
                    z_n, neib(part.d_dx), neib(part.d_dy), neib(part.d_dz),
                    cell_c, neib_c, face_c);

            // Интерполяция на грань со стороны ячейки
            Vector3d zm = face_extra.m(z_ch);

            // Интерполяция на грань со стороны соседа
            Vector3d zp;
            if (!face.is_boundary()) {
                zp = face_extra.p(z_nh);
            }
            else {
                zp = boundary_value(zm, normal, face.flag());
            }

            // Переводим в локальную систему координат
            Rotate::to_local(zm, normal);
            Rotate::to_local(zp, normal);

            // Численный поток на грани
            Vector3d loc_flux = calc_flux(zm, zp);
            Rotate::to_global(loc_flux, normal);

            // Суммируем поток
            flux += loc_flux * face.area();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        z_c -= (m_dt / cell.volume()) * flux;

        // Значение примитивных переменных на новом слое
        cell[part.next] = z_c;
    });
}

void BurgersSolver::swap(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        cell[part.init] = cell[part.next];
    });
}

Distributor BurgersSolver::distributor(const std::string& type) const {
    if (type != "const" && type != "slope") {
        throw std::runtime_error("BurgersSolver error: unknown m_distributor type '" + type + "'");
    }
    
    using mesh::Children;
    
    Distributor distr;

    // Консервативное суммирование
    distr.merge = [this](const Children &children, EuCell &parent) {
        Vector3d q_p = Vector3d::Zero();
        for (auto child: children) {
            Vector3d q_ch = child[part.init];
            q_p += q_ch * child.volume();
        }
        parent[part.init] = q_p / parent.volume();
    };

    // Снос копированием
    auto split_const = [this](const EuCell &parent, Children &children) {
        Vector3d z_p = parent[part.init];
        for (auto child: children) {
            child[part.init] = z_p;
        }
    };
    
    // Снос по градиентам
    auto split_slope = [this](const EuCell &parent, Children &children) {
        Vector3d z_p  = parent[part.init];
        Vector3d d_dx = parent[part.d_dx];
        Vector3d d_dy = parent[part.d_dy];
        Vector3d d_dz = parent[part.d_dz];

        for (auto child: children) {
            Vector3d dr = child.center() - parent.center();
            child[part.init] = z_p + d_dx * dr.x() + d_dy * dr.y() + d_dz * dr.z();
        }
    };
    
    if (type == "const") {
        // Снос копированием
        distr.split = split_const;
    }
    else {
        // Снос по градиентам
        distr.split = split_slope;
    }

    return distr;
}

void BurgersSolver::set_flags(EuMesh &mesh) const {
    if (!mesh.adaptive()) {
        return;
    }

    compute_grad(mesh);

    // Пороги (относительные) на разбиение
    const double xi_norm = 0.05;

    for (auto cell: mesh) {
        //cell.set_flag(1); continue;
        cell.set_flag(-1);

        double norm = cell[part.init].norm();

        double norm_split = xi_norm * std::abs(norm);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }

            double norm_n = face.neib(part.init).norm();

            // Большой перепад плотностей или давлений
            if (std::abs(norm_n - norm) > norm_split) {
                cell.set_flag(1);
                break;
            }

            // Пороги минимум в два раза меньше
            if (std::abs(norm_n - norm) > 0.4 * norm_split) {
                cell.set_flag(0);
            }
        }
    }
}

} // namespace zephyr::math