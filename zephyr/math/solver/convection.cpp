#include <zephyr/math/solver/convection.h>

#include <zephyr/math/cfd/face_extra.h>

using namespace zephyr::geom;
using zephyr::mesh::Children;

namespace zephyr::math {

Convection::Convection() {
    m_accuracy = 3;
    m_CFL = 0.5;
    m_dt = 1.0e+300;
    m_limiter = "minmod";
}

void Convection::add_types(EuMesh& mesh) {
    u_curr = mesh.add<double>("u1");
    du_dx = mesh.add<double>("ux");
    du_dy = mesh.add<double>("uy");
    du_dz = mesh.add<double>("uz");
    u_half = mesh.add<double>("uh");
    u_next = mesh.add<double>("u2");
}

double Convection::CFL() const {
    return m_CFL;
}

void Convection::set_CFL(double C) {
    m_CFL = std::max(0.0, std::min(C, 1.0));
}

double Convection::dt() const {
    return m_dt;
}

std::string Convection::limiter() const {
    return m_limiter.name();
}

void Convection::set_limiter(const std::string& limiter) {
    m_limiter = Limiter(limiter);
}

int Convection::accuracy() const {
    return m_accuracy;
}

void Convection::set_accuracy(int acc) {
    m_accuracy = std::max(1, std::min(acc, 3));
}

Vector3d Convection::velocity(const Vector3d& c) const {
    return Vector3d::UnitX();
}

double Convection::compute_dt(EuCell &cell) const {
    double h = cell.incircle_diameter();
    return m_CFL * h / velocity(cell.center()).norm();
}

void Convection::compute_grad(EuCell &cell, int stage) const {
    double uc = stage < 1 ? cell(u_curr) : cell(u_half);

    Vector3d grad = Vector3d::Zero();
    for (auto face: cell.faces()) {
        double un = stage < 1 ? face.neib(u_curr) : face.neib(u_half);
        grad += (0.5 * (uc + un) * face.area()) * face.normal();
    }
    grad /= cell.volume();

    cell(du_dx) = grad.x();
    cell(du_dy) = grad.y();
    cell(du_dz) = grad.z();
}

void Convection::fluxes(EuCell &cell, int stage) {
    double uc = stage < 1 ? cell(u_curr) : cell(u_half);
    double uc_dx = cell(du_dx);
    double uc_dy = cell(du_dy);
    double uc_dz = cell(du_dz);

    double fluxes = 0.0;
    for (auto &face: cell.faces()) {
        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();

        // Пока проблемки с периодичностью
        Vector3d neib_c = 2.0 * face_c - cell_c;

        // Значения в ячейках
        double un = stage < 1 ? face.neib(u_curr) : face.neib(u_half);

        // Проекции на грань
        double u_m = uc;
        double u_p = un;

        double un_dx = face.neib(du_dx);
        double un_dy = face.neib(du_dy);
        double un_dz = face.neib(du_dz);

        if (m_accuracy == 2) {
            // Второй порядок
            auto fe = FaceExtra::Simple(
                    m_limiter,
                    uc, uc_dx, uc_dy, uc_dz,
                    un, un_dx, un_dy, un_dz,
                    cell_c, neib_c, face_c);

            u_m = fe.m(uc);
            u_p = fe.p(un);

        }
        else if (m_accuracy > 2) {
            // Третий порядок
            auto fe = FaceExtra::ATvL(
                    uc, uc_dx, uc_dy, uc_dz,
                    un, un_dx, un_dy, un_dz,
                    cell_c, neib_c, face_c);

            u_m = fe.m(uc);
            u_p = fe.p(un);
        }

        double af = velocity(face.center()).dot(face.normal());
        double a_p = std::max(af, 0.0);
        double a_m = std::min(af, 0.0);

        fluxes += (a_p * u_m + a_m * u_p) * face.area();
    }

    if (stage < 1) {
        cell(u_half) = cell(u_curr) - 0.5 * m_dt * fluxes / cell.volume();
    } else {
        cell(u_next) = cell(u_curr) - m_dt * fluxes / cell.volume();
    }
}

void Convection::update(EuMesh &mesh) {
    // Определяем dt
    m_dt = mesh.min(
            [this](EuCell cell) -> double {
                return compute_dt(cell);
            }, 1.0e300);

    if (m_accuracy < 2) {
        // Схема первого порядка, простой снос значений
        // на промежуточный слой
        mesh.for_each(
                [this](EuCell cell) {
                    cell(u_half) = cell(u_curr);
                });
    } else {
        // Схема высокого порядка, считаем производные,
        // выполняем шаг предиктора

        // Считаем производные
        mesh.for_each(
                [this](EuCell cell) {
                    compute_grad(cell, 0);
                });

        // Шаг предиктора
        mesh.for_each(
                [this](EuCell cell) {
                    fluxes(cell, 0);
                });

        // Считаем производные
        mesh.for_each(
                [this](EuCell cell) {
                    compute_grad(cell, 1);
                });
    }

    // Шаг корректора
    mesh.for_each(
            [this](EuCell cell) {
                fluxes(cell, 1);
            });

    // Обновляем слои
    mesh.swap(u_curr, u_next);
}

void Convection::set_flags(EuMesh& mesh) {
    if (!mesh.adaptive()) {
        return;
    }

    mesh.for_each([this](EuCell cell) {
        double min_val = cell(u_curr);
        double max_val = cell(u_curr);

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            min_val = std::min(min_val, face.neib(u_curr));
            max_val = std::max(max_val, face.neib(u_curr));
        }

        if (max_val - min_val > 0.1) {
            cell.set_flag(+1);
        } else {
            cell.set_flag(-1);
        }
    });
}

Distributor Convection::distributor() const {
    Distributor distr;

    distr.split = [&](EuCell &parent, Children &children) {
        for (auto child: children) {
            Vector3d dr = parent.center() - child.center();
            child(u_curr) = parent(u_curr) +
                            parent(du_dx) * dr.x() +
                            parent(du_dy) * dr.y() +
                            parent(du_dz) * dr.z();
        }
    };

    distr.merge = [&](Children &children, EuCell &parent) {
        double sum = 0.0;
        for (auto child: children) {
            sum += child(u_curr) * child.volume();
        }
        parent(u_curr) = sum / parent.volume();
    };

    return distr;
}

} // namespace zephyr::math