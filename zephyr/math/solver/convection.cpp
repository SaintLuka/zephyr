#include <zephyr/math/solver/convection.h>

#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/mesh/primitives/bface.h>

namespace zephyr::math {

using mesh::AmrStorage;
using namespace geom;

static const Convection::State U = Convection::datatype();

Convection::State Convection::datatype() {
    return {};
}

Convection::Convection() {
    m_accuracy = 3;
    m_CFL = 0.5;
    m_dt = 1.0e+300;
    m_limiter = "minmod";
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

double Convection::compute_dt(Cell &cell) const {
    double h = cell.incircle_diameter();
    return m_CFL * h / velocity(cell.center()).norm();
}

void Convection::compute_grad(Cell &cell, int stage) {
    double uc = stage < 1 ? cell(U).u1 : cell(U).uh;

    Vector3d grad = Vector3d::Zero();
    for (auto face: cell.faces()) {
        double un = stage < 1 ? face.neib(U).u1 : face.neib(U).uh;
        grad += (0.5 * (uc + un) * face.area()) * face.normal();
    }
    grad /= cell.volume();

    cell(U).ux = grad.x();
    cell(U).uy = grad.y();
    cell(U).uz = grad.z();
}

void Convection::fluxes(Cell &cell, int stage) {
    auto &zc = cell(U);

    double fluxes = 0.0;
    for (auto &face: cell.faces()) {
        const auto &zn =  face.neib(U);

        Vector3d cell_c = cell.center();
        Vector3d face_c = face.center();
        // Пока проблемки с периодичностью
        Vector3d neib_c = 2.0 * face_c - cell_c;

        // Значения в ячейках
        double uc = stage < 1 ? zc.u1 : zc.uh;
        double un = stage < 1 ? zn.u1 : zn.uh;

        // Проекции на грань
        double u_m = uc;
        double u_p = un;

        if (m_accuracy == 2) {
            // Второй порядок
            auto fe = FaceExtra::Simple(
                    m_limiter,
                    uc, zc.ux, zc.uy, zc.uz,
                    un, zn.ux, zn.uy, zc.uz,
                    cell_c, neib_c, face_c);

            u_m = fe.m(uc);
            u_p = fe.p(un);

        }
        else if (m_accuracy > 2) {
            // Третий порядок
            auto fe = FaceExtra::ATvL(
                    uc, zc.ux, zc.uy, 0.0,
                    un, zn.ux, zn.uy, 0.0,
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
        zc.uh = zc.u1 - 0.5 * m_dt * fluxes / cell.volume();
    } else {
        zc.u2 = zc.u1 - m_dt * fluxes / cell.volume();
    }
}

void Convection::update(Mesh &mesh) {
    // Определяем dt
    m_dt = mesh.min(
            [this](Cell &cell) -> double {
                return compute_dt(cell);
            });

    if (m_accuracy < 2) {
        // Схема первого порядка, простой снос значений
        // на промежуточный слой
        mesh.for_each(
                [](Cell &cell) {
                    cell(U).uh = cell(U).u1;
                });
    } else {
        // Схема высокого порядка, считаем производные,
        // выполняем шаг предиктора

        // Считаем производные
        mesh.for_each(
                [this](Cell &cell) {
                    compute_grad(cell, 0);
                });

        // Шаг предиктора
        mesh.for_each(
                [this](Cell &cell) {
                    fluxes(cell, 0);
                });

        // Считаем производные
        mesh.for_each(
                [this](Cell &cell) {
                    compute_grad(cell, 1);
                });
    }

    // Шаг корректора
    mesh.for_each(
            [this](Cell &cell) {
                fluxes(cell, 1);
            });

    // Обновляем слои
    mesh.for_each(
            [this](Cell &cell) {
                cell(U).u1 = cell(U).u2;
                cell(U).uh = 0.0;
                cell(U).u2 = 0.0;
            });
}

void Convection::set_flags(Mesh& mesh) {
    mesh.for_each([](Cell &cell) {
        double min_val = cell(U).u1;
        double max_val = cell(U).u1;

        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            min_val = std::min(min_val, face.neib(U).u1);
            max_val = std::max(max_val, face.neib(U).u1);
        }

        if (max_val - min_val > 0.1) {
            cell.set_flag(+1);
        } else {
            cell.set_flag(-1);
        }
    });
}

Distributor Convection::distributor() const {
    using mesh::Children;

    Distributor distr;

    distr.split = [](AmrStorage::Item &parent, Children &children) {
        for (auto& child: children) {
            Vector3d dr = parent.center - child.center;
            child(U).u1 = parent(U).u1 +
                          parent(U).ux * dr.x() +
                          parent(U).uy * dr.y() +
                          parent(U).uz * dr.z();
        }
    };

    distr.merge = [](Children &children, AmrStorage::Item &parent) {
        double sum = 0.0;
        for (auto &child: children) {
            sum += child(U).u1 * child.volume();
        }
        parent(U).u1 = sum / parent.volume();
    };

    return distr;
}

} // namespace zephyr::math