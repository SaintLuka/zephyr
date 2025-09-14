#include <zephyr/math/solver/free_boundary.h>
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/cfd/gradient.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

using namespace geom;
using namespace smf;

using utils::threads;
using utils::mpi;


FreeBoundary::FreeBoundary(Eos::Ptr eos) : m_eos(eos) {
    m_CFL = 0.5;
    m_dt = NAN;
    m_max_dt = std::numeric_limits<double>::max();
}

FreeBoundary::Parts FreeBoundary::add_types(EuMesh& mesh) {
    part.init = mesh.add<PState>("init");
    part.next = mesh.add<PState>("next");
    return part;
}

void FreeBoundary::set_CFL(double CFL) {
    m_CFL = std::max(0.0, std::min(CFL, 1.0));
}

double FreeBoundary::CFL() const {
    return m_CFL;
}

double FreeBoundary::dt() const {
    return m_dt;
}

void FreeBoundary::set_max_dt(double dt) {
    m_max_dt = dt;
}

namespace{
    PState boundary_value(const PState &zc, const Vector3d &normal, Boundary flag) {
        if (flag != Boundary::WALL) {
            return zc;
        }

        PState zn(zc);
        Vector3d Vn = normal * zc.velocity.dot(normal);
        zn.velocity = zc.velocity - 2.0 * Vn;

        return zn;
    } 
}


void FreeBoundary::update(EuMesh &mesh) {
    // Определяем dt
    compute_dt(mesh);

    mesh.sync(part.init);
    fluxes(mesh);

    // Обновляем слои
    swap(mesh);
}

void FreeBoundary::compute_dt(EuMesh &mesh) {
    double dt = mesh.min([this](EuCell cell) -> double {
        //double c = m_eos->sound_speed_rP(cell(part.density), cell(part.pressure));
        //return cell.incircle_diameter() / (cell(part.velocity).norm() + c);
        double c = m_eos->sound_speed_rP(cell[part.init].density, cell[part.init].pressure);
        return cell.incircle_diameter() / (cell[part.init].velocity.norm() + c);
    });

    dt = std::min(m_CFL * dt, m_max_dt);
    m_dt = mpi::min(dt);
}

void FreeBoundary::fluxes(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        // Примитивный вектор в ячейке
        PState z_c = cell[part.init];

        // Консервативный вектор в ячейке
        QState q_c(z_c);

        // Переменная для потока
        Flux flux;
        for (auto &face: cell.faces()) {
            // Внешняя нормаль
            auto normal = face.normal();

            // Примитивный вектор соседа
            PState z_n;
            if (!face.is_boundary()) {
                z_n = face.neib(part.init);
            } else {
                z_n = boundary_value(z_c, normal, face.flag());
            }

            // Значение на грани со стороны ячейки
            PState zm = z_c.in_local(normal);

            // Значение на грани со стороны соседа
            PState zp = z_n.in_local(normal);

            // Численный поток на грани
            Flux loc_flux = HLLC::calc_flux(zm, zp, *m_eos);
            loc_flux.to_global(normal);

            // Суммируем поток
            // flux.arr() += loc_flux.arr() * face.area(m_axial);
            flux.arr() += loc_flux.arr();
        }

        // Обновляем значение в ячейке (консервативные переменные)
        // q_c.arr() -= (m_dt / cell.volume(m_axial)) * flux.arr();
        q_c.arr() -= m_dt * flux.arr();

        // Новое значение примитивных переменных
        cell(part.next) = PState(q_c, *m_eos);
    });
}

void FreeBoundary::swap(EuMesh &mesh) const {
    mesh.for_each([this](EuCell &cell) {
        cell[part.init] = cell(part.next);
    });
}

} // namespace zephyr::math