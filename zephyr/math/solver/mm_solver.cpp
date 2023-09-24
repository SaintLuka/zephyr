#include <zephyr/math/solver/mm_solver.h>
#include <zephyr/math/cfd/face_extra.h>
#include "zephyr/math/cfd/models.h"

namespace zephyr {
    namespace math {
        using mesh::Storage;
        using namespace geom;
        using namespace smf;

        static const MmSolver::State U = MmSolver::datatype();

        [[nodiscard]] MmSolver::State MmSolver::datatype(){
            return {};
        }

        MmSolver::MmSolver(const phys::Eos& eos, Fluxes flux=Fluxes::HLLC2) :m_eos(eos){
            m_nf = std::move(flux_from_enum(flux));
            m_CFL = 0.9;
            m_dt = std::numeric_limits<double>::max();
        }

        [[nodiscard]] double MmSolver::CFL() const {
            return m_CFL;
        }

        void MmSolver::set_CFL(double CFL) {
            m_CFL = std::max(0.0, std::min(CFL, 1.0));
        }

        [[nodiscard]] double MmSolver::dt() const {
            return m_dt;
        }

        void MmSolver::init_cells(Mesh &mesh, const phys::ClassicTest &test) {
            // Заполняем начальные данные
            for (auto cell: mesh.cells()) {
                cell(U).rho1 = test.density(cell.center());
                cell(U).v1 = test.velocity(cell.center());
                cell(U).p1 = test.pressure(cell.center());
                cell(U).e1 = m_eos.energy_rp(cell(U).rho1, cell(U).p1);
            }
        }

        double MmSolver::compute_dt(Mesh &mesh) {
            m_dt = std::numeric_limits<double>::max();
            for (auto cell: mesh.cells()) {
                // скорость звука
                double c = m_eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
                for (auto &face: cell.faces()) {
                    // Нормальная составляющая скорости
                    double vn = cell(U).v1.dot(face.normal());

                    // Максимальное по модулю СЗ
                    double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

                    // Условие КФЛ
                    m_dt = std::min(m_dt, cell.volume() / face.area() / lambda);
                }
            }
            m_dt *= m_CFL;

            return m_dt;
        }

        void MmSolver::fluxes(Mesh &mesh) {
            // Расчет по некоторой схеме
            for (auto cell: mesh.cells()) {
                // Примитивный вектор в ячейке
                PState zc(cell(U).rho1, cell(U).v1, cell(U).p1, cell(U).e1);

                // Консервативный вектор в ячейке
                QState qc(zc);

                // Переменная для потока
                Flux flux;
                for (auto &face: cell.faces()) {
                    // Внешняя нормаль
                    auto &normal = face.normal();

                    // Примитивный вектор соседа
                    PState zn(zc);

                    if (!face.is_boundary()) {
                        auto neib = face.neib();
                        zn.density = neib(U).rho1;
                        zn.velocity = neib(U).v1;
                        zn.pressure = neib(U).p1;
                        zn.energy = neib(U).e1;
                    }

                    // Значение на грани со стороны ячейки
                    PState zm = zc.in_local(normal);

                    // Значение на грани со стороны соседа
                    PState zp = zn.in_local(normal);

                    PState f = zm.vec() - zp.vec();

                    // Численный поток на грани
                    auto loc_flux = m_nf->flux(zm, zp, m_eos);
                    loc_flux.to_global(normal);

                    // Суммируем поток
                    flux.vec() += loc_flux.vec() * face.area();
                }

                // Новое значение в ячейке (консервативные переменные)
                QState Qc = qc.vec() - m_dt * flux.vec() / cell.volume();

                // Новое значение примитивных переменных
                PState Zc(Qc, m_eos);

                cell(U).rho2 = Zc.density;
                cell(U).v2 = Zc.velocity;
                cell(U).p2 = Zc.pressure;
                cell(U).e2 = Zc.energy;
            }
        }

        void MmSolver::update(Mesh &mesh) {
            for (auto cell: mesh.cells()) {
                std::swap(cell(U).rho1, cell(U).rho2);
                std::swap(cell(U).v1, cell(U).v2);
                std::swap(cell(U).p1, cell(U).p2);
                std::swap(cell(U).e1, cell(U).e2);
            }
            m_time += m_dt;
            m_step += 1;
        }

        void MmSolver::set_num_flux(Fluxes flux) {
            m_nf = std::move(flux_from_enum(flux));
        }

        [[nodiscard]] double MmSolver::get_time() const {
            return m_time;
        }

        [[nodiscard]] size_t MmSolver::getStep() const {
            return m_step;
        }

        [[nodiscard]] std::string MmSolver::get_flux_name() const{
            return m_nf->get_name();
        }
    }
}