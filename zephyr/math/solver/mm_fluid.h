#pragma once

#include <boost/format.hpp>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/gradient.h>

namespace zephyr::math {

using geom::Vector3d;
using mesh::EuCell;
using mesh::EuMesh;
using mesh::Storable;
using mesh::Direction;
using mesh::Distributor;
using gradient::GetState;
using gradient::GetBoundary;

using phys::Fractions;
using phys::VectorSet;

using namespace zephyr::math::mmf;

/// @brief Тип расщепления по направлениям
enum class DirSplit {
    NONE,    ///< Без расщепления
    SIMPLE,  ///< Схема первого порядка
    STRANG,  ///< Схема второго порядка (G. Strang 1968)
};

/// @brief Метод выбора alpha_sigma в расчетах методом CRP
enum class CrpMode {
    NONE,       ///< Не использовать CRP
    V3, V5, VS, ///< Эвристические шаблоны
    PLIC,       ///< Линейная реконструкция границы
    MUSCL,      ///< По градиенту объемных долей
};

/// @brief Multi-Material Fluid. Решатель для многоматериальной гидро- и
/// газодинамики с равновесием по давлению и температуре.
class MmFluid {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct Parts {
        Storable<PState> init;  ///< Состояние на начало шага
        Storable<PState> half;  ///< Состояние на полушаге
        Storable<PState> next;  ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        Storable<PState> d_dx, d_dy, d_dz;

        Storable<VectorSet> n;      ///< Нормаль к реконструкции границы
        Storable<VectorSet> p;      ///< Точка реконструкции границы
        Storable<VectorSet> grad_a; ///< Градиент объемных долей (для CrpMode::MUSCL)
    };

    /*
        bool is_bad1() const { return get_state().is_bad(); }

        bool is_bad2() const { return half.is_bad(); }

        bool is_bad3() const { return next.is_bad(); }

        bool is_bad() const { return is_bad1() || is_bad3(); }

        /// @brief Собрать вектор состояния на предыдущем шаге
        PState get_state() const {
            return PState(density, velocity, pressure, energy,
                          temperature, mass_frac, densities);
        }

        /// @brief Установить вектор состояния на предыдущем шаге
        void set_state(const PState &z) {
            density     = z.density;
            velocity    = z.velocity;
            pressure    = z.pressure;
            energy      = z.energy;
            temperature = z.temperature;

            mass_frac   = z.mass_frac;
            densities   = z.densities;
        }

        /// @brief Массив объемных долей
        double vol_frac(int idx) const;

        /// @brief Массив объемных долей
        Fractions vol_fracs() const;
    };

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const State &state);
    */

    Parts part;



    /// @brief Конструктор класса
    explicit MmFluid(const phys::MixturePT &eos);

    /// @brief Деструктор
    ~MmFluid() = default;

    /// @brief Добавить типы на сетку
    Parts add_types(EuMesh& mesh);

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Установить режим CRP
    void set_crp_mode(CrpMode mode);

    /// @brief Задать ограничитель
    void set_limiter(const std::string& limiter);

    /// @brief Установить метод расщепления по направлениям
    void set_splitting(DirSplit splitting);

    /// @brief Задать ускорение свободного падения
    void set_gravity(double g);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Установить максимальный шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Основная функция решателя, сделать шаг
    void update(EuMesh &mesh);

    /// @brief Сделать отсечение, построить поверхность
    EuMesh body(EuMesh& mesh, int idx) const;

    /// @brief Установить флаги адаптации
    void set_flags(EuMesh &mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

public:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(EuMesh &mesh);

    /// @brief Проинтегрировать на шаг dt, вдоль направления dir
    void integrate(EuMesh &mesh, double dt, Direction dir = Direction::ANY);

    /// @brief Лимитированный градиент вектора состояния
    void compute_grad(EuMesh &mesh, Storable<PState> U);

    /// @brief Лимитированный градиент объемных долей
    void fractions_grad(EuMesh &mesh, Storable<PState> U);

    /// @brief Подсеточная линейная реконструкция интерфейса
    void interface_recovery(EuMesh &mesh);

    /// @brief Расчёт потоков с первым порядком
    void fluxes(EuMesh &mesh, double dt, Direction dir = Direction::ANY);

    /// @brief Стадия предиктора при расчете со вторым порядком
    void fluxes_stage1(EuMesh &mesh, double dt, Direction dir = Direction::ANY);

    /// @brief Стадия корректора при расчете со вторым порядком
    void fluxes_stage2(EuMesh &mesh, double dt, Direction dir = Direction::ANY);

    /// @brief Обмен слоев
    void swap(EuMesh &mesh);

    Flux calc_crp_flux(const PState& zL, const PState& zR, double hL, double hR, int iA, double a_sig, double dt);

    Flux calc_flux(mesh::EuCell& cell, mesh::EuFace& face,
                   const PState& z_L, const PState& z_R,
                   double h_L, double h_R, double dt);

    //double alpha_sigma(Cell& cell_L, Cell& cell_R, mesh::Face& face, int idx) const;

protected:
    MixturePT mixture;  ///< Список материалов (смесь)
    NumFlux::Ptr m_nf;  ///< Метод расчёта потока
    CrpMode m_crp_mode; ///< Композитная задача Римана?
    DirSplit m_split;   ///< Расщепление по направлениям
    int m_acc = 1;      ///< Порядок точности
    Limiter m_limiter;  ///< Ограничитель для методов второго порядка
    double m_CFL;       ///< Число Куранта
    double m_g = 0.0;   ///< Ускорение свободного падения, направлено против oy
    double m_dt;        ///< Шаг интегрирования
    double m_max_dt;    ///< Максимальный шаг интегрирования
};

/*
std::ostream &operator<<(std::ostream &os, const MmFluid::State &state) {
    os << boost::format(
            "Main state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() % state.pressure %
          state.temperature % state.energy % state.mass_frac;
    os << boost::format(
            "Half state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
          state.half.density % state.half.velocity.x() % state.half.velocity.y() % state.half.velocity.z() %
          state.half.pressure % state.half.temperature % state.half.energy % state.half.mass_frac;
    os << boost::format(
            "Next state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
          state.next.density % state.next.velocity.x() % state.next.velocity.y() % state.next.velocity.z() %
          state.next.pressure % state.next.temperature % state.next.energy % state.next.mass_frac;
    return os;
}
*/

} // namespace zephyr