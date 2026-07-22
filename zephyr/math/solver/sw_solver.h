#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::mesh::Storable;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;
using zephyr::phys::Eos;

using namespace swe;

/// @brief Описание дна
class IBed {
public:
    using Ptr = std::shared_ptr<IBed>;
    using Ref = const std::shared_ptr<IBed>&;

    virtual ~IBed() = default;

    /// @brief Высота дна и градиент в зависимости от координаты
    virtual std::tuple<double, Vector2d> get(const Vector3d& v) const {
        throw std::runtime_error("IBed::get() is not implemented");
    }
};

class ConstBed : public IBed {
    std::tuple<double, Vector2d> res;
public:
    explicit ConstBed(double level) : res({level, Vector2d::Zero()}) {}

    /// @brief Высота дна в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        return res;
    }
};

/// @class SwSolver sw_solver.h
/// @brief Shallow-Water Solver. Решатель для мелкой воды.
class SwSolver {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct Parts {
        Storable<PState> init;  ///< Состояние на основном слое
        Storable<PState> half;  ///< Состояние на полушаге
        Storable<PState> next;  ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        Storable<PState> d_dx, d_dy;
    };

    Parts part;

    /// @brief Конструктор класса
    SwSolver(double bed);

    /// @brief Деструктор
    ~SwSolver() = default;

    /// @brief Добавить типы на сетку
    Parts add_types(EuMesh& mesh);

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Установить ограничитель градиента
    void set_limiter(const std::string& limiter);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Установить шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Выполнить шаг интегрирования по времени
    void update(EuMesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(EuMesh& mesh) const;

    /// @brief Распределитель данных при адаптации
    /// @param type Тип "const" или "slope" переноса при разбиении
    Distributor distributor(const std::string& type = "slope") const;


    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(EuMesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(EuMesh &mesh) const;

    /// @brief Обновление ячеек
    void swap(EuMesh &mesh) const;

    /// @brief Вычислить производные
    void compute_grad(EuMesh &mesh) const;

    /// @brief Вычислить потоки на стадии предиктора
    void fluxes_stage1(EuMesh &mesh) const;

    /// @brief Вычислить потоки на стадии корректора
    void fluxes_stage2(EuMesh &mesh) const;

protected:
    IBed::Ptr m_bed;         ///< Положение дна
    NumFlux::Ptr m_nf;       ///< Метод расчёта потока
    int m_acc = 1;           ///< Порядок точности
    Limiter m_limiter;       ///< Ограничитель градиента
    double m_CFL = 0.5;      ///< Число Куранта
    double m_dt;             ///< Шаг интегрирования
    double m_max_dt=1.e300;  ///< Максимальный шаг интегрирования
};

} // namespace zephyr::math