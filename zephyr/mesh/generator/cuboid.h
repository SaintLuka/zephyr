#pragma once

#include <array>

#include <zephyr/mesh/storage.h>

#ifdef ZEPHYR_ENABLE_MPI
#include <zephyr/network/mpi/network.h>
#endif

#include <zephyr/mesh/generator/generator.h>

namespace zephyr { namespace mesh { namespace generator {

#ifdef ZEPHYR_ENABLE_MPI
using ::zephyr::network::mpi::Network;
#endif

/// @class Cuboid. Простой класс для создания Storage для хранения
/// декартовой сетки или сетки из ячеек Воронного внутри прямоугольной области.
class Cuboid : public Generator {
public:

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit Cuboid(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    /// @param ymin, ymax Границы прямоугольника по оси y
    /// @param zmin, zmax Границы прямоугольника по оси z
    Cuboid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по осям Oy и Oz подбирается так, чтобы aspect
    /// ячеек был около единицы
    void set_nx(int nx);

    /// @brief Установить желаемое число ячеек сетки по оси Oy
    /// @details Число ячеек по осям Ox и Oz подбирается так, чтобы aspect
    /// ячеек был около единицы
    void set_ny(int ny);

    /// @brief Установить желаемое число ячеек сетки по оси Oz
    /// @details Число ячеек по осям Ox и Oy подбирается так, чтобы aspect
    /// ячеек был около единицы
    void set_nz(int nz);

    /// @brief Установить желаемое число ячеек сетки
    /// @details Число ячеек по осям координат Ox, Oy и Oz подбирается так,
    /// чтобы Nx Ny Nz ~ N и aspect ячеек был около единицы
    void set_size(int N);

    /// @brief Установить точные размеры сетки по осям Ox, Oy, Oz
    void set_sizes(int nx, int ny, int nz);

    /// @brief Установить флаги граничных условий
    void set_boundary_flags(FaceFlag left, FaceFlag right,
                            FaceFlag bottom, FaceFlag top,
                            FaceFlag back, FaceFlag front);

    void initialize(Storage &cells) final;

    double xmin() const;

    double xmax() const;

    double ymin() const;

    double ymax() const;

    double zmin() const;

    double zmax() const;

    int nx() const;

    int ny() const;

    int nz() const;

    bool periodic_along_x() const;

    bool periodic_along_y() const;

    bool periodic_along_z() const;

    Box bbox() const final;

private:

    void compute_size();

    void check_params() const final;

    int m_nx, m_ny, m_nz;
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    FaceFlag m_left_flag, m_right_flag;
    FaceFlag m_bottom_flag, m_top_flag;
    FaceFlag m_back_flag, m_front_flag;
};

} // generator
} // mesh
} // zephyr
