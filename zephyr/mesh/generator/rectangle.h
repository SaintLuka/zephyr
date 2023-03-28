#pragma once

#include <array>

#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/generator/base.h>


namespace zephyr { namespace mesh { namespace generator {


/// @class Rectangle. Простой класс для создания Storage для хранения
/// декартовой сетки или сетки из ячеек Воронного внутри прямоугольной области.
class Rectangle : public Base {
public:

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit Rectangle(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    /// @param ymin, ymax Границы прямоугольника по оси y
    /// @param voronoi Использовать ячейки Вороного
    Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi = false);

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по оси Oy подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_nx(int nx);

    /// @brief Установить желаемое число ячеек сетки по оси Oy
    /// @details Число ячеек по оси Ox подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_ny(int ny);

    /// @brief Установить желаемое число ячеек сетки
    /// @details Число ячеек по осям координат Ox и Oy подбирается так,
    /// чтобы Nx Ny ~ N и aspect ячеек был около единицы
    void set_size(int N);

    /// @brief Установить точные размеры сетки по осям Ox и Oy
    void set_sizes(int nx, int ny);

    /// @brief Установить флаги граничных условий
    void set_boundary_flags(FaceFlag left, FaceFlag right, FaceFlag bottom, FaceFlag top);

    double x_min() const;

    double x_max() const;

    double y_min() const;

    double y_max() const;

    int nx() const;

    int ny() const;

    bool periodic_along_x() const;

    bool periodic_along_y() const;

    Box bbox() const final;

    /// @brief Инициализация переданной сетки
    void initialize(Storage& storage) final;

    /// @brief Инициализация части переданной сетки
    void initialize(Storage& storage, Part part) final;

private:

    /// @brief Проверить параметры сетки перед созданием
    void check_params() const final;

    void compute_size();

    void init_classic(Storage &cells, Part part) const;

    void init_voronoi(Storage &cells, Part part) const;


    int m_nx, m_ny;
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;

    FaceFlag m_left_flag;
    FaceFlag m_right_flag;
    FaceFlag m_bottom_flag;
    FaceFlag m_top_flag;

    bool m_voronoi; ///< use voronoi cells
};

} // generator
} // mesh
} // zephyr
