#pragma once

#include <zephyr/configuration.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Простой класс для генерации декартовой сетки или сетки
/// внутри параллелепипеда.
class Cuboid : public Generator {
public:
    using Ptr = std::shared_ptr<Cuboid>;
    using Ref = const std::shared_ptr<Cuboid>&;

    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
        Boundary back   = Boundary::WALL;
        Boundary front  = Boundary::WALL;
    };

    /// @brief Конструктор класса по кофигу
    explicit Cuboid(const Json& config);

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    /// @param ymin, ymax Границы прямоугольника по оси y
    /// @param zmin, zmax Границы прямоугольника по оси z
    Cuboid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Cuboid::Ptr create(Args&&... args){
        return std::make_shared<Cuboid>(std::forward<Args>(args)...);
    }

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
    void set_boundaries(Boundaries bounds);

    /// @brief Количество ячеек сетки
    int size() const final;

    /// @brief Ограничивающий объем
    Box bbox() const final;

    /// @brief Создать сетку общего вида
    Grid make() final;


    // Далее не самые полезные get-функции

    double x_min() const;

    double x_max() const;

    double y_min() const;

    double y_max() const;

    double z_min() const;

    double z_max() const;

    /// @brief Число ячеек по оси X
    int nx() const;

    /// @brief Число ячеек по оси Y
    int ny() const;

    /// @brief Число ячеек по оси Z
    int nz() const;

    /// @brief Есть ли периодичность по оси X?
    bool periodic_along_x() const;

    /// @brief Есть ли периодичность по оси Y?
    bool periodic_along_y() const;

    /// @brief Есть ли периодичность по оси Z?
    bool periodic_along_z() const;

private:
    /// @brief Проверить параметры сетки перед созданием
    void check_params() const final;

    /// @brief Обновить число ячеек
    void compute_size();

    int m_nx, m_ny, m_nz;
    int m_size;
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_zmin, m_zmax;

    Boundaries m_bounds;
};

} // namespace zephyr::geom::generator
