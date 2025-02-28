#pragma once

#include <zephyr/configuration.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @class Rectangle. Простой класс для генерации декартовой сетки
/// или сетки из ячеек Воронного внутри прямоугольной области.
class Rectangle : public Generator {
public:
    using Ptr = std::shared_ptr<Rectangle>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
    };

#ifdef ZEPHYR_YAML
    /// @brief Конструктор класса по кофигу
    explicit Rectangle(YAML::Node config);
#endif

    /// @brief Единичный квадрат из одной ячейки
    Rectangle();

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    /// @param ymin, ymax Границы прямоугольника по оси y
    /// @param voronoi Использовать ячейки Вороного
    Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi = false);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Rectangle::Ptr create(Args&&... args){
        return std::make_shared<Rectangle>(std::forward<Args>(args)...);
    }

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
    void set_boundaries(Boundaries bounds);

    /// @brief Количество ячеек сетки
    int size() const final;

    /// @brief Ограничивающий объем
    Box bbox() const final;

    /// @brief Создать сетку общего вида
    Grid make() final;


    // Далее не самые полезные get-функции

    /// @brief Левая граница
    double x_min() const;

    /// @brief Правая граница
    double x_max() const;

    /// @brief Нижняя граница
    double y_min() const;

    /// @brief Правая граница
    double y_max() const;

    /// @brief Число ячеек по оси X
    int nx() const;

    /// @brief Число ячеек по оси Y
    int ny() const;

    /// @brief Есть ли периодичность по оси X?
    bool periodic_along_x() const;

    /// @brief Есть ли периодичность по оси Y?
    bool periodic_along_y() const;

private:
    /// @brief Проверить параметры сетки перед созданием
    void check_params() const final;

    /// @brief Обновить число ячеек
    void compute_size();

    /// @brief Создать классическую декартову сетку
    Grid create_classic() const;

    /// @brief Создать сетку из шестиугольников
    Grid create_voronoi() const;


    int m_nx, m_ny;         ///< Число ячеек по осям
    int m_size;             ///< Суммарное число ячеек
    double m_xmin, m_xmax;  ///< Границы области по оси X
    double m_ymin, m_ymax;  ///< Границы области по оси Y
    Boundaries m_bounds;    ///< Граничные условия
    bool m_voronoi;         ///< Использовать ячейки Вороного
};

} // namespace zephyr::geom::generator
