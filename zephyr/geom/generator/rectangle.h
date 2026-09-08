#pragma once

#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Простой класс для генерации декартовой сетки
/// или сетки из ячеек Вороного внутри прямоугольной области.
/// По умолчанию генерируется адаптивная сетка!
class Rectangle final : public Generator {
public:
    using Ptr = std::shared_ptr<Rectangle>;
    using Ref = const std::shared_ptr<Rectangle>&;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
    };

    /// @brief Единичный квадрат из одной ячейки
    Rectangle();

    /// @brief Конструктор класса по конфигу
    explicit Rectangle(const Json& config);

    /// @brief Конструктор класса
    /// @param x_min, x_max Границы прямоугольника по оси x
    /// @param y_min, y_max Границы прямоугольника по оси y
    /// @param voronoi Использовать ячейки Вороного
    Rectangle(double x_min, double x_max, double y_min, double y_max, bool voronoi = false);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Rectangle::Ptr create(Args&&... args){
        return std::make_shared<Rectangle>(std::forward<Args>(args)...);
    }

    /// @brief Структурированная сетка?
    bool structured() const override { return !voronoi_; }

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

    /// @brief Использовать осевую симметрию
    void set_axial(bool axial) override;

    /// @brief Использовать адаптацию
    void set_adaptive(bool adaptive) override;

    /// @brief Использовать линейную адаптацию
    void set_linear(bool) override { linear_ = true; }

    /// @brief Ограничивающий объем
    Box bbox() const override;

    /// @brief Создать сетку общего вида
    Grid make() const override;

    /// @brief Может инициализировать хранилище (для адаптивной декартовой сетки)
    bool can_make_cells() const override { return adaptive_ && !voronoi_; }

    /// @brief Инициализация SoA-хранилища сетки
    mesh::AmrCells make_cells(bool unique_nodes) const override;

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

    /// @brief Граничные условия
    Boundaries bounds() const;

    /// @brief Есть ли периодичность по оси X?
    bool periodic_along_x() const;

    /// @brief Есть ли периодичность по оси Y?
    bool periodic_along_y() const;

private:
    /// @brief Проверить параметры сетки перед созданием
    void check_params() const;

    /// @brief Обновить число ячеек
    void compute_size();

    /// @brief Создать классическую декартову сетку
    Grid create_classic() const;

    /// @brief Создать классическую декартову сетку
    Grid create_classic_amr() const;

    /// @brief Создать сетку из шестиугольников
    Grid create_voronoi() const;

    /// @brief Создать классическую декартову сетку
    void initialize_classic(mesh::AmrCells& cells) const;

    /// @brief Создать сетку из шестиугольников
    void initialize_voronoi(mesh::AmrCells& cells);

    int nx_{0}, ny_{0};     ///< Число ячеек по осям
    int size_{0};           ///< Суммарное число ячеек
    double x_min_, x_max_;  ///< Границы области по оси X
    double y_min_, y_max_;  ///< Границы области по оси Y
    Boundaries bounds_;    ///< Граничные условия
    bool voronoi_ = false; ///< Использовать ячейки Вороного
};

} // namespace zephyr::geom::generator
