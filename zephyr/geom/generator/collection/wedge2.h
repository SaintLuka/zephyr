#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @class Wedge2. Генератор для создания сетки внутри прямоугольной области
/// с полукруглым вырезом.
class Wedge2 : public BlockStructured {
public:
    using Ptr = std::shared_ptr<Wedge2>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
    };

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc -- Начало клина
    /// @param phi -- Угол наклона клина 
    Wedge2(double xmin, double xmax, double ymin, double ymax,
                  double xc, double phi);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Wedge2::Ptr create(Args&&... args){
        return std::make_shared<Wedge2>(std::forward<Args>(args)...);
    }

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по оси Oy подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_ny(int ny);

    /// @brief Установить флаги граничных условий
    void set_boundaries(Boundaries bounds);

    /// @brief Ограничивающий объем
    Box bbox() const final;

private:
    void check_params() const override;

    void init_blocks();

    // Геометрия
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xc, m_phi;

    // Флаги граничных условий
    Boundaries m_bounds;

    // Куча параметров
    // Базисные вершины для струтурированных блоков
    BaseVertex::Ptr v1, v2, v3, v4, v5, v6;

    // Ограничивающие кривые области
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;
    Curve::Ptr wedge;
};

} // namespace zephyr::geom::generator::collection
