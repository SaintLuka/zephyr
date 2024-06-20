#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @class SemicircleCutout. Генератор для создания сетки внутри прямоугольной области
/// с полукруглым вырезом.
class SemicircleCutout : public BlockStructured {
public:
    using Ptr = std::shared_ptr<SemicircleCutout>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
    };

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc -- Центр выреза
    /// @param r -- Радиус выреза
    SemicircleCutout(double xmin, double xmax, double ymin, double ymax,
                  double xc, double r, Boundaries bounds);

    /// @brief Создать указатель на класс
    template <class... Args>
    static SemicircleCutout::Ptr create(Args&&... args){
        return std::make_shared<SemicircleCutout>(std::forward<Args>(args)...);
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
    double m_xc, m_r;

    // Флаги граничных условий
    Boundaries m_bounds;

    // Куча параметров
    // Базисные вершины для струтурированных блоков
    BaseVertex::Ptr v1, v2, v3, v4;
    BaseVertex::Ptr v5, v6, v7, v8;
    BaseVertex::Ptr v9, v10, v11, v12;

    // Ограничивающие кривые области
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;
    Curve::Ptr circle;
};

} // namespace zephyr::geom::generator::collection
