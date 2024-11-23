#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @class PlaneWithHole. Генератор для создания сетки внутри прямоугольника с отверстием.
class PlaneWithHole : public BlockStructured {
public:
    using Ptr = std::shared_ptr<PlaneWithHole>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::ZOE;
        Boundary right  = Boundary::ZOE;
        Boundary bottom = Boundary::ZOE;
        Boundary top    = Boundary::ZOE;
        Boundary hole   = Boundary::WALL;
    };

#ifdef ZEPHYR_YAML
    /// @brief Конструктор класса по кофигу
    explicit PlaneWithHole(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc, yc, r Центр и радиус отверстия
    PlaneWithHole(double xmin, double xmax, double ymin, double ymax,
                  double xc, double yc, double r);

    /// @brief Создать указатель на класс
    template <class... Args>
    static PlaneWithHole::Ptr create(Args&&... args){
        return std::make_shared<PlaneWithHole>(std::forward<Args>(args)...);
    }

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по оси Oy подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_nx(int nx);

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
    double m_xc, m_yc, m_r;

    /// @brief Безразмерный параметр, соотношение радиусов "внешней оружности"
    /// и внутренней (отверстия)
    double m_xi;

    // Флаги граничных условий
    Boundaries m_bounds;

    // Куча параметров
    // Базисные вершины для струтурированных блоков

    BaseVertex::Ptr v1, v2, v3, v4;
    BaseVertex::Ptr v5, v6, v7, v8;
    BaseVertex::Ptr v9, v10, v11, v12;
    BaseVertex::Ptr v13, v14, v15, v16;
    BaseVertex::Ptr v17, v18, v19, v20;

    // Ограничивающие кривые области
    Curve::Ptr circle;
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;
};

} // namespace zephyr::geom::generator::collection
