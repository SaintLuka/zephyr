#pragma once

#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @brief Генератор сетки внутри прямоугольника с отверстием.
class PlaneWithHole final : public Generator {
public:
    using Ptr = std::shared_ptr<PlaneWithHole>;
    using Ref = const std::shared_ptr<PlaneWithHole>&;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::ZOE;
        Boundary right  = Boundary::ZOE;
        Boundary bottom = Boundary::ZOE;
        Boundary top    = Boundary::ZOE;
        Boundary hole   = Boundary::WALL;
    };

    /// @brief Конструктор класса по конфигу
    explicit PlaneWithHole(const Json& config);

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

    /// @brief Число ячеек вдоль верхней границы
    void set_nx(int nx);

    /// @brief Число ячеек вдоль левой границы
    void set_ny(int ny);

    /// @brief Установить флаги граничных условий
    void set_boundaries(Boundaries bounds) const;

    /// @brief Ограничивающий объем
    Box bbox() const override;

    /// @brief Создать сетку
    Grid make() const override;

private:
    void check_params() const;

    // Геометрия
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xc, m_yc, m_r;

    // Базисные вершины
    BaseNode::Ptr v1, v2, v3, v4;
    BaseNode::Ptr v5, v6, v7, v8;
    BaseNode::Ptr v9, v10, v11, v12;
    BaseNode::Ptr v13, v14, v15, v16;
    BaseNode::Ptr v17, v18, v19, v20;

    // Ограничивающие кривые области
    Curve::Ptr circle;
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;

    // Блочная структура
    BlockStructured m_blocks;
};

} // namespace zephyr::geom::generator::collection
