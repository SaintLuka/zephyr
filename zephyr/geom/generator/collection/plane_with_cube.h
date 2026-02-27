#pragma once

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @brief Генератор для создания сетки внутри прямоугольника
/// с квадратным отверстием со стороной r.
class PlaneWithCube final : public Generator {
public:
    using Ptr = std::shared_ptr<PlaneWithCube>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::ZOE;
        Boundary right  = Boundary::ZOE;
        Boundary bottom = Boundary::ZOE;
        Boundary top    = Boundary::ZOE;
        Boundary hole   = Boundary::WALL;
    };

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc, yc, r Центр и радиус отверстия
    PlaneWithCube(double xmin, double xmax, double ymin, double ymax,
                  double xc, double yc, double r);

    /// @brief Создать указатель на класс
    template <class... Args>
    static PlaneWithCube::Ptr create(Args&&... args){
        return std::make_shared<PlaneWithCube>(std::forward<Args>(args)...);
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
    Grid make() const override { return m_blocks.make(); }

private:
    void check_params() const;

    // Геометрия
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xc, m_yc, m_r;

    Plane::Ptr left, right, bottom, top;
    Plane::Ptr side1, side2, side3, side4;

    // Базисные вершины
    BaseNode::Ptr v01, v02, v03, v04;
    BaseNode::Ptr v05, v06, v07, v08;
    BaseNode::Ptr v09, v10, v11, v12;
    BaseNode::Ptr v13, v14, v15, v16;
    BaseNode::Ptr v17, v18, v19, v20;
    BaseNode::Ptr v21, v22, v23, v24;
    BaseNode::Ptr v25, v26, v27, v28;
    BaseNode::Ptr v29, v30, v31, v32;

    // Блочная структура
    BlockStructured m_blocks;
};

} // namespace zephyr::geom::generator::collection
