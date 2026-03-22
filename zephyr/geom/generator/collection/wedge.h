#pragma once

#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @brief Генератор для создания сетки внутри прямоугольной области
/// с отрезанным углом.
class Wedge final : public Generator {
public:
    using Ptr = std::shared_ptr<Wedge>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary left   = Boundary::WALL;
        Boundary right  = Boundary::WALL;
        Boundary bottom = Boundary::WALL;
        Boundary top    = Boundary::WALL;
    };

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xw Положение клина
    /// @param phi Угол наклона
    Wedge(double xmin, double xmax, double ymin, double ymax, double xw, double phi);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Wedge::Ptr create(Args&&... args){
        return std::make_shared<Wedge>(std::forward<Args>(args)...);
    }

    /// @brief Число ячеек вдоль верхней границы
    void set_nx(int nx);

    /// @brief Число ячеек вдоль левой границы
    void set_ny(int ny);

    /// @brief Установить флаги граничных условий
    void set_boundaries(Boundaries bounds) const;

    /// @brief Использовать осевую симметрию
    void set_axial(bool axial) override;

    /// @brief Использовать адаптацию
    void set_adaptive(bool adaptive) override;

    /// @brief Использовать линейную адаптацию
    void set_linear(bool linear) override;

    /// @brief Ограничивающий объем
    Box bbox() const override;

    /// @brief Создать сетку
    Grid make() const override;

private:
    void check_params() const;

    // Геометрия
    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xw, m_phi;

    // Базисные вершины для структурированных блоков
    BaseNode::Ptr v1, v2, v3;
    BaseNode::Ptr v4, v5, v6;

    // Ограничивающие кривые области
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;
    Curve::Ptr wedge;

    // Блочная структура
    BlockStructured m_blocks;
};

} // namespace zephyr::geom::generator::collection
