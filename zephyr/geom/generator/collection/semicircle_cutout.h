#pragma once

#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

/// @brief Генератор сетки внутри прямоугольной области с полукруглым вырезом.
class SemicircleCutout final : public Generator {
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
    SemicircleCutout(double xmin, double xmax, double ymin, double ymax, double xc, double r);

    /// @brief Создать указатель на класс
    template <class... Args>
    static SemicircleCutout::Ptr create(Args&&... args){
        return std::make_shared<SemicircleCutout>(std::forward<Args>(args)...);
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
    double m_xc, m_r;

    // Базисные вершины
    BaseNode::Ptr v1, v2, v3, v4;
    BaseNode::Ptr v5, v6, v7, v8;
    BaseNode::Ptr v9, v10, v11, v12;
    BaseNode::Ptr v13, v14, v15, v16;

    // Ограничивающие кривые области
    Curve::Ptr left;
    Curve::Ptr right;
    Curve::Ptr bottom;
    Curve::Ptr top;
    Curve::Ptr circle;

    // Блочная структура
    BlockStructured m_blocks;
};

} // namespace zephyr::geom::generator::collection
