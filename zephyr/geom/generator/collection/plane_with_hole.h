#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/mesh/storage.h>

#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block_structured.h>

namespace zephyr::geom::generator::collection {

using zephyr::geom::Boundary;

/// @class PlaneWithHole. Генератор для создания сетки внутри прямоугольника с отверстием.
class PlaneWithHole : public BlockStructured {
public:

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit PlaneWithHole(YAML::Node config);
#endif

    struct Boundaries {
        Boundary left;
        Boundary right;
        Boundary bottom;
        Boundary top;
        Boundary hole;
    };

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc, yc, r Центр и радиус отверстия
    PlaneWithHole(double xmin, double xmax, double ymin, double ymax,
                  double xc, double yc, double r);


    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по оси Oy подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_nx(int nx);

    void set_boundaries(const Boundaries& flags);

    /// @brief Создать указатель на соответствующий volumes::Box
    Box bbox() const final;

private:

    void check_params() const override;

    void init_blocks();


    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xc, m_yc, m_r;

    /// @brief Безразмерный параметр, соотношение радиусов "внешней оружности"
    /// и внутренней (отверстия)
    double m_xi;

    Boundary m_left_flag, m_right_flag;
    Boundary m_bottom_flag, m_top_flag;
    Boundary m_hole_flag;


    // Куча параметров
    // Задаем базисные вершины для струтурированных блоков
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
