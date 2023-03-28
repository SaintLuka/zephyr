#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/mesh/storage.h>

#ifdef ZEPHYR_ENABLE_MPI
#include <zephyr/network/mpi.h>
#endif

#include <zephyr/mesh/generator/base.h>
#include <zephyr/mesh/generator/vertex.h>
#include <zephyr/mesh/generator/curve/curve.h>
#include <zephyr/mesh/generator/block_structured.h>

namespace zephyr { namespace mesh { namespace generator { namespace collection {

using zephyr::geom::FaceFlag;

/// @class PlaneWithHole. Генератор для создания сетки внутри прямоугольника с отверстием.
class PlaneWithHole : public generator::Base {
public:

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit PlaneWithHole(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax, ymin, ymax Границы прямоугольника
    /// @param xc, yc, r Центр и радиус отверстия
    PlaneWithHole(double xmin, double xmax, double ymin, double ymax,
                  double xc, double yc, double r);


    /// @brief Установить желаемое число ячеек сетки по оси Ox
    /// @details Число ячеек по оси Oy подбирается так, чтобы aspect ячеек
    /// был около единицы
    void set_nx(int nx);

    /// @brief Размер сетки
    int size() const override;

    /// @brief Создать указатель на соответствующий volumes::Box
    Box bbox() const final;

private:

    void check_params() const override;

    void init_blocks();

    void initialize(Storage& storage) final;



    BlockStructured blocks;

    double m_xmin, m_xmax;
    double m_ymin, m_ymax;
    double m_xc, m_yc, m_r;

    /// @brief Безразмерный параметр, соотношение радиусов "внешней оружности"
    /// и внутренней (отверстия)
    double m_xi;

    FaceFlag m_left_flag, m_right_flag;
    FaceFlag m_bottom_flag, m_top_flag;
    FaceFlag m_hole_flag;


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

} // collection
} // generator
} // mesh
} // zephyr
