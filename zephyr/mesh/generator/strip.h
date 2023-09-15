#pragma once

#include <array>

#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/generator/generator.h>


namespace zephyr { namespace mesh { namespace generator {

/// @class Strip. Простой класс для создания Storage для хранения
/// квазиодномерной сетки (прямоугольная сетка толщиной в одну ячейку).
class Strip : public Generator {
public:

    enum class Nodes {
        UNIFORM,
        SINE,
        EXP,
        RANDOM
    };

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit Strip(YAML::Node config);
#endif

    /// @brief Конструктор класса
    /// @param xmin, xmax Границы прямоугольника по оси x
    Strip(double xmin, double xmax, Nodes nodes = Nodes::UNIFORM);

    /// @brief Установить желаемое число ячеек сетки по оси Ox
    void set_nx(int nx);

    /// @brief Установить желаемое число ячеек сетки
    void set_size(int N);

    /// @brief Установить флаги граничных условий
    void set_boundary_flags(FaceFlag left, FaceFlag right);

    double x_min() const;

    double x_max() const;

    double y_min() const;

    double y_max() const;

    int nx() const;

    bool periodic_along_x() const;

    Box bbox() const final;

    /// @brief Инициализация переданной сетки
    void initialize(Storage& storage) final;

private:

    /// @brief Проверить параметры сетки перед созданием
    void check_params() const final;

    ///@brief Соотношение сторон прямоугольника
    const double aspect = 0.1;

    /// @brief Тип задания узлов сетки
    Nodes m_nodes;

    double m_xmin, m_xmax;

    FaceFlag m_left_flag;
    FaceFlag m_right_flag;
};

} // generator
} // mesh
} // zephyr
