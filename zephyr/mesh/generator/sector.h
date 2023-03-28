#pragma once

#include <array>
#include <vector>

#include <zephyr/configuration.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/generator/base.h>


namespace zephyr { namespace mesh { namespace generator {


/// @class Sector. Простой генератор для создания Storage для хранения
/// блочно-структурированной сетки внутри сектора/круга
class Sector : public Base {
public:

#ifdef ZEPHYR_ENABLE_YAML
    /// @brief Конструктор класса по кофигу
    explicit Sector(YAML::Node config);
#endif

    /// @brief Генератор сетки внутри круга из четырехугольников
    /// @param r_max Радиус круга
    /// @param r_min Радиус внутренней части
    /// @param angle Угол раствора: от 0 до 2п
    /// @param hole Оставлять ли дырку в центре
    Sector(double r_max, double r_min, double angle, bool hole = false);

    /// @brief Установить желаемое число ячеек по углу
    void set_n_phi(int N);

    /// @brief Установить желаемое число ячеек сетки.
    /// Точное количество ячеек подбирается алгоритмически.
    void set_size(int N);

    /// @brief Установить флаги граничных условий
    /// @param outer На внешней окружности
    /// @param left На левой границе
    /// @param right На правой границе
    /// @param inner На внутренней границе
    void set_boundary_flags(FaceFlag outer, FaceFlag left,
            FaceFlag right, FaceFlag inner = FaceFlag::UNDEFINED);

    Box bbox() const final;

    void initialize(Storage &cells) final;

private:

    void init_params();

    void check_params() const final;

    Vector3d rotate(const Vector3d& v, int r) const;

    void fit_angle();

    void create_vertices();

    double m_r2;      ///< Внешний радиус
    double m_r1;      ///< Радиус внутренней части
    double m_angle;   ///< Угол раствора
    bool m_hole;      ///< Оставлять дырку

    FaceFlag m_left_flag;   ///< Флаг граничных условий
    FaceFlag m_right_flag;  ///< Флаг граничных условий
    FaceFlag m_inner_flag;  ///< Флаг граничных условий
    FaceFlag m_outer_flag;  ///< Флаг граничных условий

    double m_part_angle;    ///< Раствор угла / m_parts
    double m_part_cos;      ///< Косинус угла
    double m_part_sin;      ///< Косинус угла

    double m_zeta;    ///< Параметр при построении сетки

    int m_parts;   ///< Количество частей
    int m_nx;      ///< Разбиение больших неструктурированных квадратов
    int m_Nx;      ///< Разбиение по углу в структурированной части (Nx = 2*m_parts*nx)
    int m_ny;      ///< Разбиение по радиусу во внутренней части
    int m_Ny;      ///< Разбиение по радиусу во внешней структурированной части

    /// @brief Массив вершин, каждая строка описывает ячейку
    std::vector<std::array<Vector3d, 4>> m_vertices;
};

} // generators
} // mesh
} // zephyr
