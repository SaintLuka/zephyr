#pragma once

#include <zephyr/configuration.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @class Sector. Генератор для создания блочно-структурированной
/// сетки внутри сектора или круга.
class Sector : public Generator {
public:
    using Ptr = std::shared_ptr<Sector>;

    /// @brief Флаги граничных условий
    struct Boundaries {
        Boundary outer = Boundary::ZOE;  ///< На внешней окружности
        Boundary left  = Boundary::WALL; ///< На левой границе (если есть)
        Boundary right = Boundary::WALL; ///< На правой границе (если есть)
        Boundary inner = Boundary::WALL; ///< На внутренней границе (если есть)
    };

#ifdef ZEPHYR_YAML
    /// @brief Конструктор класса по кофигу
    explicit Sector(YAML::Node config);
#endif

    /// @brief Генератор сетки внутри круга из четырехугольников
    /// @param r_max Радиус круга
    /// @param r_min Радиус внутренней части
    /// @param angle Угол раствора: от 0 до 2п
    /// @param hole Оставлять ли дырку в центре
    Sector(double r_max, double r_min, double angle, bool hole = false);

    /// @brief Создать указатель на класс
    template <class... Args>
    static Sector::Ptr create(Args&&... args){
        return std::make_shared<Sector>(std::forward<Args>(args)...);
    }

    /// @brief Установить желаемое число ячеек по углу
    void set_n_phi(int N);

    /// @brief Установить флаги граничных условий
    void set_boundaries(Boundaries bounds);

    /// @brief Количество ячеек сетки
    int size() const final;

    /// @brief Ограничивающий объем
    Box bbox() const final;

    /// @brief Создать сетку общего вида
    Grid make() final;

private:
    /// @brief Подобрать параметры сетки исходя из геометрии
    /// и требуемого размера сетки
    void init_params();

    /// @brief Проверить на корявый ввод
    void check_params() const final;

    /// @brief Повернуть вектор несколько раз на 90 градусов
    /// @param r Число поворотов от 0 до 3 (включительно)
    Vector3d rotate(const Vector3d& v, int r) const;

    /// @brief Углы кратные pi/2 дотягиваются
    /// до точных значений
    void fit_angle();

    double m_r2;      ///< Внешний радиус
    double m_r1;      ///< Радиус внутренней части
    double m_angle;   ///< Угол раствора
    bool m_hole;      ///< Оставлять дырку

    Boundaries m_bounds;  ///< Флаги граничных условий

    double m_part_angle;  ///< Раствор угла / m_parts
    double m_part_cos;    ///< Косинус угла
    double m_part_sin;    ///< Косинус угла

    double m_zeta;    ///< Параметр при построении сетки

    int m_parts;   ///< Количество частей
    int m_nx;      ///< Разбиение больших неструктурированных квадратов
    int m_Nx;      ///< Разбиение по углу в структурированной части (Nx = 2*m_parts*nx)
    int m_ny;      ///< Разбиение по радиусу во внутренней части
    int m_Ny;      ///< Разбиение по радиусу во внешней структурированной части
    int m_size;    ///< Полное число ячеек
};

} // namespace zephyr::geom::generator
