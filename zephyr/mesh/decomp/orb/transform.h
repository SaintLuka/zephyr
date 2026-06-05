#pragma once

#include <string>

#include <zephyr/geom/box.h>

namespace zephyr::mesh::decomp {

using geom::Vector3d;
using geom::Box;

/// @brief Пределы отображения
struct Limits {
    Vector3d min, max;

    /// @brief Бесконечные пределы по умолчанию
    Limits();

    /// @brief Ограничения для выбранного типа, к примеру, для полярных r > 0.
    explicit Limits(std::string type);

    /// @brief Конструктор
    Limits(const Vector3d &vmin, const Vector3d &vmax);

    /// @brief
    void set(char c, int axes);

    /// @brief (-inf, inf) по оси axes
    void set_coord(int axes);

    /// @brief [0, inf) по оси axes
    void set_radius(int axes);

    /// @brief [-pi, pi] по оси axes
    void set_angle(int axes);
};

/// @brief Преобразование координат в прямоугольник
class Transform {
public:
    /// @brief Конструктор по типу декомпозиции
    /// @param type Полный тип, если двумерная область то 2 символа,
    /// если трёхмерная область, то три символа!
    // explicit Transform(const std::string &type = "XYZ");

    /// @brief Трехмерное отображение, бесконечная область
    Transform();

    /// @brief Конструктор по типу декомпозиции
    /// @param type Допускается неполный тип (к примеру, "X" или "R").
    /// @param domain Область декомпозиции (для проверки размерности).
    Transform(const Box& domain, const std::string& type);

    /// @brief Тип декомпозиции
    const std::string &type() const;

    /// @brief Размерность декомпозиции
    int dim() const;

    /// @brief Прямое отображение
    Vector3d operator()(const Vector3d &) const;

    /// @brief Прямое отображение
    Vector3d mapping(const Vector3d &) const;

    /// @brief Обратное отображение
    Vector3d inverse(const Vector3d &) const;

    /// @brief Границы после отображения
    const Box &box() const { return m_box; }

    /// @brief Нижняя граница по оси
    double x_min() const { return m_box.vmin.x(); }
    double y_min() const { return m_box.vmin.y(); }
    double z_min() const { return m_box.vmin.z(); }

    /// @brief Верхняя граница по оси
    double x_max() const { return m_box.vmax.x(); }
    double y_max() const { return m_box.vmax.y(); }
    double z_max() const { return m_box.vmax.z(); }

    /// @brief Ширина по оси
    double x_width() const { return m_box.sizes().x(); }
    double y_width() const { return m_box.sizes().y(); }
    double z_width() const { return m_box.sizes().z(); }

    /// @brief Доступные варианты декомпозиции
    static std::vector<std::string> available_types();

    /// @brief Тестирует преобразования
    static void check();

    /// @brief Установить границы отображения (domain - до отображения)
    void update_limits(const Box& domain);

private:
    void init_(const std::string& type);

    using map_function = std::function<Vector3d(const Vector3d&)>;

    std::string  m_type;     ///< Тип отображения/декомпозиции
    map_function m_mapping;  ///< Прямое отображение
    map_function m_inverse;  ///< Обратное отображение
    Box          m_box;      ///< Границы (после отображения)
};

} // namespace zephyr::mesh::decomp