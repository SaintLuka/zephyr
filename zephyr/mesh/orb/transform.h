#pragma once

#include <string>
#include <Vector3d>

#include <zephyr/data/type/Vector3d.h>

namespace zephyr { namespace network { namespace decomposition {

using zephyr::data::Vector3d;


/// @brief Пределы отображения
class Limits {
public:
    Vector3d min, max;

    /// @brief Конструктор по умолчанию,
    /// бесконечные границы
    Limits();

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
struct Transform {
public:
    /// @brief Конструктор по типу декомпозиции
    explicit Transform(const std::string &type);

    /// @brief Тип декомпозиции
    const std::string &type() const;

    /// @brief Размерность декомпозиции
    size_t dimension() const;

    /// @brief Границы отображения
    const Limits &limits() const;

    /// @brief Прямое отображение
    Vector3d operator()(const Vector3d &) const;

    /// @brief Прямое отображение
    Vector3d mapping(const Vector3d &) const;

    /// @brief Обратное отображение
    Vector3d inverse(const Vector3d &) const;

    /// @brief Тестирует преобразования
    static void check();

private:
    typedef Vector3d (*map_function)(const Vector3d &);

    /// @brief Доступные варианты декомпозиции
    static const std::vector<std::string> available;

    /// @brief Тип декомпозиции
    std::string m_type;

    /// @brief Границы отображения
    Limits m_limits;

    /// @brief Прямое отображение
    map_function m_mapping;

    /// @brief Обратное отображение
    map_function m_inverse;
};

} // decomposition
} // network
} // zephyr