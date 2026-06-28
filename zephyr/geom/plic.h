#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::mesh {
class EuCell;
}

namespace zephyr::geom {

/// @brief Восстановление интерфейса для многоматериальных задач.
class Plic {
public:
    enum Type : int {
        GRAD,    ///< Простой градиент, для произвольных сеток
        PnY,     ///< Parker & Youngs, для декартовых
        ELVIRA,  ///< ELVIRA, для декартовых
        CSIR,    ///< CSIR, для декартовых, учитывает размерность
        CSIR_2D, ///< CSIR_2D, для декартовых, простая формула как для 2D
        // LSF,     ///< Least-squares fit, итеративная для произвольных сеток
    };

    /// @brief Пара (p, n), где p - расстояние от центра ячейки
    /// до плоскости со знаком, n - нормаль плоскости
    using plane_t = std::tuple<double, Vector3d>;

    /// @brief Функция вытаскивания i-ой объемной доли из ячейки
    using get_fraction_t = std::function<double(const mesh::EuCell&, int)>;

    /// @brief Тривиальный конструктор
    Plic();

    /// @brief Установить все свойства
    /// @param get_vf Функция получения объемной доли из ячейки
    Plic(int dim, bool cartesian, Type type, const get_fraction_t& get_vf);

    /// @brief Реконструировать плоскость в ячейке
    plane_t plane(mesh::EuCell& cell, int idx) const;

private:
    /// @brief Основная функция, инициализируется после выбора параметров
    std::function<plane_t(mesh::EuCell&, int)> m_find_plane;
};

} // namespace zephyr::geom