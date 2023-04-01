#pragma once

#include <array>
#include <map>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/generator/generator.h>

namespace zephyr { namespace mesh { namespace generator {

using VerticesList = std::vector<std::array<Vector3d, 4>>;

/// @brief Генератор двумерной блочно-структурированной сетки из
/// четырехугольных элементов.
class BlockStructured : public Generator {
private:
    using Curve_Ptr = std::shared_ptr<Curve>;
    using Curve_Ref = const std::shared_ptr<Curve> &;

public:
    /// @brief Конструктор
    /// @param n_blocks Число блоков, нельзя изменить в дальнейшем
    explicit BlockStructured(int n_blocks);

    /// @brief Ссылка на структурированный блок по индексу
    Block &operator[](int idx);

    /// @brief Ссылка на структурированный блок по индексу
    const Block &operator[](int idx) const;

    /// @brief Связать структурированнык блоки.
    /// Добавленные блоки должны опираться на одинаковые базисные вершины.
    void link();

    /// @brief Установить точность сглаживания
    /// @param eps Безразмерная величина меньше единицы
    void set_accuracy(double eps);

    /// @brief Вывод информации при генерации
    void set_verbose(bool v = true);

    /// @brief Количество ячеек
    int size() const override;

    /// @brief Инициализировать хранилище
    void initialize(Storage& cells) override;

private:
    /// @brief Проверяет размеры блоков
    void check_consistency() const;

    /// @brief Создать вершины будущей сетки
    VerticesList get_vertices();

    /// @brief Список граничных условий
    BoundaryFlags get_boundary_flags();

    /// @brief Массив блоков
    std::vector<Block> m_blocks;

    /// @brief Точность сглаживания
    double m_epsilon;

    /// @brief Вывод информации при генерации
    bool m_verbose;
};

} // namepsace generator
} // namespace mesh
} // namespace zephyr
