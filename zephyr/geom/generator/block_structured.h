#pragma once

#include <vector>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Генератор двумерной блочно-структурированной сетки из
/// четырехугольных элементов.
class BlockStructured : public Generator {
public:
    using Ptr = std::shared_ptr<BlockStructured>;

    /// @brief Конструктор
    /// @param n_blocks Число блоков, нельзя изменить в дальнейшем
    explicit BlockStructured(int n_blocks);

    /// @brief Создать указатель на класс
    template <class... Args>
    static BlockStructured::Ptr create(Args&&... args){
        return std::make_shared<BlockStructured>(std::forward<Args>(args)...);
    }

    /// @brief Ссылка на структурированный блок по индексу
    Block &operator[](int idx);

    /// @brief Ссылка на структурированный блок по индексу
    const Block &operator[](int idx) const;

    /// @brief Связать структурированные блоки.
    /// Добавленные блоки должны опираться на одинаковые базисные вершины.
    void link();

    /// @brief Установить точность сглаживания
    /// @param eps Безразмерная величина меньше единицы
    void set_accuracy(double eps);

    /// @brief Вывод информации при генерации
    void set_verbose(bool v = true);

    /// @brief Количество ячеек
    int size() const final;

    /// @brief Инициализировать хранилище
    Grid make() override;

protected:
    /// @brief Количество ячеек
    int calc_cells() const;

    /// @brief Максимальное число узлов
    int calc_nodes() const;

    void initialize();

    /// @brief Проверяет размеры блоков
    void check_consistency() const;

    /// @brief Массив блоков
    std::vector<Block> m_blocks;

    /// @brief Точность сглаживания
    double m_epsilon;

    /// @brief Вывод информации при генерации
    bool m_verbose;
};

} // namespace zephyr::geom::generator
