#pragma once

#include <vector>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {


/// @brief Генератор двумерной блочно-структурированной сетки из
/// четырехугольных элементов.
class BlockStructured : public Generator {
public:
    /// @brief Конструктор
    BlockStructured(int n_blocks = 0);

    /// @brief Увеличить (!) число блоков
    void resize(int n_blocks);

    /// @brief Ссылка на структурированный блок по индексу
    Block &operator[](int idx);

    /// @brief Связать структурированные блоки.
    /// Добавленные блоки должны опираться на одинаковые базисные вершины.
    void link();

    /// @brief Построить блоки
    void plot() const;

    /// @brief Инициализировать хранилище
    Grid make() override;

    /// @brief Установить точность сглаживания
    /// @param eps Безразмерная величина меньше единицы
    void set_accuracy(double eps);

    /// @brief Вывод информации при генерации
    void set_verbose(bool v = true);

    /// @brief Установить условие на узлы, которые не сглаживаются
    /// @param func Функция, возвращает true для неподвижных узлов.
    void set_fixed(std::function<bool(const Vector3d&)> func) {
        m_fixed = func;
    }

public:
    enum Stage {
        EDITABLE,
        LINKED,
        OPTIMIZED,
    };

    void remove_null_blocks();

    /// @brief Количество ячеек
    int calc_cells() const;

    /// @brief Максимальное число узлов
    int calc_nodes() const;

    void initialize();
    void optimize();

    /// @brief Проверяет размеры блоков
    void check_consistency() const;

    /// @brief Стадия редактирования
    Stage m_stage = EDITABLE;

    /// @brief Массив блоков
    std::vector<Block::Ptr> m_blocks{};

    /// @brief Точность сглаживания
    double m_epsilon{1.0e-3};

    /// @brief Вывод информации при генерации
    bool m_verbose{false};

    /// @brief Условие стационарных узлов
    std::function<bool(const Vector3d&)> m_fixed;
};

} // namespace zephyr::geom::generator
