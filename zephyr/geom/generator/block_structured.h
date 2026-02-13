#pragma once

#include <vector>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom::generator {

/// @brief Опции оптимизатора
struct optimize_options {
    /// @brief Число стадий оптимизации, столько раз вызовется optimize()
    int steps{2};

    /// @brief Начальное число ячеек по оси минимальной ячейки
    int N{5};

    /// @brief Уровень вывода информационных сообщений (0 - нет вывода)
    int verbose{0};

    /// @brief Относительная погрешность (критерий остановки сглаживания)
    double eps{1.0e-4};
};

/// @brief Пара связанных базисных узлов
struct node_pair_t {
    BaseNode::WPtr v1, v2;

    const void* ptr1() const { return v1.lock().get(); }
    const void* ptr2() const { return v2.lock().get(); }

    bool operator==(const node_pair_t& other) const {
        return (ptr1() == other.ptr1() && ptr2() == other.ptr2()) ||
               (ptr1() == other.ptr2() && ptr2() == other.ptr1());
    }
};

/// @brief Пара смежных блоков
struct block_pair_t {
    Block::WPtr b1{}, b2{};
    Side side1{}, side2{};

    /// @brief На границе только один блок
    bool boundary() const { return b2.expired(); }
};

} // zephyr::geom::generator

namespace std {
template <>
struct hash<zephyr::geom::generator::node_pair_t> {
    size_t operator()(const zephyr::geom::generator::node_pair_t& pp) const noexcept {
        hash<const void*> hasher;
        size_t h1 = hasher(pp.ptr1());
        size_t h2 = hasher(pp.ptr2());
        return (h1 + h2) * (h1 + h2 + 1) / 2 + std::min(h1, h2);
    }
};
} // namespace std

namespace zephyr::geom::generator {

/// @brief Генератор двумерной блочно-структурированной сетки из
/// четырехугольных элементов.
class BlockStructured : public Generator {
public:
    /// @brief Конструктор
    BlockStructured();

    /// @brief Ссылка на структурированный блок по индексу
    Block &operator[](int idx) const;

    /// @brief Число блоков
    int size() const { return static_cast<int>(m_blocks.size()); }

    /// @brief Добавить блок
    void operator+=(const std::array<BaseNode::Ptr, 4>& base_nodes);

    /// @brief Установить границу на сторону (v1, v2)
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve);

    /// @brief Установить границу для цепочки узлов (v1, v2, v3, ...)
    void set_boundary(std::initializer_list<BaseNode::Ptr> nodes, Curve::Ref curve);

    /// @brief Установить условие на узлы, которые не сглаживаются
    /// @param func Функция, возвращает true для неподвижных узлов.
    void set_fixed(std::function<bool(const Vector3d&)> func) {
        m_fixed = func;
    }

    /// @brief Полный цикл оптимизаций. Последовательный вызов optimize на сетках
    /// различного размера, в зависимости от переданных настроек.
    void optimize(const optimize_options& opts = {});

    /// @brief Можно вызвать несколько раз подряд
    void optimize(int N, double eps = 1.0e-3, int verbose = 0);

    /// @brief Построить блоки с разбиением
    void plot() const;

    /// @brief Сгенерировать сетку
    Grid make() override;

protected:
    enum Stage {
        EDITABLE,
        LINKED,
        OPTIMIZED,
    };

    /// @brief Построить блоки
    void plot_debug() const;

    /// @brief Удалить неактуальные элементы
    void remove_redundant();

    /// @brief Связать структурированные блоки.
    /// Добавленные блоки должны опираться на одинаковые базисные вершины.
    void link_blocks();

    /// @brief Количество ячеек
    int calc_cells() const;

    /// @brief Максимальное число узлов
    int calc_nodes() const;

    /// @brief Проставить размеры блоков (от конформного модуля)
    void setup_block_sizes(int N);

    /// @brief Первый запуск оптимизации
    void optimize_init(int N, double eps = 1.0e-3, int verbose = 0);

    /// @brief Последующие запуски оптимизации
    void optimize_again(int N, double eps = 1.0e-3, int verbose = 0);


    /// @brief Стадия редактирования
    Stage m_stage = EDITABLE;

    /// @brief Массив блоков
    std::vector<Block::Ptr> m_blocks{};

    /// @brief Множество рёбер
    std::unordered_map<node_pair_t, block_pair_t> m_edges;

    /// @brief Условие стационарных узлов
    std::function<bool(const Vector3d&)> m_fixed;
};

} // namespace zephyr::geom::generator
