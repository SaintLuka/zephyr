#pragma once

#include <vector>

#include <zephyr/geom/generator/array2d.h>
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
struct BaseEdge {
    BaseNode::WPtr v1, v2;

    const void* ptr1() const { return v1.lock().get(); }
    const void* ptr2() const { return v2.lock().get(); }

    bool operator==(const BaseEdge& other) const {
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

template <>
struct std::hash<zephyr::geom::generator::BaseEdge> {
    size_t operator()(const zephyr::geom::generator::BaseEdge& pp) const noexcept {
        hash<const void*> hasher;
        size_t h1 = hasher(pp.ptr1());
        size_t h2 = hasher(pp.ptr2());
        return (h1 + h2) * (h1 + h2 + 1) / 2 + std::min(h1, h2);
    }
};

namespace zephyr::geom::generator {

template <typename T>
using Pairs = std::vector<AxisPair<T>>;

using Tables2D = std::vector<Array2D<BsVertex::Ptr>>;


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

    /// @brief Проверить размеры смежных блоков перед созданием вершин
    void check_consistency(const Pairs<int>& sizes) const;

    /// @brief Полный цикл оптимизаций. Последовательный вызов optimize на сетках
    /// различного размера, в зависимости от переданных настроек.
    void optimize(const optimize_options& opts = {});

    /// @brief Можно вызвать несколько раз подряд
    void optimize(int N, double eps = 1.0e-3, int verbose = 0);

    /// @brief Построить блоки с разбиением
    void plot() const;

    /// @brief Построить блоки с разбиением
    static void plot(const Tables2D& all_vertices);

    /// @brief Сгенерировать сетку
    Grid make() override;

    /// @brief Сгенерировать узлы сетки с нуля
    Tables2D create_vertices(const Pairs<int>& sizes) const;

    /// @brief Сгенерировать узлы сетки на основе сохраненного отображения
    Tables2D create_vertices_again(const Pairs<int>& sizes) const;

    /// @brief Склеить узлы на границах блоков
    void merge_vertices(Tables2D& all_vertices, const Pairs<int>& sizes) const;

    /// @brief Связать узлы сетки
    void link_vertices(Tables2D& all_vertices, const Pairs<int>& sizes, Pairs<double>& lambda) const;

    /// @brief Сглаживание вершин в блоке
    static void smooth_vertices(const Tables2D& all_vertices);

    /// @brief Обновить положение вершин
    /// @return Максимальный относительный сдвиг вершин
    static double update_vertices(const Tables2D& all_vertices);

    /// @brief Строка информации о параметрах
    void conformal_info(const Pairs<double>& lambda) const;


    std::vector<AxisPair<double>> get_lambda_nan(const Pairs<int>& sizes) const;




    // ---------------------------- Размеры блока -----------------------------

    /// @brief Установить число ячеек вдоль оси блока
    void set_size(Pairs<int>& sizes, int b1, Axis axis, int N) const;

    /// @brief Строка информации о разбиении блока на ячейки
    void sizes_info(const Pairs<int>& sizes) const;

    /// @brief Выставить начальные относительные размеры
    void init_rel_sizes(Pairs<double>& sizes, int b1) const;

    /// @brief Обновить размеры (из конформного модуля)
    /// @return true, если оба размера теперь заданы
    bool update_rel_sizes(Pairs<double>& sizes, Block::Ref b) const;

    /// @brief Установить относительный размер вдоль оси блока
    void set_rel_size(Pairs<double>& sizes, int b1, Axis axis, double N) const;

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
    int calc_cells(const Pairs<int>& sizes) const;

    /// @brief Максимальное число узлов
    int calc_nodes(const Pairs<int>& sizes) const;

    /// @brief Проставить размеры блоков (от конформного модуля)
    Pairs<int> auto_block_sizes(int N);

    /// @brief Первый запуск оптимизации
    void optimize_init(int N, double eps = 1.0e-3, int verbose = 0);

    /// @brief Последующие запуски оптимизации
    void optimize_again(int N, double eps = 1.0e-3, int verbose = 0);


    /// @brief Стадия редактирования
    Stage m_stage = EDITABLE;

    /// @brief Массив блоков
    std::vector<Block::Ptr> m_blocks{};

    /// @brief Множество рёбер
    std::unordered_map<BaseEdge, block_pair_t> m_edges;

    /// @brief Условие стационарных узлов
    std::function<bool(const Vector3d&)> m_fixed;
};

} // namespace zephyr::geom::generator
