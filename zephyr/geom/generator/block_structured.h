#pragma once

#include <vector>

#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/generator/block.h>

namespace zephyr::geom::generator {

/// @brief Массив пар значений, присвоенных осям
template <typename T>
using Pairs = std::vector<AxisPair<T>>;

/// @brief Несколько таблиц с узлами
using Tables2D = std::vector<Table2D>;

/// @brief Опции оптимизатора
struct optimize_options {
    /// @brief Число стадий оптимизации, столько раз вызовется optimize()
    int steps{2};

    /// @brief Начальное число ячеек по оси минимальной ячейки
    int N{5};

    /// @brief Относительная погрешность (критерий остановки сглаживания)
    double eps{1.0e-4};

    /// @brief Использовать кэширование после оптимизации
    bool cache = true;
};

/// @brief Генератор двумерной блочно-структурированной сетки из
/// четырехугольных элементов. Этапы работы с классом:
///   1. Добавить блоки (оператор +=).
///   2. Установить граничные условия.
///   3. Задать число ячеек хотя бы вдоль одного ребра.
///   4. Вызвать оптимизацию (или не вызывать).
///   5. Вызвать make() и создать сетку.
///
/// Число ячеек (этап 3) обязательно задавать после этапа 2, то есть после
/// определения всей блочной структуры.
class BlockStructured : public Generator {
    static constexpr int default_verbosity{1};
    static constexpr int default_iters_count{3000};

public:
    /// @brief Конструктор
    BlockStructured();

    /// @brief Число блоков
    int n_blocks() const;

    /// @brief Отключить вывод информационных сообщений
    void disable_verbosity();

    /// @brief Задать уровень вывода информационных сообщений
    void set_verbosity(int verbose);

    /// @brief Задать число итераций сглаживания финальной сетки
    void set_iters_count(int iters_count);

    // ------------------------- Стадия редактирования ------------------------

    /// @brief Ссылка на структурированный блок по индексу
    Block &operator[](int idx) const;

    /// @brief Добавить блок
    void operator+=(const std::array<BaseNode::Ptr, 4>& nodes);

    /// @brief Граница, на которой лежит ребро (v1, v2), может быть nullptr
    Curve::Ptr get_boundary(BaseNode::Ref v1, BaseNode::Ref v2);

    /// @brief Установить границу на сторону (v1, v2)
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve);

    /// @brief Сделать прямолинейную границу (v1, v2).
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2);

    /// @brief Установить граничное условие на сторону (v1, v2).
    /// Если кривая не задана, то будет создана прямолинейная граница.
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Boundary boundary);

    /// @brief Установить границу для цепочки узлов (v1, v2, v3, ...)
    void set_boundary(const std::vector<BaseNode::Ptr>& nodes, Curve::Ref curve);

    /// @brief Сделать прямолинейную границу для цепочки узлов (v1, v2, v3, ...)
    /// Если точки лежат на одной прямой, то будет создана единственная граница.
    void set_boundary(const std::vector<BaseNode::Ptr>& nodes);

    /// @brief Установить граничное условие для цепочки узлов (v1, v2, v3, ...)
    /// Если кривые не заданы, то будет создана по возможности единственная
    /// прямолинейная граница, если такую прямую провести нельзя, то ломаная.
    void set_boundary(const std::vector<BaseNode::Ptr>& nodes, Boundary boundary);

    /// @brief Задать желаемое число ячеек на ребре.
    /// @param v1, v2 Пара базисных вершин, образующих ребро
    /// @param size Желаемое число ячеек на ребро
    void set_size(BaseNode::Ref v1, BaseNode::Ref v2, int size);

    /// @brief Задать желаемое число ячеек на набор граней.
    /// @param nodes Список узлов, которые образуют цепочку ребер
    /// @param size Число ячеек, которое должно приходиться на всю цепочку
    void set_size(const std::vector<BaseNode::Ptr>& nodes, int size);

    // ----------------------- Оптимизация (опционально) ----------------------

    /// @brief Полный цикл оптимизаций. Последовательный вызов optimize на сетках
    /// различного размера, в зависимости от переданных настроек.
    void optimize(const optimize_options& opts = {});

    /// @brief Построить базисные вершины и блоки
    void plot_layout() const;

    /// @brief Построить блоки с разбиением
    void plot() const;

    // -------------------- Финальная сетка и сглаживание ---------------------

    /// @brief Оценить количество ячеек финальной сетки
    int calc_cells() const;

    /// @brief Оценить количество узлов финальной сетки
    int calc_nodes() const;

    /// @brief Сгенерировать финальную сетку
    Grid make() const override;

protected:
    /// @brief Стадии существования класса
    enum Stage {
        EDITABLE,  ///< Редактирование (добавление блоков, граничных условий, числа ячеек)
        LINKED,    ///< Стадия после создания связей блоков, запрещает дальнейшее редактирование
        OPTIMIZED, ///< Стадия после вызова оптимизации, не является обязательной.
    };

    // -------------------------- Стадия связывания ---------------------------

    /// @brief Удалить неактуальные элементы
    void remove_redundant();

    /// @brief Связать структурированные блоки.
    /// Добавленные блоки должны опираться на одинаковые базисные вершины.
    void link_blocks();

    // ---------------------------- Размеры блоков ----------------------------

    /// @brief Количество ячеек
    int calc_cells(const Pairs<int>& sizes) const;

    /// @brief Максимальное число узлов
    int calc_nodes(const Pairs<int>& sizes) const;

    /// @brief Установить число ячеек вдоль оси блока
    void set_size(Pairs<int>& sizes, int b1, Axis axis, int N) const;

    /// @brief Проставить размеры блоков (от конформного модуля)
    Pairs<int> auto_block_sizes(int N) const;

    /// @brief Расставить размеры с приоритетом по желаемым
    Pairs<int> wanted_block_sizes() const;

    /// @brief Выставить начальные относительные размеры
    void init_rel_sizes(Pairs<double>& sizes, int b1) const;

    /// @brief Обновить размеры (из конформного модуля)
    /// @return true, если оба размера теперь заданы
    bool update_rel_sizes(Pairs<double>& sizes, int b1) const;

    /// @brief Установить относительный размер вдоль оси блока
    void set_rel_size(Pairs<double>& sizes, int b1, Axis axis, double N) const;

    /// @brief Распространить размеры на все блоки (должны быть определены
    /// модули, хотя бы у одного блока должен быть выставлен размер по одной
    /// из осей, блоки должны быть связаны)
    void propagate_sizes(Pairs<double>& sizes) const;

    // ---------------------- Создание внутренних узлов -----------------------

    /// @brief Можно вызвать несколько раз подряд
    void optimize(int N, double eps = 1.0e-3);

    /// @brief Шаг оптимизации.
    ///   Расстановка размеров, создание внутренних узлов и связей,
    ///   итерации сглаживания, установка конформного отображения.
    void optimization_step(int N, double eps = 1.0e-3);

    /// @brief Посчитать хэш сетки (вызывается до оптимизации)
    size_t get_hash(const optimize_options& opts) const;

    /// @brief Сохранить кэш-файл с оптимизированной сеткой
    void save_cache(const optimize_options& opts) const;

    /// @brief Загрузить оптимизированную сетку из кэша
    void load_cached();

    /// @brief Выставить консистентные конформные модули (после пересчета)
    void correct_modulus();

    // ------------------------ Информационные выводы -------------------------

    /// @brief Построить блоки, внутренние узлы и прочее
    void plot_debug() const;

    /// @brief Построить таблицы вершин
    static void plot(const Tables2D& all_vertices);

    /// @brief Строка информации о разбиении блока на ячейки
    void sizes_info(const Pairs<int>& sizes) const;

    /// @brief Строка информации о параметрах блоков
    void conformal_info(const Pairs<double>& lambda) const;

    // ----------------------------- Поля класса ------------------------------

    /// @brief Стадия редактирования
    Stage m_stage = EDITABLE;

    /// @brief Уровень вывода информационных сообщений (0 - нет вывода)
    int m_verbosity = default_verbosity;

    /// @brief Число итераций сглаживания финальной сетки
    int m_iters_count = default_iters_count;

    /// @brief Массив блоков
    std::vector<std::shared_ptr<Block>> m_blocks{};

    /// @brief Множество уникальных рёбер
    std::unordered_map<BaseEdge, BlockPair> m_edges;

    /// @brief Связанная цепочка вершин и желаемое число ячеек на эту цепочку
    struct wanted_size_t {
        std::vector<BaseNode::Ptr> nodes;
        int size;
    };

    /// @brief Желаемое число ячеек на рёбрах и цепочках рёбер.
    /// Имеет значение только при создании финальной сетки (вызове make()).
    std::vector<wanted_size_t> m_wanted_sizes{};

    // Внутренний узел + индексы и направления инцидентных блоков
    struct inner_node_t {
        BaseNode::Ptr node;
        std::vector<std::tuple<int, bool>> blocks;
    };

    /// @brief Внутренние узлы
    std::vector<inner_node_t> m_inner_nodes;
};

} // namespace zephyr::geom::generator
