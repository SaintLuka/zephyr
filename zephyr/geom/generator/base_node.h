#pragma once

#include <memory>
#include <vector>
#include <set>

#include <zephyr/geom/vector.h>

namespace zephyr::geom::generator {

class Block;

/// @brief Базисный узел блока. Через базисные узлы проходят границы области,
/// на базисных узлах строятся структурированные блоки.
class BaseNode {
    using Block_Ptr = std::shared_ptr<Block>;
    using Block_Ref = const std::shared_ptr<Block>&;
    using Block_WPtr = std::weak_ptr<Block>;
public:
    using Ptr = std::shared_ptr<BaseNode>;
    using Ref = const std::shared_ptr<BaseNode>&;
    using WPtr = std::weak_ptr<BaseNode>;

    /// @brief Конструктор
    /// @param fixed Неподвижный узел?
    BaseNode(const Vector3d &v, bool fixed);

    /// @brief Создать умный указатель
    /// @param fixed Неподвижный узел?
    static BaseNode::Ptr create(const Vector3d &v, bool fixed);

    /// @brief Создать умный указатель
    /// @param x, y Координаты точки
    /// @param fixed Неподвижный узел?
    static BaseNode::Ptr create(double x, double y, bool fixed);

    /// @brief В редактируемом состоянии?
    bool editable() const { return m_editable; }

    /// @brief В завершенном состоянии?
    bool finalized() const { return !m_editable; }

    /// @brief Неподвижный узел?
    bool is_fixed() const { return m_fixed; }

    /// @brief X-координата узла
    double x() const { return m_pos.x(); }

    /// @brief Y-координата узла
    double y() const { return m_pos.y(); }

    /// @brief Ссылка на содержимое в виде Vector3d
    const Vector3d &pos() const { return m_pos; }

    /// @brief Изменить состояние (в редактируемом состоянии)
    void set_fixed(bool fixed);

    /// @brief Установить положение (в редактируемом состоянии)
    void set_pos(const Vector3d &pos);

    /// @brief Установить положение (в редактируемом состоянии)
    void set_pos(double x, double y);

    /// @brief Перевести в редактируемое состояние
    void reset();

    /// @brief ???
    void finalize(const std::set<std::shared_ptr<Block>>& blocks);

    /// @brief Вершина лежит на границе?
    /// Проверка по списку смежных блоков
    bool is_boundary() const;

    /// @brief Внутренняя вершина считается регулярной, если имеет четыре
    /// смежных вершины и четыре смежных блока. Вершина на границе является
    /// регулярной, если имеет три смежных вершины и два смежных блока.
    bool regular() const;

    /// @brief Сингулярная вершина?
    bool singular() const { return !regular(); }

    /// @brief Степень вершины (число смежных узлов)
    /// @details Для базисной вершины внутри области на границе области нормальным является
    /// наличие двух смежных блоков (одного блока для угловой вершины),
    /// для базисной вершины внутри области нормальное число смежных блоков
    /// равно четырем. В этих случаях вершина считается регулярной,
    /// иначе - сингулярной.
    int degree() const;

    /// @brief Число смежных блоков
    int n_adjacent_blocks() const;

    /// @brief Массив смежных блоков
    /// Блоки обходятся против часовой стрелки внутри области. i-ый блок располагается
    /// между i-ой и (i+1) смежной вершиной. Если смежные узлы расположены по осям
    /// координат, тогда блоки нумеруются как квадранты плоскости.
    const std::vector<Block_WPtr>& adjacent_blocks() const;

    /// @brief Число смежных узлов (= degree)
    int n_adjacent_nodes() const;

    /// @brief Массив смежных базисных вершин
    /// Узлы обходятся против часовой области внутри области. Если есть смежный
    /// узел в положительном направлении оси Ox, то он будет первым.
    const std::vector<BaseNode::WPtr>& adjacent_nodes() const;

private:
    // @brief В редактируемом состоянии
    bool m_editable{true};

    /// @brief Точка может быть подвижной или фиксированной. Угловые точки
    /// области фиксируются автоматически. Можно дополнительно фиксировать
    /// точки на границе и внутри области.
    bool m_fixed{false};

    /// @brief Положение узла
    Vector3d m_pos{Vector3d::Zero()};

    /// @brief Внутренний узел или внешний?
    bool m_boundary{false};

    /// @brief Смежные блоки (в завершенном состоянии)
    std::vector<Block_WPtr> m_adjacent_blocks{};

    /// @brief Смежные вершины (в завершенном состоянии)
    std::vector<BaseNode::WPtr> m_adjacent_nodes{};
};

} // namespace zephyr::geom::generator
