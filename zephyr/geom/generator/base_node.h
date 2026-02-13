#pragma once

#include <memory>
#include <vector>
#include <set>

#include <zephyr/geom/vector.h>

namespace zephyr::geom::generator {

class Block;

/// @brief Базисная вершина блока. Через базисные вершины проходят границы
/// области, на базисных вершинах строятся структурированные блоки.
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

    /// @brief Положение узла в виде Vector3d
    const Vector3d &pos() const { return m_pos; }

    /// @brief X-координата узла
    double x() const { return m_pos.x(); }

    /// @brief Y-координата узла
    double y() const { return m_pos.y(); }

    /// @brief Изменить состояние (доступно в редактируемом состоянии)
    void set_fixed(bool fixed);

    /// @brief Установить положение (доступно в редактируемом состоянии)
    void set_pos(double x, double y);

    /// @brief Установить положение (доступно в редактируемом состоянии)
    void set_pos(const Vector3d &pos);

    /// @brief Очистить списки и перевести узел в редактируемое состояние
    void clear();

    /// @brief Передать множество смежных блоков и завершить редактирование
    /// базисной вершины. Функция достраивает массивы смежных вершин и блоков,
    /// упорядочивает их в правильном порядке и т.д.
    void finalize(const std::set<Block_Ptr>& blocks);

    // --------------- Функции доступны только после finalize() ---------------

    /// @brief Вершина лежит на границе?
    bool is_boundary() const;

    /// @brief Внутренняя вершина считается регулярной, если имеет четыре
    /// смежных вершины и четыре смежных блока. Вершина на границе является
    /// регулярной, если имеет три смежных вершины и два смежных блока.
    bool regular() const;

    /// @brief Сингулярная вершина? (любая не регулярная)
    bool singular() const { return !regular(); }

    /// @brief Степень вершины (число смежных узлов).
    int degree() const;

    /// @brief Число смежных узлов (== degree)
    int n_adjacent_nodes() const;

    /// @brief Число смежных блоков
    int n_adjacent_blocks() const;

    /// @brief Массив смежных базисных вершин. Узлы обходятся против часовой
    /// стрелки внутри области. Если есть смежный узел в положительном
    /// направлении оси Ox, то он будет первым.
    const std::vector<BaseNode::WPtr>& adjacent_nodes() const;

    /// @brief Массив смежных блоков. Число смежных блоков меньше или равно
    /// числу смежных узлов. Блоки обходятся против часовой стрелки внутри
    /// области: i-ый блок располагается между i-ой и (i+1)-ой смежными
    /// вершинами. Если смежные узлы расположены по осям координат, тогда
    /// блоки нумеруются как квадранты плоскости.
    const std::vector<Block_WPtr>& adjacent_blocks() const;

private:
    // @brief Узел в редактируемом состоянии
    bool m_editable{true};

    /// @brief Точка может быть подвижной или фиксированной. Угловые точки
    /// области фиксируются автоматически. Можно дополнительно фиксировать
    /// точки на границе и внутри области.
    bool m_fixed{false};

    /// @brief Положение узла
    Vector3d m_pos{Vector3d::Zero()};

    /// @brief Узел на границе сетки? Проверка по списку смежных блоков
    bool m_boundary{false};

    /// @brief Смежные блоки (в завершенном состоянии)
    std::vector<Block_WPtr> m_adjacent_blocks{};

    /// @brief Смежные вершины (в завершенном состоянии)
    std::vector<BaseNode::WPtr> m_adjacent_nodes{};
};

} // namespace zephyr::geom::generator
