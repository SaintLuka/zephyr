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
    BaseNode(const Vector3d &v, bool fixed = false);

    /// @brief Создать умный указатель
    /// @param fixed Неподвижный узел?
    static BaseNode::Ptr create(const Vector3d &v, bool fixed = false);

    /// @brief Создать умный указатель
    /// @param x, y Координаты вершины
    /// @param fixed Неподвижный узел?
    static BaseNode::Ptr create(double x, double y, bool fixed = false);

    /// @brief В редактируемом состоянии?
    bool editable() const { return m_editable; }

    /// @brief В завершенном состоянии?
    bool finalized() const { return !m_editable; }

    /// @brief Неподвижная вершина?
    bool fixed() const { return m_fixed; }

    /// @brief Положение вершины
    const Vector3d &pos() const { return m_pos; }

    /// @brief X-координата
    double x() const { return m_pos.x(); }

    /// @brief Y-координата
    double y() const { return m_pos.y(); }

    /// @brief Изменить состояние (доступно в редактируемом состоянии)
    void set_fixed(bool fixed);

    /// @brief Установить положение (доступно в редактируемом состоянии)
    void set_pos(double x, double y);

    /// @brief Установить положение (доступно в редактируемом состоянии)
    void set_pos(const Vector3d &pos);

    /// @brief Очистить списки и перевести вершину в редактируемое состояние
    void clear();

    /// @brief Передать множество смежных блоков и завершить редактирование
    /// базисной вершины. Функция достраивает массивы смежных вершин и блоков,
    /// упорядочивает их в правильном порядке и т.д.
    void finalize(const std::set<Block_Ptr>& blocks);

    // --------------- Функции доступны только после finalize() ---------------

    /// @brief Вершина лежит на границе?
    bool boundary() const;

    /// @brief Внутренняя вершина считается регулярной, если имеет четыре
    /// смежных вершины и четыре смежных блока. Вершина на границе является
    /// регулярной, если имеет три смежных вершины и два смежных блока.
    bool regular() const;

    /// @brief Сингулярная вершина? (любая не регулярная)
    bool singular() const { return !regular(); }

    /// @brief Степень вершины (число смежных вершин).
    int degree() const;

    /// @brief Число смежных вершин (== degree)
    int n_adjacent_nodes() const;

    /// @brief Число смежных блоков
    int n_adjacent_blocks() const;

    /// @brief Массив смежных базисных вершин. Вершины обходятся против часовой
    /// стрелки внутри области. Если есть смежная вершина в положительном
    /// направлении оси Ox, то она будет первым.
    const std::vector<BaseNode::WPtr>& adjacent_nodes() const;

    /// @brief Массив смежных блоков. Число смежных блоков меньше или равно
    /// числу смежных вершин. Блоки обходятся против часовой стрелки внутри
    /// области: i-ый блок располагается между i-ой и (i+1)-ой смежными
    /// вершинами. Если смежные вершины расположены по осям координат, тогда
    /// блоки нумеруются как квадранты плоскости.
    const std::vector<Block_WPtr>& adjacent_blocks() const;

private:
    // @brief Вершина в редактируемом состоянии?
    bool m_editable{true};

    /// @brief Вершина может быть подвижной или фиксированной. Угловые точки
    /// области фиксируются автоматически. Можно дополнительно фиксировать
    /// точки на границе или внутри области.
    bool m_fixed{false};

    /// @brief Положение вершины
    Vector3d m_pos{Vector3d::Zero()};

    /// @brief Вершина на границе сетки? (по смежным блокам)
    bool m_boundary{false};

    /// @brief Смежные блоки (в завершенном состоянии)
    std::vector<Block_WPtr> m_adjacent_blocks{};

    /// @brief Смежные вершины (в завершенном состоянии)
    std::vector<BaseNode::WPtr> m_adjacent_nodes{};
};

/// @brief Пара смежных базисных вершин, образующих ребро
struct BaseEdge {
    BaseNode::WPtr v1, v2;

    const void* ptr1() const { return v1.lock().get(); }
    const void* ptr2() const { return v2.lock().get(); }

    /// @brief Сравнение пар вершин, порядок не важен
    bool operator==(const BaseEdge& other) const {
        return (ptr1() == other.ptr1() && ptr2() == other.ptr2()) ||
               (ptr1() == other.ptr2() && ptr2() == other.ptr1());
    }
};

} // namespace zephyr::geom::generator

/// @brief Хэш-функция для BaseEdge (для хранения в unordered_map)
template <>
struct std::hash<zephyr::geom::generator::BaseEdge> {
    size_t operator()(const zephyr::geom::generator::BaseEdge& pp) const noexcept {
        hash<const void*> hasher;
        size_t h1 = hasher(pp.ptr1());
        size_t h2 = hasher(pp.ptr2());
        return (h1 + h2) * (h1 + h2 + 1) / 2 + std::min(h1, h2);
    }
};
