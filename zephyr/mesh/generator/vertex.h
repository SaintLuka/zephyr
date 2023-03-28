#pragma once

#include <memory>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr { namespace mesh { namespace generator {

using zephyr::geom::Vector3d;

class Block;
class Curve;

/// @brief Базисная вершина.
/// Через базисные вершины проходят границы области, на базисных вершинах
/// строятся базовые структурированные блоки.
class BaseVertex {
public:
    using Ptr = std::shared_ptr<BaseVertex>;
    using Ref = const std::shared_ptr<BaseVertex> &;

    /// @brief Конструктор вершины
    /// @param fixed Фиксированная вершина?
    BaseVertex(const Vector3d &v, bool fixed);

    /// @brief Конструктор вершины
    /// @param fixed Фиксированная вершина?
    static BaseVertex::Ptr create(const Vector3d &v, bool fixed);

    /// @brief Конструктор вершины
    /// @param x, y Координаты точки
    /// @param fixed Фиксированная вершина?
    static BaseVertex::Ptr create(double x, double y, bool fixed);

    /// @brief Является ли точка неподвижной?
    bool is_fixed() const;

    /// @brief Ссылка на содержимое в виде zephyr::Vector3d
    const Vector3d &v() const;

    /// @brief Неявное приведение к типу zephyr::Vector3d
    operator const Vector3d&() const {
        return m_v;
    }

    /// @brief Число смежных блоков
    /// @details Для базисных вершин большое значение играет число смежных
    /// блоков. Для базисной вершины на границе области нормальным является
    /// наличие двух смежных блоков (одного блока для угловой вершины),
    /// для базисной вершины внутри области нормальное число смежных блоков
    /// равно четырем. В этих случаях вершина считается регулярной,
    /// иначе - сингулярной.
    int degree() const;

    /// @brief Добавить смежный блок
    /// @details Производится проверка, что такой блок отсутствует
    void add_adjacent_block(Block *block);

    /// @brief Массив смежных блоков
    const std::vector<Block *> &adjacent_blocks() const;

private:

    /// @brief Положение
    Vector3d m_v;

    /// @brief Точка может быть подвижной или фиксированной (определяется
    /// только в конструкторе). В качестве фиксированных точек следует выбирать
    /// угловые точки области. Но при необходимости можно фиксировать
    /// произвольные точки на границе или даже внутри области.
    bool m_fixed;

    /// @brief Смежные блоки
    std::vector<Block *> m_adjacent_blocks;
};

/// @brief Внутренний тип вершины для класса Block
struct Vertex {
public:
    using Ptr = std::shared_ptr<Vertex>;
    using Ref = const std::shared_ptr<Vertex> &;

    /// @brief Конструктор
    explicit Vertex(const Vector3d &v);

    /// @brief Конструктор
    explicit Vertex(double x, double y);

    /// @brief Создание умного указателя
    static Vertex::Ptr create(const Vector3d &v);

    /// @brief Создание умного указателя
    static Vertex::Ptr create(double x, double y);

    /// @brief Число смежных вершин, больше единицы.
    /// @details Для вершины на границе области нормальным является наличие
    /// трех смежных вершин, для вершины внутри области нормальное число
    /// смежных вершин равно четырем. В этих случаях вершина считается регулярной,
    /// иначе - сингулярной. Сингулярные вершины встречаются только в углах
    /// структурированных блоков.
    int n_adjacent() const;

    /// @brief Смежные вершины
    const std::vector<Vertex *> &adjacent_vertices() const;

    /// @brief Установить соседние вершины. Вершины массива переставляются
    /// в соответствии с правилами обхода для смежных вершин. Первая вершина
    /// остается в начале.
    void set_adjacent_vertices(const std::vector<Vertex::Ptr> &vertices);

    /// @brief Указатель на границу
    Curve *boundary() const;

    /// @brief Установить границу
    void set_boundary(Curve *boundary);


    Vector3d v1;  ///< Основное положение
    Vector3d v2;  ///< Дополнительное положение

private:
    /// @brief Кривая, на которой лежит вершина.
    /// Если кривая отсутствует (nullptr), то вершина считается внутренней.
    Curve *m_boundary;

    /// @brief Список смежных вершин. Используется строгий порядок заполнения.
    /// Вершины обходятся против часовой стрелки внутри области. Таким образом,
    /// для граничной вершины крайние вершины внутри списка m_adjacent также
    /// должны быть граничными.
    std::vector<Vertex *> m_adjacent;
};

} // namespace generator
} // namespace mesh
} // namespace zephyr
