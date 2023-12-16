#pragma once

#include <zephyr/geom/primitives/element.h>
#include <zephyr/geom/primitives/bvertices.h>
#include <zephyr/geom/primitives/bfaces.h>
#include <zephyr/geom/primitives/bnodes.h>

namespace zephyr::geom {

/// @class Обязательные данные ячейки сетки
class AmrCell : public Element {
public:
    // Геометрия ячейки

    int dim;       ///< Размерность ячейки
    bool adaptive;  ///< Адаптивная ячейка?
    bool linear;    ///< Линейная ячейка?
    Vector3d center;    ///< Барицентр ячейки
    double size;      ///< Линейный размер ячейки
    BVertices vertices;  ///< Вершины ячейки
    BFaces faces;     ///< Список граней

    /// @brief Индексы узлов в глобальном массиве вершин.
    /// Необязательное поле, используется только в некоторых алгоритмах,
    /// обычно заполнено неопределенными индексами (-1)
    BNodes nodes;

    // Данные AMR

    int b_idx;  ///< Индекс среди базовых ячеек
    int z_idx;  ///< Индекс ячейки на z-кривой
    int level;  ///< Уровень адаптации (0 для базовой)
    int flag;   ///< Желаемый флаг адаптации


    /// @brief Двумерная простая
    explicit AmrCell(const Quad &quad);

    /// @brief Двумерная криволинейная
    explicit AmrCell(const SqQuad &quad);

    /// @brief Трехмерная простая
    explicit AmrCell(const Cube &cube);

    /// @brief Трехмерная криволинейная ячейка
    explicit AmrCell(const SqCube &cube);

    /// @brief Двумерная полигональная ячейка. Не адаптивная ячейка, может
    /// представлять четырехугольник, но вершины упорядочены иначе.
    explicit AmrCell(const Polygon &poly);

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    double volume() const;

    /// @brief Скорпировать вершины в полигон (двумерные ячейки)
    /// Для нелинейных AMR-ячеек возвращает до 8 граней.
    PolygonS<8> polygon() const;

    /// @brief Диаметр вписаной окружности.
    /// @details Для AMR-ячейки представляет собой минимальное расстояние между
    /// противоположными гранями. Для полигона --- диаметр вписаной окружности
    /// для правильного многоугольника аналогичной площади.
    /// Величину удобно использовать совместно с условием Куранта.
    /// Для двумерных расчетов на прямоугольных сетках совпадает с минимальной
    /// стороной прямоугольной ячейки.
    double incircle_radius() const;

    /// @brief Оценка объемной доли, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @details Относительно быстрая функция, проверяет функцию inside только
    /// на узлах ячейки, позволяет быстро выяснить, содержит ли ячейка
    /// границу двух областей. Если ячейка внутри тела, то возвращает строго
    /// единицу 1.0, если снаружи -- строго ноль 0.0.
    double approx_vol_fraction(const std::function<double(const Vector3d &)> &inside) const;

    /// @brief Объемная доля, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @param n_points Число тестовых точек, для которых проверяется функция
    /// inside, погрешность определения объемной доли ~ 1/N.
    double volume_fraction(const std::function<double(const Vector3d &)> &inside, int n_points) const;

    /// @brief Помечает индексы в массиве nodes. Выставляет индексы актуальных
    /// вершин равными flag, для остальных вершин выставляет значение -1.
    /// Актуальными считаем вершины, на которые ссылаются грани.
    void mark_actual_nodes(int mark);

    /// @brief Вывести информацию о ячейке
    void print_info() const;

    /// @brief Вывести информацию о ячейке в виде python скрипта
    /// для визуализации
    void visualize(std::string filename) const;

    /// @brief Проверить базовую геометрию ячейки
    /// @return -1 для плохой ячейки
    int check_geometry() const;

    /// @brief Проверить ориентацию граней
    /// @return -1 для плохой ячейки
    int check_base_face_orientation() const;

    /// @brief Проверить порядок вершин
    /// @return -1 для плохой ячейки
    int check_base_vertices_order() const;

    /// @brief Проверить сложные грани
    /// @return -1 для плохой ячейки
    int check_complex_faces() const;

};

} // namespace zephyr::geom