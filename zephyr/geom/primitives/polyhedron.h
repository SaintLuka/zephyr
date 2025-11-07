#pragma once

#include <vector>
#include <iostream>
#include <functional>

#include <zephyr/geom/box.h>
#include <zephyr/geom/vector.h>
#include <zephyr/geom/cell_type.h>

namespace zephyr::geom {

/// @addtogroup geom-primitives
/// @{

/// @brief Многогранник общего вида.
class Polyhedron {
public:
    /// @brief Пустой многогранник (заглушка)
    Polyhedron() = default;

    /// @brief Основной конструктор многогранника
    /// @param vertices Список вершин
    /// @param face_indices Индексы вершин
    Polyhedron(const std::vector<Vector3d>& vertices,
               const std::vector<std::vector<int>>& face_indices);

    /// @brief Конструктор многогранника
    /// @param ctype Тип многогранника, один из VTK типов
    ///     CellType::TETRA        // Тетраэдр
    ///     CellType::PYRAMID,     // Пирамида с четырехугольным основанием
    ///     CellType::WEDGE,       // Призма с треугольным основанием
    ///     CellType::HEXAHEDRON,  // Топологический куб
    /// @param vertices Вершины должны быть упорядочены как соответствующие
    /// типы ячеек в VTK формате
    Polyhedron(CellType ctype, const std::vector<Vector3d>& vertices);

    /// @brief Создать пустой многогранник
    static Polyhedron Empty() { return Polyhedron(); }

    /// @brief Пустой многогранник?
    inline bool empty() const { return verts.empty(); }

    /// @brief Число вершин
    int n_verts() const { return verts.size(); }

    /// @brief Число граней
    int n_faces() const { return faces.size(); }

    /// @brief Вершина по индексу
    const Vector3d& vertex(int idx) const { return verts[idx]; }

    /// @brief Индексы вершин грани
    const std::vector<int>& face_indices(int idx) const { return faces[idx]; }

    /// @brief Ограничивающий прямоугольник
    Box bbox() const;

    /// @brief Центр полигона (среднее вершин)
    Vector3d center() const;

    /// @brief Радиус описаной окружности.
    /// На самом деле не совсем, возвращает максимальное расстояние от центра
    /// многогранника до его вершин
    double excircle_radius() const;

    /// @brief Площадь грани
    double face_area(int idx) const;

    /// @brief Центр грани (среднее вершин)
    Vector3d face_center(int idx) const;

    /// @brief Внешняя нормаль к грани
    Vector3d face_normal(int idx) const;

    /// @brief Внешняя нормаль к грани, умноженная на площадь
    Vector3d face_area_n(int idx) const;

    /// @brief Точка внутри многогранника?
    bool inside(const Vector3d& p) const;

    /// @brief Объем произвольного многогранника
    double volume() const;

    /// @brief Барицентр произвольного многогранника.
    /// @param vol Объем многогранника (если известен)
    /// @details На данный момент приближение -- просто ценр
    Vector3d centroid(double vol = 0.0) const;

    /// @brief Возвращает 'true', если одна из граней имеет более
    /// max_vertices вершин
    bool need_simplify(int max_vertices) const;

    /// @brief Разбить грани таким образом, чтобы каждая грань содержала
    /// не более max_vertices вершин. При этом вершины не добавляются.
    void simplify_faces(int max_vertices);

    /// @brief Объем внутри многогранника на пересечении с характеристической
    /// функцией inside. Вычисляется приближенно с точностью ~ 1 / N
    /// @param inside Характеристическая функция, возвращает true, если точка p
    /// находится внутри области, иначе -- false.
    /// @param N Число пробных точек
    double clip_volume(const std::function<bool(const Vector3d& p)>& inside,
                       int N = 10000) const;

    /// @brief Отсечь от полигона часть с помощью плоскости с внешней нормалью n,
    /// проходящей через точку p
    /// @return Отсеченный многогранник
    Polyhedron clip(const Vector3d& p, const Vector3d& n) const;

    /// @brief Объем многогранника, отсекаемого прямой с внешней нормалью n,
    /// проходящей через точку p.
    /// Равносильно вызову polyhedron.clip(p, n).volume(), но быстрее.
    double clip_volume(const Vector3d& p, const Vector3d& n) const;

    /// @brief Находит отсечение от многогранника с заданной объемной долей
    /// @param n Внешняя нормаль плоскости
    /// @param alpha Объемная доля
    /// @return Точка плоскости
    Vector3d find_section(const Vector3d& n, double alpha) const;

    /// @brief Посчитать объемную долю, которая отсекается от ячейки некоторым
    /// телом, точки которого определяются характеристической функцией inside
    /// @param inside Характеристическая функция: true, если точка находится
    /// внутри тела, иначе false
    /// @param n_points Число тестовых точек, погрешность ~ 1/N.
    double volume_fraction(const std::function<bool(const Vector3d&)>& inside,
                           int n_points = 10000) const;

    /// @brief Сдвинуть вершины на указанный вектор
    void move(const Vector3d& shift);

    /// @brief Проверяет многоугольник, выдает код неисправности.
    /// Если 0, то всё в порядке.
    int checkout() const;

    /// @brief Вывести в поток
    void print(std::ostream& os) const;

    /// @brief Вывести в стандартный поток
    void print() const { print(std::cout); }

    // ========================================================================
    //                          Несколько пресетов
    // ========================================================================

    /// @brief Единичный куб с центром в начале координат
    static Polyhedron Cube() { return Polyhedron::Cuboid(1.0, 1.0, 1.0); }

    /// @brief Прямоугольный параллелепипед со сторонами [a, b, c] и центром
    /// в начале координат
    static Polyhedron Cuboid(double a, double b, double c);

    /// @brief Четырехугольная пирамида, вписана в единичную сферу
    static Polyhedron Pyramid();

    /// @brief Треугольная призма, вписана в единичную сферу
    static Polyhedron Wedge();

    /// @brief Правильный тетраэдр, вписан в единичную сферу
    static Polyhedron Tetrahedron();

    /// @brief Правильный октаэдр, вписан в единичную сферу
    static Polyhedron Octahedron();

    /// @brief Правильный додекаэдр, вписан в единичную сферу
    static Polyhedron Dodecahedron();

    /// @brief Правильный икосаэдр, вписан в единичную сферу
    static Polyhedron Icosahedron();

    /// @brief Усеченный куб, вписан в единичный куб
    static Polyhedron TruncatedCube();

protected:
    void build(const std::vector<Vector3d>& vertices,
               const std::vector<std::vector<int>>& face_indices);


    /// @brief Упрощает грань, возвращает массивы индексов, на которых можно
    /// построить новые грани
    std::vector<std::vector<int>> simplified_faces(
            int face_idx, int max_vertices) const;

    /// @brief Заменить вершины грани
    void replace_face(int face_idx, const std::vector<int>& vs);

    /// @brief Добавить грань на существующих вершинах
    void add_face(const std::vector<int>& vs);

    // Базовые поля
    std::vector<Vector3d> verts;          ///< Массив вершин
    std::vector<std::vector<int>> faces;  ///< Индексы вершин граней

    // Вычисляемые
    Vector3d m_center;                    ///< Центр многогранника
    std::vector<Vector3d> faces_c;        ///< Центры граней
    std::vector<Vector3d> faces_s;        ///< Нормаль с площадью
};

/// @brief Вывод многогранника
static std::ostream& operator<<(std::ostream& os, const Polyhedron& poly) {
    poly.print(os);
    return os;
}

/// @}

} // namespace zephyr::geom