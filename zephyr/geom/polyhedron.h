#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/box.h>
#include <zephyr/geom/vector.h>
#include <zephyr/geom/line.h>
#include <zephyr/geom/primitives/cell_type.h>

namespace zephyr::geom {

/// @brief Представление многогранника
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

    /// @brief Пустой многогранник?
    inline bool empty() const { return vs.empty(); }

    /// @brief Число вершин
    int n_verts() const { return vs.size(); }

    /// @brief Число граней
    int n_faces() const { return fs.size(); }

    /// @brief Вершина по индексу
    const Vector3d& vertex(int idx) const {
        return vs[idx];
    }

    /// @brief Индексы вершин грани
    const std::vector<int>& face_indices(int idx) const {
        return fs[idx];
    }

    /// @brief Ограничивающий прямоугольник
    Box bbox() const;

    /// @brief Центр полигона (среднее вершин)
    Vector3d center() const;

    /// @brief Площадь грани
    double face_area(int idx) const;

    /// @brief Центр грани (среднее вершин)
    Vector3d face_center(int idx) const;

    /// @brief Внешняя нормаль к грани
    Vector3d face_normal(int idx) const;

    /// @brief Точка внутри многогранника?
    bool inside(const Vector3d& p) const;

    /// @brief Объем произвольного многогранника
    double volume() const;

    /// @brief Барицентр произвольного многогранника.
    /// @param vol Объем многогранника (если известен)
    /// @details На данный момент приближение -- просто ценр
    Vector3d centroid(double vol = 0.0) const;

    /// @brief Объем внутри многогранника на пересечении с характеристической
    /// функцией inside. Вычисляется приближенно с точностью ~ 1 / N
    /// @param inside Характеристическая функция, возвращает true, если точка p
    /// находится внутри области, иначе -- false.
    /// @param N Число пробных точек
    double clip_volume(const std::function<bool(const Vector3d& p)>& inside, int N = 10000) const;

    /// @brief Отсечь от полигона часть с помощью плоскости с внешней нормалью n,
    /// проходящей через точку p
    /// @return Отсеченный многогранник
    Polyhedron clip(const Vector3d& p, const Vector3d& n) const;

    /// @brief Объем многогранника, отсекаемого прямой с внешней нормалью n,
    /// проходящей через точку p.
    /// Равносильно вызову polyhedron.clip(p, n).volume(), но быстрее.
    double clip_volume(const Vector3d& p, const Vector3d& n) const;

    /// @brief Находит отсечение от многогранника с заданой объемной долей
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

protected:
    void build(const std::vector<Vector3d>& vertices,
               const std::vector<std::vector<int>>& face_indices);

    std::vector<Vector3d> vs;          ///< Массив вершин
    std::vector<std::vector<int>> fs;  ///< Индексы вершин граней
    std::vector<Vector3d> fcs;         ///< Центры граней
    std::vector<Vector3d> fns;         ///< Внешние нормали граней
    Vector3d m_center;                 ///< Центр полигона
};

std::ostream& operator<<(std::ostream& os, const Polyhedron& poly);

} // namespace zephyr::geom