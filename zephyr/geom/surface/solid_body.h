#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::mesh {
class EuCell; // forward declaration
}

namespace zephyr::geom {

class SolidBody {
public:
    /// @brief Базовый конструктор
    SolidBody();

    /// @brief Конструктор с указанием центра тела
    explicit SolidBody(const Vector3d& center);

    /// @brief Виртуальный деструктор
    virtual ~SolidBody() = default;

    /// @brief Положение центра тела
    Vector3d center() const;

    /// @brief Сдвинуть тело на вектор
    void move(const Vector3d& shift);

    /// @brief Поворот в локальной системе координат
    /// @param R Матрица поворота
    void rotation_local(const Matrix3d& R);

    /// @brief Поворот в глобальной системе координат
    /// @param R Матрица поворота
    void rotation_global(const Matrix3d& R);

    /// @brief Поворот относительно точки
    /// @param R Матрица поворота
    /// @param c Точка, вокруг которой выполняется вращение
    void rotation_relative(const Matrix3d& R, const Vector3d& c);


    /// @brief Точка находится внутри тела?
    bool inside(const Vector3d& v) const;

    /// @brief Объемная доля тела внутри ячейки сетки
    /// @param cell Ячейка сетки
    /// @param eps Точность определения объемной доли
    virtual double volume_fraction(const mesh::EuCell& cell, double eps) const;

    /// @brief Объем тела внутри ячейки
    /// @param cell Ячейка сетки
    /// @param eps Точность определения объемной доли
    virtual double volume_inside(const mesh::EuCell& cell, double eps) const;

protected:
    /// @brief Координаты точки в локальной системе координат,
    /// связанной с телом
    Vector3d in_local(const Vector3d& v);

    Vector3d m_center;   ///< Точка внутри тела/центр масс
    Matrix3d m_rotation; ///< Поворот относительно центра тела

    /// @brief Характеристическая функция области
    std::function<bool(const Vector3d& v)> m_inside;
};

} // namespace zephyr::geom
