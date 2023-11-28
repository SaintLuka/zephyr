#pragma once

#include <zephyr/mesh/storage.h>
#include <zephyr/geom/primitives/bface.h>


namespace zephyr::mesh {

class LaCell;
using zephyr::geom::BFace;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::geom::Direction;

/// @brief Обертка для типа geom::BFace, реализует интерфейс грани,
/// необходимый для работы. Также содержит несколько новых функций.
class LaFace {
public:

    /// @brief Внешняя нормаль
    const Vector3d &normal() const;

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    LaCell neib() const;

    /// @brief Ссылка на данные соседней ячейки через грань
    /// @details При отсутствии соседа возвращает данные текущей ячейки.
    template <class T>
    const T& neib(const T&) const {
        return *reinterpret_cast<const T *>(neib_data());
    }

    /// @brief Флаг граничных условий
    Boundary flag() const;

    /// @brief Является ли грань граничной
    bool is_boundary() const;

    /// @brief Установить флаг грани
    void set_boundary(Boundary flag);

    /// @brief Площадь (в 3D) или длина (в 2D) грани
    double area() const;

    /// @brief Центр грани
    Vector3d center() const;

    /// @brief Координата x центра грани
    inline double x() const { return m_face->center.x(); }

    /// @brief Координата y центра грани
    inline double y() const { return m_face->center.y(); }

    /// @brief Координата z центра грани
    inline double z() const { return m_face->center.z(); }

    /// @brief Доступ к вершине грани
    const Vector3d& vs(int idx) const;


    inline const geom::Adjacent &adjacent() const {
        return m_face->adjacent;
    }

    LaFace &operator++() {
        do {
            m_face += 1;
        } while (m_face < m_end && m_face->to_skip(m_dir));
        return *this;
    }

    bool operator!=(const LaFace &iface) {
        return m_face != iface.m_face;
    }

    LaFace &operator*() {
        return *this;
    }

protected:
    friend class LaFaces;

    /// @brief Создать грань по ячейке
    /// @param cell Ячейка сетки
    /// @param self Указатель на геометрию текущей грани
    /// @param end Указатель на область за границей массива
    /// @param dir Выбирать грани с определенным направлением
    LaFace(const LaCell &cell, BFace* self, BFace* end,
           Direction dir = Direction::ANY);

    const Byte* neib_data() const;

    const LaCell &m_cell;  ///< Родительская ячейка
    BFace *m_face;         ///< Указатель на геометрию грани
    BFace *m_end;          ///< Указатель на грань за границей массива
    Direction m_dir;       ///< Направление
};


/// @brief Интерфейс для итерирования по граням ячейки
class LaFaces {
public:

    LaFaces(const LaCell& cell,
            Direction dir = Direction::ANY);

    LaFace begin() const;

    LaFace end() const;

private:
    const LaCell &m_cell;  ///< Родительская ячейка
    Direction     m_dir;   ///< Направления граней
};

} // namespace zephyr::mesh