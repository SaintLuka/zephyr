#pragma once

#include <zephyr/geom/primitives/bface.h>
#include <zephyr/mesh/euler/amr_storage.h>

namespace zephyr::mesh {

class EuCell;
using zephyr::geom::BFace;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::geom::Direction;

/// @brief Обертка для типа geom::BFace, реализует интерфейс грани,
/// необходимый для работы. Также содержит несколько новых функций.
class EuFace {
public:

    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    EuFace(const EuCell &cell, geom::Side side);

    /// @brief Внешняя нормаль
    const Vector3d &normal() const;

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    EuCell neib() const;

    /// @brief Ссылка на данные соседней ячейки через грань
    /// @details При отсутствии соседа возвращает данные текущей ячейки.
    template <class T>
    inline const T& neib(const T&) const {
        return *reinterpret_cast<const T *>(neib_data());
    }

    /// @brief Ссылка на данные соседней ячейки через грань
    /// @details При отсутствии соседа возвращает данные текущей ячейки.
    template <class T>
    inline const T& neib(const VarExtra<T>& var) const {
        return *reinterpret_cast<const T *>(neib_data() + var.offset);
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
    inline const Vector3d& center() const { return m_face->center; }

    /// @brief Координата x центра грани
    inline double x() const { return center().x(); }

    /// @brief Координата y центра грани
    inline double y() const { return center().y(); }

    /// @brief Координата z центра грани
    inline double z() const { return center().z(); }

    /// @brief Количество вершин грани
    inline int size() const { return m_face->size(); }

    /// @brief Доступ к вершине грани
    const Vector3d& vs(int idx) const;


    inline const geom::Adjacent &adjacent() const {
        return m_face->adjacent;
    }

    EuFace &operator++() {
        do {
            m_face += 1;
        } while (m_face < m_end && m_face->to_skip(m_dir));
        return *this;
    }

    bool operator!=(const EuFace &iface) {
        return m_face != iface.m_face;
    }

    EuFace &operator*() {
        return *this;
    }

protected:
    friend class EuFaces;

    /// @brief Создать грань по ячейке
    /// @param cell Ячейка сетки
    /// @param self Указатель на геометрию текущей грани
    /// @param end Указатель на область за границей массива
    /// @param dir Выбирать грани с определенным направлением
    EuFace(const EuCell &cell, BFace* self, BFace* end,
           Direction dir = Direction::ANY);

    const Byte* neib_data() const;

    const EuCell &m_cell;  ///< Родительская ячейка
    BFace *m_face;         ///< Указатель на геометрию грани
    BFace *m_end;          ///< Указатель на грань за границей массива
    Direction m_dir;       ///< Направление
};


/// @brief Интерфейс для итерирования по граням ячейки
class EuFaces {
public:

    EuFaces(const EuCell& cell,
            Direction dir = Direction::ANY);

    EuFace begin() const;

    EuFace end() const;

private:
    const EuCell &m_cell;  ///< Родительская ячейка
    Direction m_dir;       ///< Направления граней
};

} // namespace zephyr::mesh