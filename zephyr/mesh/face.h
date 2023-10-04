#pragma once

#include <zephyr/geom/face.h>


namespace zephyr { namespace mesh {

class ICell;
using zephyr::geom::Face;
using zephyr::geom::FaceFlag;
using zephyr::geom::Vector3d;
using zephyr::geom::Direction;

/// @brief Обертка для типа geom::Face, реализует интерфейс грани,
/// необходимый для работы. Также содержит несколько новых функций.
class IFace {
public:
    /// @brief Создать грань по ячейке
    /// @param cell Ячейка сетки
    /// @param fid Индекс грани в ячейке
    /// @param dir Выбирать грани с определенным направлением
    IFace(const ICell &cell, int fid,
            Direction dir = Direction::ANY);

    //Face& face();

    const Face &geom() const;

    /// @brief Внешняя нормаль
    const Vector3d &normal() const;

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    ICell neib() const;

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    ICell neighbor() const;

    /// @brief Флаг граничных условий
    FaceFlag flag() const;

    /// @brief Является ли грань граничной
    bool is_boundary() const;

    /// @brief Установить флаг грани
    void set_boundary(FaceFlag flag);

    /// @brief Площадь (в 3D) или длина (в 2D) грани
    double area() const;

    /// @brief Центр грани
    Vector3d center() const;

    /// @brief Координата x центра грани
    inline double x() const { return center().x(); }

    /// @brief Координата y центра грани
    inline double y() const { return center().y(); }

    /// @brief Координата z центра грани
    inline double z() const { return center().z(); }

    IFace &operator++() {
        do {
            ++face_idx;
        } while (face_idx < geom::Faces::max_size && geom().to_skip(dir));
        return *this;
    }

    bool operator!=(const IFace &iface) {
        return face_idx != iface.face_idx;
    }

    IFace &operator*() {
        return *this;
    }

public:
    const ICell &m_cell;  ///< Родительская ячейка
    int face_idx;         ///< Индекс грани в ячейке
    Direction dir;        ///< Направление
};


/// @brief Интерфейс для итерирования по граням ячейки
class IFaces {
public:

    IFaces(const ICell& cell,
            Direction dir = Direction::ANY);

    IFace begin() const;

    IFace end() const;

private:
    const ICell &m_cell;  ///< Родительская ячейка
    Direction dir;        ///< Направления граней
};

} // namespace mesh
} // namespace zephyr