#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/boundary.h>
#include <zephyr/geom/primitives/adjacent.h>

namespace zephyr::geom {

/// @brief Перечисление используется для выбора граней с определенным
/// направлением нормалей
enum class Direction : int {
    ANY = 0,  // Любое направление нормали
    X   = 1,
    Y   = 2,
    Z   = 3
};

/// @brief Базовый класс грани произвольной ячейки, хранит индексы
/// вершин (узлов), которые хранятся в массиве в самой ячейке.
/// BFace -- Base/Basic face
class BFace {
public:
    /// @brief Максимальное число вершин грани
    static const int max_size = 4;

    /// @brief Тип граничного условия
    Boundary boundary = Boundary::UNDEFINED;

    /// @brief Составной индекс смежной ячейки
    Adjacent adjacent = {};

    /// @brief Площадь грани
    double area = NAN;

    /// @brief Внешняя нормаль грани
    Vector3d normal;

    /// @brief Центр грани
    Vector3d center;

    /// @brief Список индексов вершин в массиве вершин ячейки,
    /// которой принадлежит грань
    std::array<int, max_size> vertices = {-1, -1, -1, -1};

    /// @brief Конструктор по умолчанию
    BFace() = default;

    /// @brief Является ли грань граничной?
    inline bool is_boundary() const {
        return boundary != Boundary::ORDINARY &&
               boundary != Boundary::PERIODIC &&
               boundary != Boundary::UNDEFINED;
    }

    /// @brief Является ли грань актуальной?
    inline bool is_actual() const {
	    return boundary != Boundary::UNDEFINED;
	}

	/// @return 'true', если грань не актуальна
	inline bool is_undefined() const {
	    return boundary == Boundary::UNDEFINED;
	}

	/// @brief Установить неопределенную грань
    inline void set_undefined() {
        boundary = Boundary::UNDEFINED;
    }

    /// @brief Число вершин грани
    /// Для граней двумерных ячеек равно двум,
    /// для граней трехмерных ячеек: 3 или 4.
    inline int size() const {
        static_assert(max_size == 4);
        return vertices[2] < 0 ? 2 : (vertices[3] < 0 ? 3 : 4);
    }

    /// @brief Вершина содержит указаный индекс?
    inline bool contain(int idx) const {
        static_assert(max_size == 4);
        return vertices[0] == idx || vertices[1] == idx ||
               vertices[2] == idx || vertices[3] == idx;
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает с выбраным направлением
    inline bool to_skip(Direction dir) const {
	    if (boundary == Boundary::UNDEFINED) {
            return true;
	    }
	    switch (dir) {
	        case Direction::ANY:
                return false;
	        case Direction::X:
                return std::abs(normal.x()) < 0.8;
	        case Direction::Y:
                return std::abs(normal.y()) < 0.8;
	        case Direction::Z:
                return std::abs(normal.z()) < 0.8;
            default:
                return false;
	    }
	}
};

} // namespace zephyr::geom