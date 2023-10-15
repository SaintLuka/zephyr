#pragma once

#include <string>
#include <algorithm>

#include <zephyr/geom/primitives/base.h>
#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Составной индекс смежности ячейки.
/// @details Короткое объяснение:
///               local neib           remote neib
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < aliens.size
struct Adjacent {

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    int rank;

    /// @brief Индекс смежной ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    int index;

    /// @brief Индекс смежной ячейки в массиве aliens.
    int ghost;

    /// @brief Конструктор по-умолчанию.
    Adjacent() : rank(0), ghost(0), index(0) { }

    /// @brief Оператор сравнения
    bool operator!=(const Adjacent& adj) const {
        return adj.rank != rank || adj.index != index || adj.ghost != ghost;
    }
};

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
template<int max_size>
class BaseFace {
public:

    /// @brief Площадь грани
    double area;

    /// @brief Список индексов вершин в массиве вершин ячейки,
    /// которой принадлежит грань
    std::array<short, max_size> vertices;

    /// @brief Внешняя нормаль грани
    Vector3d normal;

    /// @brief Тип граничного условия
    Boundary boundary;

    /// @brief Составной индекс смежной ячейки
    Adjacent adjacent;

    /// @brief Конструктор по умолчанию
    BaseFace() = default;

	/// @brief Является ли грань граничной?
	bool is_boundary() const {
        return boundary != Boundary::ORDINARY &&
               boundary != Boundary::PERIODIC &&
               boundary != Boundary::UNDEFINED;
    }

    /// @brief Является ли грань актуальной?
    bool is_actual() const {
	    return boundary != Boundary::UNDEFINED;
	}

	/// @return 'true', если грань не актуальна
	bool is_undefined() const {
	    return boundary == Boundary::UNDEFINED;
	}

	/// @brief Установить неопределенную грань
    void set_undefined() {
        area = -1.0;
        boundary = Boundary::UNDEFINED;
        vertices.fill(-1);
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает с выбраным направлением
    bool to_skip(Direction dir) const {
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