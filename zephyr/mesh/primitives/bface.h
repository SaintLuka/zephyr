#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/mesh/primitives/adjacent.h>

namespace zephyr::mesh {

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
    using Vector3d = zephyr::geom::Vector3d;
    using Boundary = zephyr::geom::Boundary;
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
        adjacent.rank  = -1;
        adjacent.index = -1;
        adjacent.alien = -1;
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

	// Расширение массива индексов для хранения полигональных граней,
	// Полная дичь, мне было скучно. В последний индекс vertices[3]
	// можно сунуть 6 индексов вместо одного.

    // У полигональной грани первые два бита индекса vertices[3] равны 0x10...
    // Для неопределенной вершины (vertices[3] == -1) все биты равны единице,
    // Для обычной определенной вершины первые два бита равны нулю.

    /// @brief Установить неопределенную полигональную грань
	inline void set_polygonal() {
        vertices[3] = 0b10111111111111111111111111111111;
    }

	/// @brief Полигональная грань? Проверяет первые биты
    inline bool poly_simple() const {
        return (vertices[3] & 0b11000000000000000000000000000000) !=
                              0b10000000000000000000000000000000;
    }

	/// @brief Число вершин полигональной грани
	/// Максимальное количество равно девяти.
	int poly_size() const {
        if (poly_simple()) { return size(); }
        for (int i = 2; i < 9; ++i) {
            if (get_poly_vertex(i) < 0) {
                return i;
            }
        }
        return 9;
    }

    /// @brief Получить индекс вершины
    int get_poly_vertex(int i) const {
        if (poly_simple() || i < 3) { return vertices[i]; }

        // Единички на нужном месте
        int mask = 0b11111 << (5 * (i - 3));
        int res = (vertices[3] & mask) >> (5 * (i - 3));
        return res < 31 ? res : -1;
    }

    /// @brief Установить индекс вершины, v_idx --- от 0 до 30
    /// фактически, максимальное число вершин BVertices::max_count = 24.
    void set_poly_vertex(int i, int v_idx) {
        if (poly_simple() || i < 3) {
            vertices[i] = v_idx;
            return;
        }

        // Зануляет значения на нужном месте
        vertices[3] &= ((0b11111 << (5 * (i - 3))) ^ 0xFFFFFFFF);

        // Помещает value на нужное место
        int shifted_value = (0b11111 & v_idx) << (5 * (i - 3));

        // По сути ставит value на обуленное место
        vertices[3] += shifted_value;
    }
};

} // namespace zephyr::mesh