#pragma once

#include <string>

#include <zephyr/geom/base.h>
#include <zephyr/geom/vertices.h>

namespace zephyr { namespace geom {

namespace secret_namespace {

typedef short secretFaceFlagType;

enum secretFaceFlag : secretFaceFlagType {
    UNDEFINED = 0,  ///< Не определено.
    ORDINARY  = 1,  ///< Обычное соседство (внутри сетки).
    WALL      = 2,  ///< Непроницаемая стенка.
    ZOE       = 3,  ///< Отображение данных запрашиваемой ячейки (простой снос).
    PERIODIC  = 4,  ///< Периодичность.
    FUNCTION  = 5,  ///< Кастомная
};
}

/// @brief Флаг граничных условий
/// @details Работает как enum class, то есть обращение к элементам
/// перечисления FaceFlag::UNDEFINED, FaceFlag::ORDINARY и т.д.
/// Однако, отсутствует необходимость использовать static_cast, элементы
/// перечисления автоматически конвертируются в FaceFlagType
using FaceFlag = secret_namespace::secretFaceFlag;

/// @brief String representation of boundary condition
inline std::string boundary_to_string(FaceFlag type) {
    switch (type) {
        case FaceFlag::ORDINARY:
            return "ordinary";
        case FaceFlag::WALL:
            return "wall";
        case FaceFlag::ZOE:
            return "zoe";
        case FaceFlag::PERIODIC:
            return "periodic";
        case FaceFlag::FUNCTION:
            return "function";
        default:
            return "undefined";
    }
}

inline FaceFlag boundary_from_string(std::string flag) {
    std::transform(flag.begin(), flag.end(), flag.begin(), ::tolower);

    if (flag == "ordinary") {
        return FaceFlag::ORDINARY;
    }
    else if (flag == "wall") {
        return FaceFlag::WALL;
    }
    else if (flag == "zoe") {
        return FaceFlag::ZOE;
    }
    else if (flag == "periodic") {
        return FaceFlag::PERIODIC;
    }
    else if (flag == "function") {
        return FaceFlag::FUNCTION;
    }
    else if ("undefined") {
        return FaceFlag::UNDEFINED;
    }
    else {
        throw std::runtime_error("Unknown boundary flag '" + flag + "'");
    }
}

/// @brief Составной индекс смежности ячейки.
/// @details Короткое объяснение:
///               local neib           remote neib
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < aliens.size
struct Adjacent {

    /// @brief Номер датасета, в котором находится смежная ячейка
    /// при распределенном расчете.
    int rank;

    /// @brief Индекс смежной ячейки в локальном датасете, где находятся
    /// копии приграничных ячеек при распределенном счете.
    int ghost;

    /// @brief Индекс смежной ячейки в основном локальном датасете
    int index;

    /// @brief Конструктор по-умолчанию.
    Adjacent() : rank{0}, ghost{0}, index{0} { }

    /// @brief Оператор сравнения
    bool operator!=(const Adjacent& adj) const {
        return adj.index != index || adj.ghost != ghost || adj.rank != rank;
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

struct Face;

/// @brief Найти центр простой грани.
/// @param face Ссылка на грань.
/// @param verts Соответствующий список вершин.
template <int dim>
static Vector3d face_center(const Face& face, const Vertices& verts);

/// @brief Грань ячейки
struct Face {

    /// @brief Площадь грани
    double area;

    /// @brief Список индексов вершин в массиве Vertices
    std::array<short, 4> vertices;

    /// @brief Внешняя нормаль грани
    Vector3d normal;

    /// @brief Тип граничного условия
    FaceFlag boundary;

    /// @brief Составной индекс смежной ячейки
    Adjacent adjacent;

    /// @brief Конструктор по умолчанию
    Face() = default;

	/// @brief Является ли грань граничной
	bool is_boundary() const {
        return boundary != FaceFlag::ORDINARY &&
               boundary != FaceFlag::PERIODIC &&
               boundary != FaceFlag::UNDEFINED;
    }

    /// @brief Является ли грань актуальной.
    bool is_actual() const {
	    return boundary != FaceFlag::UNDEFINED;
	}

	/// @return 'true', если грань не актуальна
	bool is_undefined() const {
	    return boundary == FaceFlag::UNDEFINED;
	}

	/// @brief Установить неопределенную грань
    void set_undefined() {
        area = -1.0;
        boundary = FaceFlag::UNDEFINED;
        vertices = { -1, -1, -1, -1 };
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает по направлению
    bool to_skip(Direction dir) const {
	    if (boundary == FaceFlag::UNDEFINED) {
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

    /// @brief Центр грани
    template <int dim>
    inline Vector3d center(const Vertices& verts) const {
        if (dim < 3) {
            return 0.5 * (verts[vertices[0]] + verts[vertices[1]]);
        } else {
            return 0.25 * (verts[vertices[0]] + verts[vertices[1]] +
                           verts[vertices[2]] + verts[vertices[3]]);
        }
    }

    /// @brief Центр грани
    inline Vector3d center(const Vertices& verts, int dim) const {
	    return dim < 3 ? center<2>(verts) : center<3>(verts);
    }
};

/// @brief Центр грани
template <int dim>
inline Vector3d face_center(const Face& face, const Vertices& verts) {
    if (dim < 3) {
        return 0.5 * (verts[face.vertices[0]] + verts[face.vertices[1]]);
    } else {
        return 0.25 * (verts[face.vertices[0]] + verts[face.vertices[1]] +
                       verts[face.vertices[2]] + verts[face.vertices[3]]);
    }
}

} // geom
} // zephyr