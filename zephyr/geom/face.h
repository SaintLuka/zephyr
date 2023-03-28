#pragma once

#include <string>
#include <limits>
#include <array>
#include <algorithm>

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
///    index :    < locals.size  |    < decomposition(rank).locals.size
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

    /// @brief Центр грани
    template <int dim>
    inline Vector3d center(const Vertices& verts) const {
        if (dim < 3) {
            return (verts[vertices[0]] + verts[vertices[1]]) / 2.0;
        }
        else {
            return (verts[vertices[0]] + verts[vertices[1]] +
                    verts[vertices[2]] + verts[vertices[3]]) / 4.0;
        }
    }

    /// @brief Центр грани
    inline Vector3d center(const Vertices& verts, int dim) const {
        if (dim < 3) {
            return (verts[vertices[0]] + verts[vertices[1]]) / 2.0;
        }
        else {
            return (verts[vertices[0]] + verts[vertices[1]] +
                    verts[vertices[2]] + verts[vertices[3]]) / 4.0;
        }
    }
};

} // geom
} // zephyr