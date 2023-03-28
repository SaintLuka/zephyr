/// @file Файл содержит набор определений и реализацию некоторых простых функций,
/// используемых в алгоритмах адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/data/storage.h>
#include <zephyr/performance/timer/stopwatch.h>

#define SCRUTINY 0          ///< Тщательная проверка вычислений
#define CHECK_PERFORMANCE 0 ///< Выводить производительность частей алгоритма
#define FAST_BALANCING 1    ///< Быстрая функция балансировки флагов

#if SCRUTINY
/// @brief Бросает исключение, если условие не выполняется
#define scrutiny_check(condition, message) if (!(condition)) throw std::runtime_error(message);
#else
/// @brief Опция SCRUTINY отключена, ничего не происходит
#define scrutiny_check(...)
#endif


namespace zephyr { namespace mesh { namespace impl {

using namespace ::zephyr::data;
using ::zephyr::data::amrData;

/**
	\brief
        \~russian Число вершин на простой грани.
		\~english Vertices number per simple face.
		\~
*/
inline constexpr unsigned int VpF(unsigned int dim) { return dim < 3 ? 2 : 4; }

/**
	\brief
        \~russian Число вершин простой ячейки.
		\~english Vertices number per simple cell.
		\~
*/
inline constexpr unsigned int VpC(unsigned int dim) { return dim < 3 ? 4 : 8; }

/**
	\brief
        \~russian Число граней у простой ячейки.
		\~english Faces number per simple cell.
		\~
*/
inline constexpr unsigned int FpC(unsigned int dim) { return dim < 3 ? 4 : 6; }

/**
	\brief
        \~russian Количество дочерних ячеек.
		\~english Children number per cell.
		\~
*/
inline constexpr unsigned int CpC(unsigned int dim) { return dim < 3 ? 4 : 8; }

/**
	\brief
        \~russian Количество подграней.
		\~english Number of subfaces.
		\~
*/
inline constexpr unsigned int FpF(unsigned int dim) { return dim < 3 ? 2 : 4; }

/// @brief Required numbering of faces
static_assert(Side::::LEFT   == 0, "Left   != 0");
static_assert(Side::::RIGHT  == 1, "Right  != 1");
static_assert(Side::::BOTTOM == 2, "Bottom != 2");
static_assert(Side::::TOP    == 3, "Top    != 3");
static_assert(Side::::BACK   == 4, "Back   != 4");
static_assert(Side::::FRONT  == 5, "Front  != 5");

/**
	\brief
        \~russian Локальные индексы ячеек, прилегающих к стороне.
		\~english Local indices of children adjoined to a side.
		\~
*/
template<unsigned int dim>
inline constexpr std::array<std::array<int, VpF(dim)>, FpC(dim)> get_children_by_side();

/**
	\brief
        \~russian Локальные индексы ячеек, прилегающих к стороне (2D).
		\~english Local indices of children adjoined to a side (2D).
		\~
*/
template<>
inline constexpr std::array<std::array<int, 2>, 4> get_children_by_side<2>() {
    return {
            std::array<int, 2>({0, 2}),
            std::array<int, 2>({1, 3}),
            std::array<int, 2>({0, 1}),
            std::array<int, 2>({2, 3})
    };
}

/**
	\brief
        \~russian Локальные индексы ячеек, прилегающих к стороне (3D).
		\~english Local indices of children adjoined to a side (3D).
		\~
*/
template<>
inline constexpr std::array<std::array<int, 4>, 6> get_children_by_side<3>() {
    return {
            std::array<int, 4>({0, 2, 4, 6}),
            std::array<int, 4>({1, 3, 5, 7}),
            std::array<int, 4>({0, 1, 4, 5}),
            std::array<int, 4>({2, 3, 6, 7}),
            std::array<int, 4>({0, 1, 2, 3}),
            std::array<int, 4>({4, 5, 6, 7})
    };
}

/**
	\brief
        \~russian Индексы подграней на стороне в списке faces.
		\~english Indices of subfaces in faces.
		\~
*/
template<unsigned int dim>
inline std::array<int, FpF(dim)> subface_sides(int side);

/**
	\brief
        \~russian Индексы подграней на стороне в списке faces (2D).
		\~english Indices of subfaces in faces (2D).
		\~
*/
template<>
inline std::array<int, 2> subface_sides<2>(int side) {
    return {side, side + 6};
}

/**
	\brief
        \~russian Индексы подграней на стороне в списке faces (3D).
		\~english Indices of subfaces in faces (3D).
		\~
*/
template<>
inline std::array<int, 4> subface_sides<3>(int side) {
    return {side, side + 6, side + 12, side + 18};
}

/**
	\brief
        \~russian Найти центр простой грани.
		\~english Find center of a simple face.
		\~
	\param face
        \~russian Ссылка на грань.
        \~english Reference to the face.
    \param verts
        \~russian Соответствующий список вершин.
        \~english Reference to corresponding vertices.
*/
template <unsigned int dim>
inline Vector3d face_center(const _face_& face, const _vertices_& verts) {
    Vector3d center(0.0, 0.0, 0.0);
    for(unsigned int i = 0; i < VpF(dim); ++i) {
        center.x += verts.list[face.vertices[i]].x;
        center.y += verts.list[face.vertices[i]].y;
        center.z += verts.list[face.vertices[i]].z;
    }
    return center / VpF(dim);
}

/**
	\brief
        \~russian Копирует данные из одной ячейки в другую.
		\~english Copy data from the one cell to another.
		\~
	\param src
        \~russian Указатель на источник данных.
        \~english
    \param dest
        \~russian Указатель на область, куда надо записать данные.
        \~english
*/
inline void copy_data(Storage::iterator src, Storage::iterator dest) {
    _element_ elem = dest[element];
    _coords_ coords = dest[coords];
    _size_ size = dest[size];
    _vertices_ verts = dest[vertices];
    Faces faces = dest[faces];
    _amrData_ amr_data = dest[amrData];

    dest[data::item] = src[data::item];

    dest[element] = elem;
    dest[coords] = coords;
    dest[size] = size;
    dest[vertices] = verts;
    dest[faces] = faces;
    dest[amrData] = amr_data;
}

} // namespace impl
} // namespace mesh
} // namespace zephyr
