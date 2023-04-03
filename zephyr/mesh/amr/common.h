/// @file Файл содержит набор определений и реализацию некоторых простых функций,
/// используемых в алгоритмах адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/storage.h>

#define SCRUTINY 1          ///< Тщательная проверка вычислений
#define CHECK_PERFORMANCE 1 ///< Выводить производительность частей алгоритма
#define FAST_BALANCING 1    ///< Быстрая функция балансировки флагов

#if CHECK_PERFORMANCE
#include <zephyr/utils/stopwatch.h>
#endif

#if SCRUTINY
/// @brief Бросает исключение, если условие не выполняется
#define scrutiny_check(condition, message) if (!(condition)) throw std::runtime_error(message);
#else
/// @brief Опция SCRUTINY отключена, ничего не происходит
#define scrutiny_check(...)
#endif

namespace zephyr { namespace mesh { namespace amr {

using namespace zephyr::geom;

/// @brief Число вершин на простой грани.
inline constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число вершин простой ячейки.
inline constexpr int VpC(int dim) { return dim < 3 ? 4 : 8; }

/// @brief Число граней у простой ячейки.
inline constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
inline constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }

/// @brief Количество подграней.
inline constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Required numbering of faces
static_assert(Side::LEFT   == 0, "Left   != 0");
static_assert(Side::RIGHT  == 1, "Right  != 1");
static_assert(Side::BOTTOM == 2, "Bottom != 2");
static_assert(Side::TOP    == 3, "Top    != 3");
static_assert(Side::BACK   == 4, "Back   != 4");
static_assert(Side::FRONT  == 5, "Front  != 5");

/// @brief Локальные индексы ячеек, прилегающих к стороне.
template<int dim>
inline constexpr std::array<std::array<int, VpF(dim)>, FpC(dim)> get_children_by_side();

/// @brief Локальные индексы ячеек, прилегающих к стороне (2D).
template<>
inline constexpr std::array<std::array<int, 2>, 4> get_children_by_side<2>() {
    return {
            std::array<int, 2>({0, 2}),
            std::array<int, 2>({1, 3}),
            std::array<int, 2>({0, 1}),
            std::array<int, 2>({2, 3})
    };
}

/// @brief Локальные индексы ячеек, прилегающих к стороне (3D).
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

/// @brief Индексы подграней на стороне в списке faces.
template<int dim>
inline std::array<int, FpF(dim)> subface_sides(int side);

/// @brief  Индексы подграней на стороне в списке faces (2D).
template<>
inline std::array<int, 2> subface_sides<2>(int side) {
    return {side, side + 6};
}

/// @brief Индексы подграней на стороне в списке faces (3D).
template<>
inline std::array<int, 4> subface_sides<3>(int side) {
    return {side, side + 6, side + 12, side + 18};
}

/// @brief Найти центр простой грани.
/// @param face Ссылка на грань.
/// @param verts Соответствующий список вершин.
template <int dim>
inline Vector3d face_center(const Face& face, const Vertices& verts) {
    Vector3d center(0.0, 0.0, 0.0);
    for(int i = 0; i < VpF(dim); ++i) {
        center += verts[face.vertices[i]];
    }
    return center / VpF(dim);
}

/// @brief Копирует данные из одной ячейки в другую.
/// @param src Указатель на источник данных.
/// @param dest Указатель на область, куда надо записать данные.
inline void copy_data(Storage::Item src, Storage::Item dest) {
    std::memcpy(dest.geom_ptr(), src.geom_ptr(), sizeof(Cell));
    std::memcpy(dest.data_ptr(), src.data_ptr(), dest.datasize());
}

} // namespace amr
} // namespace mesh
} // namespace zephyr
