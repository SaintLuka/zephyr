// Базовые настройки адаптации. Сразу делаем include и using для всех классов,
// которые используются внутри namespace amr.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

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

#include <set>
#include <array>
#include <iomanip>
#include <numeric>

#include <zephyr/configuration.h>

#include <zephyr/mesh/side.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/utils/threads.h>

#ifdef ZEPHYR_MPI
#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/euler/tourism.h>
#endif

#ifdef CHECK_PERFORMANCE
#include <zephyr/utils/stopwatch.h>
#endif

namespace zephyr::mesh::amr {

using utils::threads;
using mesh::AmrCells;
using mesh::Children;
using geom::Vector3d;
using geom::Line;
using geom::Quad;
using geom::SqQuad;
using geom::Cube;
using geom::SqCube;
using geom::Boundary;

#ifdef ZEPHYR_MPI
using utils::mpi;
using mesh::Tourism;
#endif

#ifdef CHECK_PERFORMANCE
using utils::Stopwatch;

/// @brief Частота вывода сообщений о производительности
static const size_t check_frequency = 100;
#endif

/// @brief Локальные индексы ячеек, прилегающих к стороне.
template<int dim>
constexpr std::array<std::array<int, FpF(dim)>, FpC(dim)> get_children_by_side() {
    if constexpr (dim == 2) {
        return {
            std::array<int, 2>({0, 2}),
            std::array<int, 2>({1, 3}),
            std::array<int, 2>({0, 1}),
            std::array<int, 2>({2, 3})
        };
    }
    else {
        return {
            std::array<int, 4>({0, 2, 4, 6}),
            std::array<int, 4>({1, 3, 5, 7}),
            std::array<int, 4>({0, 1, 4, 5}),
            std::array<int, 4>({2, 3, 6, 7}),
            std::array<int, 4>({0, 1, 2, 3}),
            std::array<int, 4>({4, 5, 6, 7})
        };
    }
}

} // namespace zephyr::mesh::amr
