// Базовые настройки адаптации. Сразу делаем include и using для всех классов,
// которые используются внутри namespace amr.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#define SCRUTINY 1           // Тщательная проверка вычислений
#define CHECK_PERFORMANCE 1  // Выводить производительность частей алгоритма
#define FAST_BALANCING 1     // Быстрая функция балансировки флагов

#if SCRUTINY
// Бросает исключение, если условие не выполняется
#define scrutiny_check(condition, message) if (!(condition)) throw std::runtime_error(message);
#else
// Опция SCRUTINY отключена, ничего не происходит
#define scrutiny_check(...)
#endif

#include <array>

#include <zephyr/configuration.h>

#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/mesh/side.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/tourism.h>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/mpi.h>

#ifdef CHECK_PERFORMANCE
#include <zephyr/utils/stopwatch.h>
#endif

namespace zephyr::mesh::amr {

using geom::Vector3d;
using geom::Line;
using geom::Quad;
using geom::SqQuad;
using geom::Cube;
using geom::SqCube;
using geom::Boundary;
using mesh::AmrCells;
using mesh::Children;

using utils::threads;
using utils::mpi;

#ifdef ZEPHYR_MPI
using mesh::Tourism;
#endif

#ifdef CHECK_PERFORMANCE
using utils::Stopwatch;

// Частота вывода сообщений о производительности
static constexpr size_t check_frequency = 100;
#endif

// Локальные индексы ячеек, прилегающих к стороне.
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