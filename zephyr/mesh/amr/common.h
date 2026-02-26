// Базовые настройки адаптации. Сразу делаем include и using для всех классов,
// которые используются внутри namespace amr.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#define SCRUTINY 0           // Тщательная проверка вычислений
#define CHECK_PERFORMANCE 0  // Выводить производительность частей алгоритма
#define FAST_BALANCING 1     // Быстрая функция балансировки флагов

#if SCRUTINY
// Бросает исключение, если условие не выполняется
#define scrutiny_check(condition, message) if (!(condition)) throw std::runtime_error(message);
#else
// Опция SCRUTINY отключена, ничего не происходит
#define scrutiny_check(...)
#endif

#include <array>
#include <bitset>

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

inline index_t pack_children(index_t next, std::bitset<8> children) {
    z_assert(next < (1 << 23), "Too large next");
    return static_cast<index_t>((children.to_ulong() << 24) + next);
}

inline std::tuple<index_t, std::bitset<8>> unpack_children(index_t next) {
    return {static_cast<std::make_unsigned_t<index_t>>(next) % (1 << 24), std::bitset<8>(next >> 24) };
}

// coded_children - закодированные дочерние ячейки, формата bitset<8> [01010101]
// z_ch локальный индекс дочерней ячейки (внутри родительской)
inline int child_next(index_t coded_next, int z_ch) {
    auto[next, children] = unpack_children(coded_next);
    scrutiny_check(children[z_ch], "Has no such children");
    // Узнать номер nc внутри cset
    for (int i = 0; i < z_ch; ++i) {
        if (children[i]) ++next;
    }
    return next;
}

} // namespace zephyr::mesh::amr