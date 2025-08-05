// Функция balancing верхнего уровня, выполняет балансировку флагов адаптации.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_slow.h>
#include <zephyr/mesh/amr/balancing_fast.h>

namespace zephyr::mesh::amr {

/// @brief Функция балансировки флагов верхнего уровня
inline void balance_flags(AmrCells &cells, int max_level) {
    if (cells.empty()) return;

#if FAST_BALANCING
    if (cells.dim() < 3) {
        amr::balance_flags_fast<2>(cells, max_level);
    } else {
        amr::balance_flags_fast<3>(cells, max_level);
    }
#else
    if (cells.dim < 3) {
        amr::balance_flags_slow<2>(cells, max_level);
    } else {
        amr::balance_flags_slow<3>(cells, max_level);
    }
#endif

#if SCRUTINY
    amr::check_flags(cells, max_level);
#endif
}

#ifdef ZEPHYR_MPI

/// @brief Функция балансировки флагов верхнего уровня
inline void balance_flags(Tourism& tourism, AmrCells &locals, AmrCells& aliens, int max_level) {
    if (locals.empty()) {
        amr::balance_flags_slow<0>(tourism, locals, aliens, max_level);
    }
    else {
        if (locals.dim() < 3) {
            amr::balance_flags_slow<2>(tourism, locals, aliens, max_level);
        } else {
            amr::balance_flags_slow<3>(tourism, locals, aliens, max_level);
        }
    }

#if SCRUTINY
    amr::check_flags(locals, aliens, max_level);
#endif
}

#endif

} // namespace zephyr::mesh::amr