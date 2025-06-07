/// @file Файл содержит реализацию простой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для
/// разработчиков. Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/siblings.h>
#include <zephyr/mesh/amr2/balancing_slow.h>
#include <zephyr/mesh/amr2/balancing_fast.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr::mesh::amr2 {

/// @brief Функция балансировки флагов верхнего уровня
inline void balance_flags(AmrCells &cells, int max_level) {
    if (cells.empty()) return;

#if FAST_BALANCING
    if (cells.dim() < 3) {
        amr2::balance_flags_fast<2>(cells, max_level);
    } else {
        amr2::balance_flags_fast<3>(cells, max_level);
    }
#else
    if (cells.dim < 3) {
        amr2::balance_flags_slow<2>(cells, max_level);
    } else {
        amr2::balance_flags_slow<3>(cells, max_level);
    }
#endif

#if SCRUTINY
    amr2::check_flags(cells, max_level);
#endif
}

#ifdef ZEPHYR_MPI
/// @brief Функция балансировки флагов верхнего уровня блээээт
void balance_flags(AmrStorage &locals, AmrStorage& aliens, int max_level, EuMesh& mesh) {
    if (locals.empty())
        return;

    auto dim = locals[0].dim;

    if (dim < 3) {
        amr2::balance_flags_slow<2>(locals, aliens, max_level, mesh);
    } else {
        amr2::balance_flags_slow<3>(locals, aliens, max_level, mesh);
    }

#if SCRUTINY
    amr2::check_flags(locals, aliens, max_level);
#endif
}
#endif

} // namespace zephyr::mesh::amr2