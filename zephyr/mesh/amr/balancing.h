/// @file Файл содержит реализацию простой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для
/// разработчиков. Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_slow.h>
#include <zephyr/mesh/amr/balancing_fast.h>

namespace zephyr::mesh::amr {

/// @brief Функция балансировки флагов верхнего уровня блээээт
void balance_flags(AmrStorage &locals, int max_level) {
    if (locals.empty())
        return;

    AmrStorage aliens;

    auto dim = locals[0].dim;

#if FAST_BALANCING
    if (dim < 3) {
        amr::balance_flags_fast<2>(locals, max_level);
    } else {
        amr::balance_flags_fast<3>(locals, max_level);
    }
#else
    if (dim < 3) {
        amr::balance_flags_slow<2>(locals, aliens, max_level);
    } else {
        amr::balance_flags_slow<3>(locals, aliens, max_level);
    }
#endif

#if SCRUTINY
    amr::check_flags(locals, aliens, max_level);
#endif
}

#ifdef ZEPHYR_MPI
/// @brief Функция балансировки флагов верхнего уровня блээээт
void balance_flags(AmrStorage &locals, AmrStorage& aliens, int max_level, EuMesh& mesh) {
    if (locals.empty())
        return;

    auto dim = locals[0].dim;

    if (dim < 3) {
        amr::balance_flags_slow<2>(locals, aliens, max_level, mesh);
    } else {
        amr::balance_flags_slow<3>(locals, aliens, max_level, mesh);
    }

#if SCRUTINY
    amr::check_flags(locals, aliens, max_level);
#endif
}
#endif

} // namespace zephyr::mesh::amr