// Функция apply верхнего уровня, производит адаптацию сетки в соответствии
// с флагами адаптации. Предполагается, что флаги адаптации сбалансированы.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/statistics.h>
#include <zephyr/mesh/amr/setup_positions.h>
#include <zephyr/mesh/amr/setup_geometry.h>
#include <zephyr/mesh/euler/tourism.h>

namespace zephyr::mesh::amr {

/// @brief Функция выполняет непосредственную адаптацию ячеек в хранилище в
/// соответствии с флагами адаптации. Предполагается, что флаги адаптации
/// сбалансированы.
/// @param locals Ссылка на хранилище ячеек
/// @param op Осуществляет распределение данных
/// @details Детали алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т. д.
/// Этап 2. Создание обменных списков. Подразумевается, что новые ячейки из
/// конца хранилища будут перемещены на места удаленных ячеек.
/// Этап 3. Определение позиций для создания новых ячеек, а также финальные
/// позиции всех ячеек, включая новые. Индексы выставляются в поле NEXT.
/// Новые ячейки будут созданы за границами исходного хранилища.
/// В исходном хранилище NEXT указывает:
///   - флаг = 0: NEXT это финальное положение ячейки в хранилище.
///   - флаг < 0: NEXT это положение родительской (за границей исходного).
///   - флаг > 0: NEXT это положение первой дочерней (за границей исходного).
/// В расширенной части хранилища NEXT всегда указывает на финальное
/// положение ячейки. Таким образом, NEXT позволяет точно определить куда
/// в результате будут перемещены все ячейки, в том числе новые.
/// Этап 4. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища. Все связи выставляются точно,
/// исходя из финальных позиций всех ячеек (поле NEXT).
/// Этап 5. На начале этапа все ячейки правильно связаны с указанием финальных
/// индексов, но внутри хранилища часть старых ячеек (не листовых) являются
/// неопределенными. При этом за пределами исходного хранилища созданы новые
/// ячейки. Необходимо выполнить цикл по обменным спискам и сделать перемещение
/// ячеек из конца хранилища на места неопределенных ячеек.
template<int dim>
void apply_impl(AmrCells &locals, const Distributor& op) {
    static Stopwatch count_timer;
    static Stopwatch swap_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch remove_timer;

    /// Этап 1. Сбор статистики
    count_timer.resume();
    const Statistics count(locals.flag, locals.dim());
    //count.print();
    count_timer.stop();

    // Нечего адаптировать
    if (count.n_cells_large <= count.n_cells) return;
    scrutiny_check(count.n_cells_large != count.n_cells, "Strange things");

    // Этап 2. Определение обменных списков
    swap_timer.resume();
    SwapLists swap_list(count, locals.flag);
    swap_timer.stop();

    /// Этап 3. Распределяем места для новых ячеек
    positions_timer.resume();
    locals.resize_amr(count.n_cells_large);
    setup_positions<dim>(locals, count, swap_list);
    positions_timer.stop();

    /// Этап 4. Создание новых ячеек
    geometry_timer.resume();
    setup_geometry<dim>(locals, count, op);
    geometry_timer.stop();

    /// Этап 5. Перемещение готовых ячеек
    remove_timer.resume();
    swap_list.move_elements(locals);
    locals.resize_amr(count.n_cells_short);
    remove_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Statistics:       " << std::setw(9) << count_timer.milliseconds() << " ms\n";
        mpi::cout << "    Create SwapList:  " << std::setw(9) << swap_timer.milliseconds() << " ms\n";
        mpi::cout << "    Setup Positions:  " << std::setw(9) << positions_timer.milliseconds() << " ms\n";
        mpi::cout << "    Setup Geometry:   " << std::setw(9) << geometry_timer.milliseconds() << " ms\n";
        mpi::cout << "    Remove undefined: " << std::setw(9) << remove_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Автоматический выбор размерности
inline void apply(AmrCells &cells, const Distributor& op) {
    if (cells.empty()) return;

    if (cells.dim() < 3) {
        amr::apply_impl<2>(cells, op);
    }
    else {
        amr::apply_impl<3>(cells, op);
    }
}

#ifdef ZEPHYR_MPI
/// @brief Многопроцессорная версия функции, выполняет непосредственную
/// адаптацию ячеек в хранилище в соответствии с флагами адаптации.
/// Предполагается, что флаги адаптации сбалансированы.
/// @param tourism Обменные слои
/// @param locals Ссылка на хранилище ячеек
/// @param op Осуществляет распределение данных
/// @details Алгоритм практически совпадает с однопроцессорной версией,
/// поскольку не допускается огрубление ячеек, у которых сиблинги находятся
/// на различных процессах. Детали алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т. д.
/// Этап 2. Создание обменных списков. Подразумевается, что новые ячейки из
/// конца хранилища будут перемещены на места удаленных ячеек.
/// Этап 3а. Выполняется для locals ячеек.
/// Определение позиций для создания новых ячеек, а также финальные
/// позиции всех ячеек, включая новые. Индексы выставляются в поле NEXT.
/// Новые ячейки будут созданы за границами исходного хранилища.
/// В исходном хранилище NEXT указывает:
///   - флаг = 0: NEXT это финальное положение ячейки в хранилище.
///   - флаг < 0: NEXT это положение родительской (за границей исходного).
///   - флаг > 0: NEXT это положение первой дочерней (за границей исходного).
/// В расширенной части хранилища NEXT всегда указывает на финальное
/// положение ячейки. Таким образом, NEXT позволяет точно определить куда
/// в результате будут перемещены все ячейки, в том числе новые.
/// Этап 3б. Выполняется для aliens ячеек. Также проставляет параметр NEXT,
/// но только для border и aliens ячеек. Правил индексации NEXT отличаются.
/// Используется другой алгоритм формирования списка border. Ячейки не
/// перемещаются из конца на неопределенные места, а располагаются на тех же
/// местах последовательно. Для split ячейки её дети (которые принадлежат
/// border границе), будут располагаться все вместе в итоговом border на
/// том же месте, где была ячейка (ну разве что ячейку сдвинут).
/// Итого, в aliens массивах NEXT указывает:
///   - флаг = 0: NEXT это финальное положение ячейки в aliens-хранилище.
///   - флаг < 0: NEXT это финальное положение ячейки в aliens-хранилище.
///   - флаг > 0: NEXT закодированное положение первой дочерней ячейки в
///                    aliens-хранилище + закодированные дети.///
/// Этап 4. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища. Все связи выставляются точно,
/// исходя из финальных позиций всех ячеек (поле NEXT). А также с правильными
/// ссылками на финальную версию aliens.
/// Этап 5. На начале этапа все ячейки правильно связаны с указанием финальных
/// индексов, но внутри хранилища часть старых ячеек (не листовых) являются
/// неопределенными. При этом за пределами исходного хранилища созданы новые
/// ячейки. Необходимо выполнить цикл по обменным спискам и сделать перемещение
/// ячеек из конца хранилища на места неопределенных ячеек.
/// В конце этапа все хранилища масштабируются под финальные размеры.
/// Этап 6. Обменные слои все сделаны корректно, необходимо только запаковать
/// и отправить геометрию из border в aliens.
/// Этап 7. Проставить индексы adjacent.index для alien-ячеек.
template<int dim>
void apply_impl(AmrCells &locals, const Distributor& op, Tourism& tourism) {
    static Stopwatch count_timer;
    static Stopwatch swap_timer;
    static Stopwatch positions_timer1;
    static Stopwatch positions_timer2;
    static Stopwatch geometry_timer;
    static Stopwatch remove_timer;
    static Stopwatch exchange_timer;
    static Stopwatch restore_timer;

    // Этап 1. Сбор статистики
    count_timer.resume();
    const Statistics count(locals.flag, locals.dim());
    count_timer.stop();

    // Этап 2. Определение обменных списков
    swap_timer.resume();
    const SwapLists swap_list(count, locals.flag);
    swap_timer.stop();

    /// Этап 3а. Распределяем места для новых ячеек
    positions_timer1.resume();
    locals.resize_amr(count.n_cells_large);
    setup_positions<dim>(locals, count, swap_list);
    positions_timer1.stop();

    // Этап 3б. Сделать setup_positions для alien-ячеек
    positions_timer2.resume();
    tourism.setup_positions<dim>(locals.next);
    positions_timer2.stop();

    /// Этап 4. Создание новых ячеек
    geometry_timer.resume();
    setup_geometry<dim>(locals, tourism, count, op);
    geometry_timer.stop();

    /// Этап 5. Перемещение готовых ячеек
    remove_timer.resume();
    swap_list.move_elements(locals);
    locals.resize_amr(count.n_cells_short);
    tourism.resize_border();
    tourism.resize_aliens();
    remove_timer.stop();

    // Этап 6. Пересылка геометрии locals -> aliens
    exchange_timer.resume();
    tourism.send_geometry(locals);
    exchange_timer.stop();

    // Как будто особо и не нужно
    restore_timer.resume();
    tourism.restore_indices(locals);
    restore_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Statistics:       " << std::setw(9) << count_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Create SwapList:  " << std::setw(9) << swap_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Local Positions:  " << std::setw(9) << positions_timer1.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Alien Positions:  " << std::setw(9) << positions_timer2.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Setup Geometry:   " << std::setw(9) << geometry_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Remove undefined: " << std::setw(9) << remove_timer.milliseconds() << " ms\n";
        mpi::cout << "    Send geometry:    " << std::setw(9) << exchange_timer.milliseconds() << " ms\n";
        mpi::cout << "    Restore indices:  " << std::setw(9) << restore_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Специализация для пустых хранилищ
template<>
inline void apply_impl<0>(AmrCells &locals, const Distributor& op, Tourism& tourism) {
    throw std::runtime_error("MPI amr::apply<0> not implemented #1");
}

/// @brief Автоматический выбор размерности
inline void apply(AmrCells &locals, const Distributor& op, Tourism& tourism) {
    if (locals.empty()) {
        amr::apply_impl<0>(locals, op, tourism);
    } else {
        if (locals.dim() < 3) {
            amr::apply_impl<2>(locals, op, tourism);
        } else {
            amr::apply_impl<3>(locals, op, tourism);
        }
    }
}
#endif

} // namespace zephyr::mesh::amr