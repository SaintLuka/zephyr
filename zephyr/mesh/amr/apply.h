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
#include <zephyr/mesh/amr/restore_connections.h>
#include <zephyr/mesh/amr/remove_undefined.h>
#include <zephyr/mesh/euler/tourism.h>

namespace zephyr::mesh::amr {

/// @brief Функция выполняет непосредственную адаптацию ячеек в хранилище в
/// соответствии с флагами адаптации. Предполагается, что флаги адаптации
/// сбалансированы.
/// @param locals Ссылка на хранилище ячеек
/// @param op Осуществляет распределение данных
/// @details Описание деталей алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т. д.
/// Этап 2. Определение позиций для создания новых ячеек (все ячейки будут
/// созданы за границами исходного хранилища)
/// Этап 3. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища.
/// Этап 5. На начале этапа все ячейки правильно связаны, но внутри хранилища
/// часть старых ячеек (не листовых) являются неопределенными.
/// Алгоритм осуществляет удаление данных ячеек.
template<int dim>
void apply_impl(AmrCells &locals, const Distributor& op) {
    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch remove_timer;

    static Stopwatch create_swap_timer;
    static Stopwatch set_mapping_timer;
    static Stopwatch change_adjacent_timer;
    static Stopwatch swap_elements_timer;

    static AmrCells aliens;

    /// Этап 1. Сбор статистики
    count_timer.resume();
    const Statistics count(locals);
    //count.print();
    count_timer.stop();

    // Нечего адаптировать
    if (count.n_cells_large <= count.n_cells) return;

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    locals.resize_amr(count.n_cells_large);
    setup_positions<dim>(locals, count);
    positions_timer.stop();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry<dim>(locals, aliens, count, op);
    geometry_timer.stop();

    /// Этап 4. Удаление неопределенных ячеек
    remove_timer.resume();
    scrutiny_check(count.n_cells_large != count.n_cells, "Strange things");

    create_swap_timer.resume();
    SwapLists swap_list(count, locals.flag);
    create_swap_timer.stop();

    set_mapping_timer.resume();
    swap_list.set_mapping(locals);
    set_mapping_timer.stop();

    change_adjacent_timer.resume();
    change_adjacent<dim>(locals);
    change_adjacent_timer.stop();

    swap_elements_timer.resume();
    swap_list.move_elements(locals);
    swap_elements_timer.stop();

    locals.resize_amr(count.n_cells_short);
    remove_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Statistics:       " << std::setw(9) << count_timer.milliseconds() << " ms\n";
        mpi::cout << "    Setup Positions:  " << std::setw(9) << positions_timer.milliseconds() << " ms\n";
        mpi::cout << "    Setup Geometry:   " << std::setw(9) << geometry_timer.milliseconds() << " ms\n";
        mpi::cout << "    Remove undefined: " << std::setw(9) << remove_timer.milliseconds() << " ms\n";
        mpi::cout << "      Create SwapList: " << std::setw(8) << create_swap_timer.milliseconds() << " ms\n";
        mpi::cout << "      Setup mapping:   " << std::setw(8) << set_mapping_timer.milliseconds() << " ms\n";
        mpi::cout << "      Change adjacent: " << std::setw(8) << change_adjacent_timer.milliseconds() << " ms\n";
        mpi::cout << "      Swap elements:   " << std::setw(8) << swap_elements_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Автоматический выбор размерности
inline void apply(AmrCells &cells, const Distributor& op) {
    if (cells.empty())
        return;

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
/// @details Алгоритм практически совпадает с однопроцессорной версией,
/// поскольку не допускается огрубление ячеек, у которых сиблинги находятся
/// на различных процессах. Детали алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т.д.
/// Этап 2. Определение позиций для создания новых ячеек (все ячейки будут
/// созданы за границами исходного хранилища)
/// Этап 3. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища.
/// Этап 4. Восстановление связей. На начале данного этапа все ячейки, включая
/// новые, ссылаются на старые индексы, а новые индексы указаны в amrData.next.
/// На данном этапе происходит выставление face.adjacent, которые указывают
/// на новые ячейки. Для граней между процессами ничего не происходит, на данном
/// этапе ссылки указывают на старые ячейки (до адаптации)
/// Этап 5. На начале этапа все ячейки правильно связаны, но внутри
/// хранилища часть старых ячеек (не листовых) являются неопределенными.
/// Алгоритм осуществляет удаление данных ячеек.
template<int dim>
void apply_impl(Tourism& tourism, AmrCells &locals, AmrCells &aliens, const Distributor& op) {
    using utils::Stopwatch;

    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch connections_timer;
    static Stopwatch remove_timer;
    static Stopwatch exchange_timer;
    static Stopwatch alien1_timer;
    static Stopwatch alien2_timer;

    int rank = mpi::rank();

    /// Этап 1. Сбор статистики
    count_timer.resume();
    Statistics count(locals);
    count_timer.stop();

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    setup_positions<dim>(locals, count);
    positions_timer.stop();

    // Отправить / получить aliens.next.
    exchange_timer.resume();
    tourism.prepare<MpiTag::NEXT>(locals);
    auto send_next = tourism.isend<MpiTag::NEXT>();
    auto recv_next = tourism.irecv<MpiTag::NEXT>(aliens);
    exchange_timer.stop();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry<dim>(locals, aliens, rank, count, op);
    geometry_timer.stop();

    // Отправить / получить aliens.index. Если index < 0, то ячейка считается
    // неактуальной, то есть она переместилась, так мы понимаем, кто удалится.
    exchange_timer.resume();
    tourism.prepare<MpiTag::INDEX>(locals);
    auto send_idx = tourism.isend<MpiTag::INDEX>();
    auto recv_idx = tourism.irecv<MpiTag::INDEX>(aliens);

    // Дождаться (aliens.index, aliens.next).
    send_next.wait();
    recv_next.wait();
    send_idx.wait();
    recv_idx.wait();
    exchange_timer.stop();

    /// Этап 4. Восстановление соседства
    connections_timer.resume();
    restore_connections<dim>(locals, aliens, rank, count);
    connections_timer.stop();

    // Частичное восстановление border/aliens. Возможно, поскольку
    // в restore_connections верно проставлены (rank, index).
    // Восстановление нужно, чтобы в дальнейшем получить alien.next.
    alien1_timer.resume();
    tourism.build_aliens_basic(locals, aliens);
    alien1_timer.stop();

    /// Этап 5. Удаление неопределенных ячеек
    /// Сделать внутри MPI запрос чтобы узнать, куда переместились ячейки !
    remove_timer.resume();
    remove_undefined<dim>(tourism, locals, aliens, count);
    remove_timer.stop();

    // Здесь точно нужна, в отличае от предыдущей
    alien2_timer.resume();
    tourism.build_aliens(locals, aliens);
    alien2_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Statistics:     " << std::setw(11) << count_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Positions:      " << std::setw(11) << positions_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Geometry:       " << std::setw(11) << geometry_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Send/Recv:      " << std::setw(11) << exchange_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Connections:    " << std::setw(11) << connections_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Build aliens 1: " << std::setw(11) << alien1_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Build aliens 2: " << std::setw(11) << alien2_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Remove undef:   " << std::setw(11) << remove_timer.milliseconds_mpi() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Специализация для пустых хранилищ
template<>
inline void apply_impl<0>(Tourism& tourism, AmrCells &locals, AmrCells &aliens, const Distributor& op) {
    throw std::runtime_error("MPI amr::apply<0> not implemented #1");
}

/// @brief Автоматический выбор размерности
inline void apply(Tourism& tourism, AmrCells &locals, AmrCells &aliens, const Distributor& op) {
    if (locals.empty()) {
        amr::apply_impl<0>(tourism, locals, aliens, op);
    } else {
        if (locals.dim() < 3) {
            amr::apply_impl<2>(tourism, locals, aliens, op);
        } else {
            amr::apply_impl<3>(tourism, locals, aliens, op);
        }
    }
}
#endif

} // namespace zephyr::mesh::amr