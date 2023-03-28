/// @file Файл содержит реализацию основной функции адаптации apply, которая производит
/// адаптацию сетки в соответствии с флагами адаптации. При этом предполагается, что
/// флаги адаптации сбалансированы.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/statistics.h>
#include <zephyr/mesh/refiner/impl/setup_positions.h>
#include <zephyr/mesh/refiner/impl/setup_geometry.h>
#include <zephyr/mesh/refiner/impl/restore_connections.h>
#include <zephyr/mesh/refiner/impl/remove_undefined.h>
#include <zephyr/mesh/refiner/impl/sorting.h>
#include <zephyr/mesh/refiner/impl/link_aliens.h>

namespace zephyr { namespace mesh { namespace impl {

using namespace ::zephyr::data;
using amrData;
#ifdef ZEPHYR_ENABLE_MULTITHREADING
using ::zephyr::multithreading::dummy_pool;
#endif

/// @brief Функция выполняет непостредственную адаптацию ячеек в хранилище в
/// соответствии с флагами адапатции. Предполагается, что флаги адаптации
/// сбалансированы.
/// @param cells Ссылка на хранилище ячеек
/// @param op Осуществляет распределение данных
/// @details Описание деталей алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т. д.
/// Этап 2. Определение позиций для создания новых ячеек (все ячейки будут
/// созданы за границами исходного хранилища)
/// Этап 3. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища.
/// Этап 4. Восстановление связей. На начале данного этапа все ячейки, включая
/// новые, ссылаются на старые индексы, а новые индексы указаны в amrData.next.
/// На данном этапе происходит выставление face.adjacent, которые указывают
/// на новые ячейки.
/// Этап 5. На начале этапа все ячейки правильно связаны, но внутри хранилища
/// часть старых ячеек (не листовых) являются неопределенными.
/// Алгоритм осуществляет удаление данных ячеек.
template<unsigned int dim = 0>
void apply(
        Storage &cells,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads = dummy_pool))
{
    using ::zephyr::performance::timer::Stopwatch;

    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch connections_timer;
    static Stopwatch clean_timer;
    static Stopwatch sort_timer;

    /// Этап 1. Сбор статистики
    count_timer.resume();
    Statistics<dim> count(cells if_multithreading(, threads));
    count_timer.stop();

    // Нечего адаптировать
    if (count.n_cells_large <= count.n_cells) return;

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    setup_positions(cells, count if_multithreading(, threads));
    positions_timer.stop();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry(cells, count, op if_multithreading(, threads));
    geometry_timer.stop();

    /// Этап 4. Восстановление соседства
    connections_timer.resume();
    restore_connections(cells, count if_multithreading(, threads));
    connections_timer.stop();

    /// Этап 5. Удаление неопределенных ячеек
    clean_timer.resume();
    remove_undefined(cells, count if_multithreading(, threads));
    clean_timer.stop();

    /// Этап 6. Сортировка ячеек по уровням (не обязательно)
    //sort_timer.resume();
    //sorting(cells);
    //sort_timer.stop();

#if CHECK_PERFORMANCE
    std::cout << "    Statistics elapsed: " << count_timer.times().wall() << "\n";
    std::cout << "    Positions elapsed: " << positions_timer.times().wall() << "\n";
    std::cout << "    Geometry elapsed: " << geometry_timer.times().wall() << "\n";
    std::cout << "    Connections elapsed: " << connections_timer.times().wall() << "\n";
    std::cout << "    Clean elapsed: " << clean_timer.times().wall() << "\n";
    std::cout << "    Sorting elapsed: " << sort_timer.times().wall() << "\n";
#endif
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
template<>
void apply<0>(
        Storage &cells,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads))
{
    if (cells.empty())
        return;

    auto dim = cells[0][element].dimension;

    if (dim < 3) {
        impl::apply<2>(cells, op if_multithreading(, threads));
    }
    else {
        impl::apply<3>(cells, op if_multithreading(, threads));
    }
}

#ifdef ZEPHYR_ENABLE_MPI
/// @brief Многопроцессорная версия функции, выполняет непостредственную
/// адаптацию ячеек в хранилище в соответствии с флагами адапатции.
/// Предполагается, что флаги адаптации сбалансированы.
/// @details Алгоритм практически совпадает с однопроцессорной версией,
/// поскольку недопускается огрубление ячеек, у которых сиблинги находятся
/// на различных процессах. Детали алгоритма:
/// Этап 1. Сбор данных о количестве ячеек для огрубления, разбиения и т. д.
/// Этап 2. Определение позиций для создания новых ячеек (все ячейки будут
/// созданы за границами исходного хранилища)
/// Этап 3. Создание геометрии ячеек, ячейки создаются на выделенных для них
/// местах за границами исходного хранилища.
/// Этап 4. Восстановление связей. На начале данного этапа все ячейки, включая
/// новые, ссылаются на старые индексы, а новые индексы указаны в amrData.next.
/// На данном этапе происходит выставление face.adjacent, которые указывают
/// на новые ячейки. Для граней между процессами ничего не происходит, на данном
/// этапе ссылки указвают на старые ячейки (до адаптации)
/// Этап 5. На начале этапа все ячейки правильно связаны, но внутри
/// хранилища часть старых ячеек (не листовых) являются неопределенными.
/// Алгоритм осуществляет удаление данных ячеек.
/// Этап 6. ???
template<unsigned int dim = 1234>
void apply(
        Decomposition &decomposition,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads = dummy_pool))
{
    using ::zephyr::performance::timer::Stopwatch;

    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch connections_timer;
    static Stopwatch clean_timer;
    static Stopwatch sort_timer;
    static Stopwatch link_timer;

    Storage& locals = decomposition.inner_elements();
    Storage& aliens = decomposition.outer_elements();
    unsigned int rank = decomposition.network().rank();

    /// Этап 1. Сбор статистики
    count_timer.resume();
    Statistics<dim> count(locals if_multithreading(, threads));
    //count.print(decomposition.network());
    count_timer.stop();

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    setup_positions(locals, count if_multithreading(, threads));
    positions_timer.stop();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry(locals, aliens, rank, count, op if_multithreading(, threads));
    geometry_timer.stop();

    /// Этап 4. Восстановление соседства
    connections_timer.resume();
    restore_connections(locals, aliens, rank, count if_multithreading(, threads));
    connections_timer.stop();

    /// Этап 5. Удаление неопределенных ячеек
    clean_timer.resume();
    remove_undefined(locals, count if_multithreading(, threads));
    clean_timer.stop();

    /// Этап 6. Пересылка и линковка alien ячеек
    link_timer.resume();
    link_aliens<dim>(decomposition if_multithreading(, threads));
    link_timer.stop();

#if CHECK_PERFORMANCE
    std::cout << "    Statistics elapsed: " << count_timer.times().wall() << "\n";
    std::cout << "    Positions elapsed: " << positions_timer.times().wall() << "\n";
    std::cout << "    Geometry elapsed: " << geometry_timer.times().wall() << "\n";
    std::cout << "    Connections elapsed: " << connections_timer.times().wall() << "\n";
    std::cout << "    Clean elapsed: " << clean_timer.times().wall() << "\n";
    std::cout << "    Sorting elapsed: " << sort_timer.times().wall() << "\n";
#endif
}

/// @brief Специализация для пустых хранилищ
template<>
void apply<0>(
        Decomposition &decomposition,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads))
{
    link_aliens<0>(decomposition if_multithreading(, threads));
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
template<>
void apply<1234>(
        Decomposition &decomposition,
        const DataDistributor& op
        if_multithreading(, ThreadPool& threads))
{
    Storage& cells = decomposition.inner_elements();

    if (cells.empty()) {
        impl::apply<0>(decomposition, op if_multithreading(, threads));
    }
    else {
        auto dim = cells[0][element].dimension;

        if (dim < 3) {
            impl::apply<2>(decomposition, op if_multithreading(, threads));
        } else {
            impl::apply<3>(decomposition, op if_multithreading(, threads));
        }
    }
}

#endif

} // namespace impl
} // namespace mesh
} // namespace zephyr