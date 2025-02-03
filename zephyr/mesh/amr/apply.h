/// @file Файл содержит реализацию основной функции адаптации apply, которая производит
/// адаптацию сетки в соответствии с флагами адаптации. При этом предполагается, что
/// флаги адаптации сбалансированы.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <iomanip>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/statistics.h>
#include <zephyr/mesh/amr/setup_positions.h>
#include <zephyr/mesh/amr/setup_geometry.h>
#include <zephyr/mesh/amr/restore_connections.h>
#include <zephyr/mesh/amr/remove_undefined.h>
#include <zephyr/mesh/amr/sorting.h>
#include <zephyr/mesh/amr/link_aliens.h>

namespace zephyr::mesh {
// Сейчас нужно для MPI, там функция обменов
class EuMesh;
}

namespace zephyr::mesh::amr {

using utils::Stopwatch;

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
template<int dim>
void apply_impl(AmrStorage &cells, const Distributor& op) {
    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch connections_timer;
    static Stopwatch clean_timer;
    static Stopwatch sort_timer;

    AmrStorage aliens;

    /// Этап 1. Сбор статистики
    count_timer.resume();
    Statistics count(cells);
    count_timer.stop();

    // Нечего адаптировать
    if (count.n_cells_large <= count.n_cells) return;

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    setup_positions<dim>(cells, count);
    positions_timer.stop();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry<dim>(cells, count, op);
    geometry_timer.stop();

    /// Этап 4. Восстановление соседства
    connections_timer.resume();
    restore_connections<dim>(cells, aliens, 0, count);
    connections_timer.stop();

    /// Этап 5. Удаление неопределенных ячеек
    clean_timer.resume();
    remove_undefined<dim>(cells, count);
    clean_timer.stop();

    for (int idx = 0; idx < cells.size(); ++idx) {
        cells[idx].index = idx;
    }

    /// Этап 6. Сортировка ячеек по уровням (не обязательно)
    //sort_timer.resume();
    //sorting(cells);
    //sort_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "    Statistics elapsed:  " << std::setw(10) << count_timer.milliseconds() << " ms\n";
        std::cout << "    Positions elapsed:   " << std::setw(10) << positions_timer.milliseconds() << " ms\n";
        std::cout << "    Geometry elapsed:    " << std::setw(10) << geometry_timer.milliseconds() << " ms\n";
        std::cout << "    Connections elapsed: " << std::setw(10) << connections_timer.milliseconds() << " ms\n";
        std::cout << "    Clean elapsed:       " << std::setw(10) << clean_timer.milliseconds() << " ms\n";
        std::cout << "    Sorting elapsed:     " << std::setw(10) << sort_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Автоматический выбор размерности
void apply(AmrStorage &cells, const Distributor& op) {
    if (cells.empty())
        return;

    auto dim = cells[0].dim;

    if (dim < 3) {
        amr::apply_impl<2>(cells, op);
    }
    else {
        amr::apply_impl<3>(cells, op);
    }
}

#ifdef ZEPHYR_MPI
/// @brief Многопроцессорная версия функции, выполняет непостредственную
/// адаптацию ячеек в хранилище в соответствии с флагами адапатции.
/// Предполагается, что флаги адаптации сбалансированы.
/// @details Алгоритм практически совпадает с однопроцессорной версией,
/// поскольку не допускается огрубление ячеек, у которых сиблинги находятся
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
///
/// \param mesh Сейчас используется, поскольку нужна операция обменов,
/// потом желательно выкинуть
template<int dim>
void apply_impl(
        AmrStorage &locals, AmrStorage &aliens,
        const Distributor& op,
        EuMesh& mesh)
{
    using zephyr::utils::Stopwatch;

    static Stopwatch count_timer;
    static Stopwatch positions_timer;
    static Stopwatch geometry_timer;
    static Stopwatch connections_timer;
    static Stopwatch clean_timer;
    static Stopwatch sort_timer;
    static Stopwatch link_timer;

    int rank = mpi::rank();

    // Я не знаю флаги, надо получить.
    mesh.exchange();

    /// Этап 1. Сбор статистики
    count_timer.resume();
    Statistics count(locals);
    //mpi::for_each([&]() {
    //    std::cout << "rank " << mpi::rank() << "\n";
    //    count.print();
    //});
    count_timer.stop();

    /// Этап 2. Распределяем места для новых ячеек
    positions_timer.resume();
    setup_positions<dim>(locals, count);
    positions_timer.stop();

    // Отправить MPI locals.next (!!!!)
    mesh.exchange();

    /// Этап 3. Восстановление геометрии
    geometry_timer.resume();
    setup_geometry<dim>(locals, aliens, rank, count, op);
    geometry_timer.stop();

    // Получить по MPI aliens.index, aliens.next
    // Если index < 0, то ячейка считается неактуальной, то есть
    // она переместилась, так мы отличаем тех, кто остался от тех, кто удалится.
    // Ну и флаги нужны, конечно. Но они и так есть, вроде как.
    mesh.exchange();

    /// Этап 4. Восстановление соседства
    connections_timer.resume();
    restore_connections<dim>(locals, aliens, rank, count);
    connections_timer.stop();

    /// Этап 5. Удаление неопределенных ячеек
    /// Сделать внутри MPI запрос чтобы узнать, куда переместились ячейки !
    clean_timer.resume();
    remove_undefined<dim>(locals, aliens, count, mesh);
    clean_timer.stop();

    mesh.build_aliens();
    mesh.exchange();

    /// Этап 6. Пересылка и линковка alien ячеек
    //link_timer.resume();
    //link_aliens<dim>(decomposition if_multithreading(, threads));
    //link_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "    Statistics elapsed: " << count_timer.times().wall() << "\n";
        std::cout << "    Positions elapsed: " << positions_timer.times().wall() << "\n";
        std::cout << "    Geometry elapsed: " << geometry_timer.times().wall() << "\n";
        std::cout << "    Connections elapsed: " << connections_timer.times().wall() << "\n";
        std::cout << "    Clean elapsed: " << clean_timer.times().wall() << "\n";
        std::cout << "    Sorting elapsed: " << sort_timer.times().wall() << "\n";
    }
    ++counter;
#endif
}

/// @brief Специализация для пустых хранилищ
template<>
void apply_impl<0>(
        AmrStorage &locals, AmrStorage &aliens,
        const Distributor& op,
        EuMesh& mesh)
{
    // TODO LINK ALIENS, WHAT??
    //link_aliens<0>(decomposition);
}

/// @brief Автоматический выбор размерности
void apply(
        AmrStorage &locals, AmrStorage &aliens,
        const Distributor& op,
        EuMesh& mesh)
{
    if (locals.empty()) {
        amr::apply_impl<0>(locals, aliens, op, mesh);
    }
    else {
        auto dim = locals[0].dim;

        if (dim < 3) {
            amr::apply_impl<2>(locals, aliens, op, mesh);
        } else {
            amr::apply_impl<3>(locals, aliens, op, mesh);
        }
    }
}
#endif

} // namespace zephyr::mesh::amr