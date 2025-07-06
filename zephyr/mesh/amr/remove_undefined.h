/// @file Файл содержит реализацию функции remove_undefined, которая очищает хранилище
/// от неопределенных (не листовых ячеек).
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/utils/stopwatch.h>

namespace zephyr::mesh::amr {

using zephyr::utils::Stopwatch;

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для одной ячейки. Предполагается, что в поле element.next записан индекс,
/// куда ячейка будет перемещена в дальнейшем
/// @param cell Целевая ячейка
/// @param locals Хранилище ячеек
inline void change_adjacent_one1(index_t ic, AmrCells& locals) {
    if (locals.is_undefined(ic)) { return; }

    for (auto iface: locals.faces_range(ic)) {
        if (locals.faces.is_undefined(iface)) {
            continue;
        }
        // Граничные указывают на ячейку
        if (locals.faces.is_boundary(iface)) {
            locals.faces.adjacent.index[iface] = locals.next[ic];
            locals.faces.adjacent.basic[iface] = locals.next[ic];
            continue;
        }

        scrutiny_check(locals.faces.adjacent.alien[iface] < 0,
                       "change_adjacent_one() error: only local neighbors")

        scrutiny_check(locals.faces.adjacent.index[iface] >= 0 &&
                       locals.faces.adjacent.index[iface] < locals.size(),
                       "change_adjacent_partial() error: adjacent.index out of range");

        locals.faces.adjacent.index[iface] = locals.next[locals.faces.adjacent.index[iface]];
        locals.faces.adjacent.basic[iface] = locals.next[ic];
    }
}

/*
/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для одной ячейки. Предполагается, что в поле element.next записан индекс,
/// куда ячейка будет перемещена в дальнейшем
/// @param cell Целевая ячейка
inline void change_adjacent_one2(AmrStorage::Item& cell, AmrStorage& locals, AmrStorage& aliens) {
    if (cell.is_undefined()) { return; }

    for (BFace &face: cell.faces) {
        if (face.is_undefined()) {
            continue;
        }
        // Граничные указывают на ячейку
        if (face.is_boundary()) {
            face.adjacent.index = cell.next;
            continue;
        }

        if (face.adjacent.alien < 0) {
            // Локальный сосед
            scrutiny_check(face.adjacent.index >= 0 && face.adjacent.index < locals.size(),
                           "change_adjacent_partial() error: adjacent.index out of range");

            face.adjacent.index = locals[face.adjacent.index].next;
        } else {
            // Удаленный сосед
            scrutiny_check(face.adjacent.alien >= 0 && face.adjacent.alien < aliens.size(),
                           "change_adjacent_partial() error: adjacent.alien out of range");

            face.adjacent.index = aliens[face.adjacent.alien].next;
        }
    }
}
*/

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в однопоточном режиме. Предполагается, что в поле
/// element.next записан индекс, куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
void change_adjacent(AmrCells& cells) {
    threads::parallel_for(
            index_t{0}, index_t{cells.size()},
            change_adjacent_one1,
            std::ref(cells));
}

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в однопоточном режиме. Предполагается, что в поле
/// element.next записан индекс, куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
void change_adjacent2(AmrCells& cells) {
    throw std::runtime_error("change_adjacent() error: Not implemented");
    /*
    threads::for_each<5>(
            cells.begin(), cells.end(),
            change_adjacent_one2,
            std::ref(cells));
            */
}

/// @brief Выполняет перемещение элемента в соответствии с индексом
/// в поле element.next для одной ячейки.
/// @details Данные актуальной ячейки перемещаются на место неактуальной
/// ячейки, индексы смежности не изменяются, они должны быть выставлены
/// заранее.
inline void move_cell(index_t ic, AmrCells& cells) {
    cells.move_item(ic);
    cells.copy_data(ic, cells.next[ic]);
}

/// @brief Вспомогательная структура, содержит два списка индексов: неопределенные
/// ячейки и актуальные ячейки. После обмена местами неопределенных и актуальных
/// ячеек, все неактуальные ячейки оказываются в конце хранилища.
struct SwapLists {
    /// @brief Список индеков неопределенных ячеек, начиная с начала хранилища
    /// @details Массив может содержать не все неопределенные ячейки, которые
    /// есть в хранилище
    std::vector<index_t> undefined_cells;

    /// @brief Список индексов актуальных ячеек, начиная с конца хранилища,
    /// размер массива обязательно совпадает с размером массива неопределенных
    /// ячеек. При перестановке элементов хранилища с индексами undefined_cells[i]
    /// и actual_cells[i] для всех i, все неопределенные ячейки хранилища
    /// должны оказаться в конце.
    std::vector<index_t> actual_cells;

    /// @brief Конструктор
    /// @param locals Хранилище ячеек
    /// @param max_index Максимальный индекс, который будет в массиве неопределенных
    /// ячеек, правильное задание - размер хранилища после удаления всех неопределенных
    /// элементов. Правильное задание параметра гарантирует, что будет выполняться
    /// свойство undefined_cells[i] < actual_cells[i] для любых i
    /// @param max_swap_count Максимальное число элементов, для которых может
    /// потребоваться перестановка, размеры списков ограничены данным числом.
    // TODO: Многопоточная версия конструктора
    SwapLists(AmrCells& locals, index_t max_index, index_t max_swap_count) {
        if (max_swap_count < 1) { return; }

        undefined_cells.reserve(max_swap_count);
        for (index_t ic = 0; ic < max_index; ++ic) {
            if (locals.is_undefined(ic)) {
                undefined_cells.push_back(ic);
                if (undefined_cells.size() >= max_swap_count) {
                    break;
                }
            }
        }

        actual_cells.reserve(undefined_cells.size());
        for (index_t jc = locals.size() - 1; jc >= 0; --jc) {
            if (locals.is_actual(jc)) {
                actual_cells.push_back(jc);
                if (actual_cells.size() >= undefined_cells.size()) {
                    break;
                }
            }
        }
    }

#if SCRUTINY
    /// @brief Проверяет свойства перестановки: перестановка состоит только из
    /// транспозций пар элементов, при этом один элемент в паре должен быть актуальным,
    /// а один неопределенным, после перестановки актуальный элемент всегда должен
    /// оказываться ближе к началу списка, чем неопределенный
    static void check_mapping(AmrCells& locals) {
        for (index_t i = 0; i < locals.size(); ++i) {
            auto j = locals.next[i];
            if (i != locals.next[j]) {
                // Перестановка не является транспозицией
                std::cout << i << " " << locals.next[i] << "\n";
                std::cout << j << " " << locals.next[j] << "\n";
                throw std::runtime_error("Only swaps are allowed");
            }
            if (i != j) {
                if (locals.is_actual(i) && locals.is_actual(j)) {
                    // Попытка обменять две актуальные ячейки
                    throw std::runtime_error("Swap two actual cells");
                }
                if (i < j && locals.is_actual(i)) {
                    // Актуальня ячейка переносится только в начало
                    throw std::runtime_error("Wrong swap #1");
                }
                if (i > j && locals.is_undefined(i)) {
                    // Неопределенная ячейка переносится только в конец
                    throw std::runtime_error("Wrong swap #2");
                }
            }
        }
    }
#endif

    /// @brief Устанавливает тождественную перестановку для части ячеек
    /// @param cells Хранилище ячеек
    /// @param from, to Диапазон ячеек в хранилище
    static void set_identical_mapping_partial(AmrCells& cells, index_t from, index_t to) {
        for (index_t ic = from; ic < to; ++ic) {
            cells.next[ic] = ic;
        }
    }

    /// @brief Устанавливает следующую позицию для части ячеек
    /// @param locals Ссылка на хранилище ячеек
    /// @param from, to Диапазон индексов в массивах actual_cells и undefined_cells
    void set_swap_mapping_partial(AmrCells& locals, const index_t from, const index_t to) const {
        for (index_t i = from; i < to; ++i) {
            scrutiny_check(i < actual_cells.size(),    "swap_mapping out of range #1")
            scrutiny_check(i < undefined_cells.size(), "swap_mapping out of range #2")

            const index_t ai = actual_cells[i];
            const index_t ui = undefined_cells[i];

            scrutiny_check(ai < locals.size(), "swap_mapping out of range #3")
            scrutiny_check(ui < locals.size(), "swap_mapping out of range #4")

            locals.next[ui] = ai;
            locals.next[ai] = ui;
        }
    }

    /// @brief Устанавливает следующие позиции ячеек в однопоточном режиме
    /// @details После выполнения операции поле element.next у ячеек указывает
    /// на следующее положение ячейки в хранилище
    void set_mapping(AmrCells& locals) const {
        set_identical_mapping_partial(locals, 0, locals.size());

        if (actual_cells.empty()) { return; }

        set_swap_mapping_partial(locals, 0, size());
#if SCRUTINY
        check_mapping(locals);
#endif
    }

#ifdef ZEPHYR_MULTITHREADING
    /// @brief Устанавливает следующие позиции ячеек в многопоточном режиме
    /// @details После выполнения операции поле element.index у ячеек указывает
    /// на следующее положение ячейки в хранилище
    void set_mapping(AmrStorage& cells, ThreadPool& threads) const {
        auto num_tasks = threads.size();
        if (num_tasks < 2) {
            set_mapping(cells);
            return;
        }
        std::vector<std::future<void>> results(num_tasks);

        index_t bin = cells.size() / num_tasks + 1;
        index_t pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(
                    &SwapLists::set_identical_mapping_partial,
                    this, std::ref(cells),
                    pos, std::min(pos + bin, cells.size())
            );
            pos += bin;
        }

        for (auto &result: results) result.get();

        bin = size() / num_tasks + 1;
        pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(
                    &SwapLists::set_swap_mapping_partial,
                    this, std::ref(cells),
                    pos, std::min(pos + bin, size())
            );
            pos += bin;
        }

        for (auto &result: results) result.get();
#if SCRUTINY
        check_mapping(cells);
#endif
    }
#endif

    /// @brief Выполняет перестановку элементов в соответствии с индексом
    /// в поле element.next в однопоточном режиме.
    /// @details Данные актуальной ячейки перемещаются на место неактуальной
    /// ячейки, индексы смежности не изменяются, они должны быть выставлены
    /// заранее.
    void move_elements(AmrCells &cells) const {
        threads::for_each<5>(
                actual_cells.begin(), actual_cells.end(),
                move_cell, std::ref(cells));
    }

    /// @brief Число ячеек для перестановки
    inline index_t size() const {
        return actual_cells.size();
    }
};

/// @brief Осуществляет удаление неопределенных ячеек из хранилища
/// после операции адаптации (!) с сохранением смежности.
/// @param cells Ссылка на хранилище
/// @param count Статистика адаптации
/// @details Функця успешно работает только после выполнения операции
/// адаптации. Параметр count необходим для получения данных о размерах
/// хранилища после адаптации. Удаление ячеек осуществляется следующим образом:
/// 1. На первом этапе составляются списки неопределенных ячеек, начиная с начала
/// хранилища, а также актуальных ячеек, начиная с конца хранилища. На данном
/// этапе используется конструктор класса SwapLists
/// 2. Далее для ячеек изменяется поле element.next, которое после выполнения
/// функции set_mapping равно индексу, куда необходимо переместить ячейку.
/// 3. Выполняется функция change_adjacent, которая меняет индексы смежности
/// в соответствии с element.next
/// 4. Актуальные ячейки перемещаются на место неопределенных.
/// 5. Хранилище меняет размер, все неопределенные ячейки остаются за пределами
/// хранилища.
template<int dim>
void remove_undefined(AmrCells &cells, const Statistics &count) {
    static Stopwatch create_swap_timer;
    static Stopwatch set_mapping_timer;
    static Stopwatch change_adjacent_timer;
    static Stopwatch move_elements_timer;

    if (count.n_cells_large == count.n_cells_short) {
        return;
    }

    create_swap_timer.resume();
    index_t max_swap_count = count.n_cells_large - count.n_cells_short;
    SwapLists swap_list(cells, count.n_cells_short, max_swap_count);
    create_swap_timer.stop();

    set_mapping_timer.resume();
    swap_list.set_mapping(cells);
    set_mapping_timer.stop();

    change_adjacent_timer.resume();
    change_adjacent(cells);
    change_adjacent_timer.stop();

    move_elements_timer.resume();
    swap_list.move_elements(cells);
    move_elements_timer.stop();

    cells.resize_amr(count.n_cells_short);

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "      Create SwapList elapsed: " << std::setw(10) << create_swap_timer.milliseconds() << " ms\n";
        std::cout << "      Set mapping elapsed:     " << std::setw(10) << set_mapping_timer.milliseconds() << " ms\n";
        std::cout << "      Change adjacent elapsed: " << std::setw(10) << change_adjacent_timer.milliseconds() << " ms\n";
        std::cout << "      Move elements elapsed:   " << std::setw(10) << move_elements_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/* ============================================================================
                                    MPI VERSION
 ==============================================================================
template<int dim>
void remove_undefined(AmrCells& cells, const Statistics &count, EuMesh& mesh) {
    static Stopwatch create_swap_timer;
    static Stopwatch set_mapping_timer;
    static Stopwatch change_adjacent_timer;
    static Stopwatch move_elements_timer;

    create_swap_timer.resume();
    index_t max_swap_count = count.n_cells_large - count.n_cells_short;
    SwapLists swap_list(cells, count.n_cells_short, max_swap_count);
    create_swap_timer.stop();

    set_mapping_timer.resume();
    swap_list.set_mapping(cells);
    set_mapping_timer.stop();

    // получить element.next у alien ячеек
    mesh.build_aliens();
    mesh.sync();

    change_adjacent_timer.resume();
    change_adjacent(cells);
    change_adjacent_timer.stop();

    move_elements_timer.resume();
    swap_list.move_elements(cells);
    move_elements_timer.stop();

    cells.resize(count.n_cells_short);

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "      Create SwapList elapsed: " << std::setw(10) << create_swap_timer.milliseconds() << " ms\n";
        std::cout << "      Set mapping elapsed:     " << std::setw(10) << set_mapping_timer.milliseconds() << " ms\n";
        std::cout << "      Change adjacent elapsed: " << std::setw(10) << change_adjacent_timer.milliseconds() << " ms\n";
        std::cout << "      Move elements elapsed:   " << std::setw(10) << move_elements_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}
 */

} // namespace zephyr::mesh::amr