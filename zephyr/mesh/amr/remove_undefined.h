// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr::mesh::amr {

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для одной ячейки. Предполагается, что в поле element.next записан индекс,
/// куда ячейка будет перемещена в дальнейшем
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
                       "change_adjacent_one() error: adjacent.index out of range");

        locals.faces.adjacent.index[iface] = locals.next[locals.faces.adjacent.index[iface]];
        locals.faces.adjacent.basic[iface] = locals.next[ic];
    }
}

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для одной ячейки. Предполагается, что в поле element.next записан индекс,
/// куда ячейка будет перемещена в дальнейшем
inline void change_adjacent_one2(index_t ic, AmrCells& locals, AmrCells& aliens) {
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

        if (locals.faces.adjacent.alien[iface] < 0) {
            // Локальный сосед
            scrutiny_check(locals.faces.adjacent.index[iface] >= 0 &&
                           locals.faces.adjacent.index[iface] < locals.size(),
                           "change_adjacent_one() error: adjacent.index out of range");

            locals.faces.adjacent.index[iface] = locals.next[locals.faces.adjacent.index[iface]];
        } else {
            // Удаленный сосед
            scrutiny_check(locals.faces.adjacent.alien[iface] < aliens.size(),
                           "change_adjacent_one() error: adjacent.alien out of range");

            locals.faces.adjacent.index[iface] = aliens.next[locals.faces.adjacent.alien[iface]];
        }
        locals.faces.adjacent.basic[iface] = locals.next[ic];
    }
}

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в однопоточном режиме. Предполагается, что в поле
/// element.next записан индекс, куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
template<int dim>
void change_adjacent(AmrCells& cells) {
    const index_t* index = cells.index.data();
    const index_t* next = cells.next.data();

    const Boundary* boundary = cells.faces.boundary.data();
    index_t* adj_index = cells.faces.adjacent.index.data();
    index_t* adj_basic = cells.faces.adjacent.basic.data();

    auto change_adjacent_fast = [index, next, boundary, adj_index, adj_basic](index_t ic) {
        if (index[ic] < 0) { return; }
        constexpr index_t face_count = FpC(dim) * FpF(dim);
        index_t face_beg = face_count * ic;
        index_t face_end = face_count * (ic + 1);
        for (index_t iface = face_beg; iface < face_end; ++iface) {
            if (boundary[iface] != Boundary::UNDEFINED) {
                adj_basic[iface] = next[ic];
                // Граничные указывают на ячейку
                if (boundary[iface] != Boundary::ORDINARY &&
                    boundary[iface] != Boundary::PERIODIC) {
                    adj_index[iface] = next[ic];
                    }
                else {
                    adj_index[iface] = next[adj_index[iface]];
                }
            }
        }
    };

    threads::parallel_for(
            index_t{0}, index_t{cells.size()},
            change_adjacent_fast);
}

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в однопоточном режиме. Предполагается, что в поле
/// element.next записан индекс, куда ячейка будет перемещена в дальнейшем
inline void change_adjacent(AmrCells& locals, AmrCells& aliens) {
    threads::parallel_for(
            index_t{0}, index_t{locals.size()},
            change_adjacent_one2,
            std::ref(locals), std::ref(aliens));
}


/// @brief Вспомогательная структура, содержит два списка индексов: неопределенные
/// ячейки и актуальные ячейки. После обмена местами неопределенных и актуальных
/// ячеек, все неактуальные ячейки оказываются в конце хранилища.
struct SwapLists {
    /// @brief Список индексов неопределенных ячеек, начиная с начала хранилища
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
    /// @param count Статистика адаптации
    /// @param flag Массив флагов ячеек (можно расширенный или нет, то есть flag.size() >= n_cells)
    SwapLists(const Statistics& count, const std::vector<int> flag) {
        index_t max_swap_count = count.n_cells_large - count.n_cells_short;
        if (max_swap_count < 1) { return; }

        undefined_cells.reserve(max_swap_count);
        for (index_t ic = 0; ic < count.n_cells_short; ++ic) {
            if (flag[ic] != 0) {
                undefined_cells.push_back(ic);
                if (undefined_cells.size() >= max_swap_count) {
                    break;
                }
            }
        }
        index_t swap_count = undefined_cells.size();

        actual_cells.reserve(swap_count);
        index_t jc = count.n_cells_large - 1;
        while (actual_cells.size() < swap_count && jc >= count.n_cells) {
            actual_cells.push_back(jc);
            --jc;
        }
        while (actual_cells.size() < swap_count && jc >= 0) {
            if (flag[jc] == 0) {
                actual_cells.push_back(jc);
            }
            --jc;
        }

        scrutiny_check(actual_cells.size() == undefined_cells.size(), "actual/undefined size mismath");
    }

#if SCRUTINY
    /// @brief Проверяет свойства перестановки: перестановка состоит только из
    /// транспозиций пар элементов, при этом один элемент в паре должен быть
    /// актуальным, а один неопределенным, после перестановки актуальный элемент
    /// всегда должен оказываться ближе к началу списка, чем неопределенный
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
                    // Актуальная ячейка переносится только в начало
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

    /// @brief Устанавливает следующие позиции ячеек в однопоточном режиме
    /// @details После выполнения операции поле element.next у ячеек указывает
    /// на следующее положение ячейки в хранилище
    void set_mapping(AmrCells& locals) const {
        // Устанавливает тождественную перестановку для всех ячеек
        threads::parallel_for(
            index_t{0}, index_t{locals.size()},
            [&locals](index_t ic) {
                locals.next[ic] = ic;
            });

        if (actual_cells.empty()) { return; }

        // Устанавливает следующую позицию для части ячеек
        // Индекс проходит по массивам actual_cells, undefined_cells.
        threads::parallel_for(
            index_t{0}, index_t{size()},
            [this, &locals](index_t i) {
                scrutiny_check(i < actual_cells.size(),    "swap_mapping out of range #1")
                scrutiny_check(i < undefined_cells.size(), "swap_mapping out of range #2")

                const index_t ai = actual_cells[i];
                const index_t ui = undefined_cells[i];

                scrutiny_check(ai < locals.size(), "swap_mapping out of range #3")
                scrutiny_check(ui < locals.size(), "swap_mapping out of range #4")

                locals.next[ui] = ai;
                locals.next[ai] = ui;
            });
#if SCRUTINY
        check_mapping(locals);
#endif
    }

    /// @brief Выполняет перестановку элементов в соответствии с next.
    /// @details Данные актуальной ячейки перемещаются на место неактуальной
    /// ячейки, индексы смежности должны быть выставлены ранее.
    void move_elements(AmrCells &cells) const {
        /// @brief Выполняет перемещение элемента в соответствии с индексом
        /// в поле element.next для одной ячейки.
        /// @details Данные актуальной ячейки перемещаются на место неактуальной
        /// ячейки, индексы смежности не изменяются, они должны быть выставлены
        /// заранее.
        auto move_cell = [this, &cells](index_t i) {
            index_t from = actual_cells[i];
            index_t to   = undefined_cells[i];
            cells.move_item(from, to);
            cells.copy_data(from, to);
        };

        threads::parallel_for(index_t{0}, index_t{size()}, move_cell);
    }

    /// @brief Число ячеек для перестановки
    index_t size() const {
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
    static Stopwatch swap_elements_timer;

    if (count.n_cells_large == count.n_cells_short) {
        return;
    }

    create_swap_timer.resume();
    SwapLists swap_list(count, cells.flag);
    create_swap_timer.stop();

    set_mapping_timer.resume();
    swap_list.set_mapping(cells);
    set_mapping_timer.stop();

    change_adjacent_timer.resume();
    change_adjacent<dim>(cells);
    change_adjacent_timer.stop();

    swap_elements_timer.resume();
    swap_list.move_elements(cells);
    swap_elements_timer.stop();

    cells.resize_amr(count.n_cells_short);

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "      Create SwapList: " << std::setw(8) << create_swap_timer.milliseconds() << " ms\n";
        mpi::cout << "      Setup mapping:   " << std::setw(8) << set_mapping_timer.milliseconds() << " ms\n";
        mpi::cout << "      Change adjacent: " << std::setw(8) << change_adjacent_timer.milliseconds() << " ms\n";
        mpi::cout << "      Swap elements:   " << std::setw(8) << swap_elements_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

// ============================================================================
//                                 MPI VERSION
// ============================================================================
#ifdef ZEPHYR_MPI
template<int dim>
void remove_undefined(Tourism& tourism, AmrCells& locals, AmrCells& aliens, const Statistics &count) {
    static Stopwatch create_swap_timer;
    static Stopwatch set_mapping_timer;
    static Stopwatch change_adjacent_timer;
    static Stopwatch swap_elements_timer;
    static Stopwatch send_next_timer;

    create_swap_timer.resume();
    index_t max_swap_count = count.n_cells_large - count.n_cells_short;
    SwapLists swap_list(locals, count.n_cells_short, max_swap_count);
    create_swap_timer.stop();

    set_mapping_timer.resume();
    swap_list.set_mapping(locals);
    set_mapping_timer.stop();

    // Перед remove_undefined сделан частичный build_aliens,
    // получим актуальные aliens.next.
    send_next_timer.resume();
    tourism.sync<MpiTag::NEXT>(locals, aliens);
    send_next_timer.stop();

    change_adjacent_timer.resume();
    change_adjacent(locals, aliens);
    change_adjacent_timer.stop();

    swap_elements_timer.resume();
    swap_list.move_elements(locals);
    swap_elements_timer.stop();

    locals.resize_amr(count.n_cells_short);

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "      Create SwapList: " << std::setw(8) << create_swap_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "      Setup mapping:   " << std::setw(8) << set_mapping_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "      Send/Recv next:  " << std::setw(8) << send_next_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "      Change adjacent: " << std::setw(8) << change_adjacent_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "      Swap elements:   " << std::setw(8) << swap_elements_timer.milliseconds_mpi() << " ms\n";
    }
    ++counter;
#endif
}
#endif
} // namespace zephyr::mesh::amr