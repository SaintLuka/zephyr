/// @file Файл содержит реализацию функции remove_undefined, которая очищает хранилище
/// от неопределенных (не листовых ячеек).
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr { namespace mesh { namespace amr {

using zephyr::utils::Stopwatch;

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для части ячеек. Предполагается, что в поле element.index записан индекс,
/// куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
/// @param from, to Диапазон ячеек, по которому осуществляется обход
void change_adjacent_partial(Storage& cells, int from, int to) {
    int n_cells = cells.size();
    for (int ic = from; ic < to; ++ic) {
        for (Face &face: cells[ic].geom().faces) {
            if (face.is_undefined() or face.is_boundary()) continue;
            if (face.adjacent.ghost < 0) {
                // Локальный сосед
#if SCRUTINY
                if (face.adjacent.index >= n_cells) {
                    throw std::runtime_error("change_adjacent_partial() error: adjacent.index out of range");
                }
#endif
                face.adjacent.index = cells[face.adjacent.index].geom().index;
            }

        }
    }
}

/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в однопоточном режиме. Предполагается, что в поле
/// element.index записан индекс, куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
void change_adjacent(Storage& cells) {
    change_adjacent_partial(cells, 0, cells.size());
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Функция меняет индексы соседей через грань на новые, выполняется
/// для всех ячеек в многопоточном режиме. Предполагается, что в поле
/// element.index записан индекс, куда ячейка будет перемещена в дальнейшем
/// @param cells Хранилище ячеек
void change_adjacent(Storage& cells, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        change_adjacent(cells);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = cells.size() / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(change_adjacent_partial,
                              std::ref(cells),
                              pos, std::min(pos + bin, cells.size())
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }
}
#endif

/// @struct Вспомогательная структура, содержит два списка индексов: неопределенные
/// ячейки и актуальные ячейки. После обмена местами неопределенных и актуальных
/// ячеек, все неактуальные ячейки оказываются в конце хранилища.
struct SwapLists {
    /// @brief Список индеков неопределенных ячеек, начиная с начала хранилища
    /// @details Массив может содержать не все неопределенные ячейки, которые
    /// есть в хранилище
    std::vector<int> undefined_cells;

    /// @brief Список индексов актуальных ячеек, начиная с конца хранилища,
    /// размер массива обязательно совпадает с размером массива неопределенных
    /// ячеек. При перестановке элементов хранилища с индексами undefined_cells[i]
    /// и actual_cells[i] для всех i, все неопределенные ячейки хранилища
    /// должны оказаться в конце.
    std::vector<int> actual_cells;

    /// @brief Конструктор
    /// @param cells Хранилище ячеек
    /// @param max_index Максимальный индекс, который будет в массиве неопределенных
    /// ячеек, правильное задание - размер хранилища после удаления всех неопределенных
    /// элементов. Правильное задание параметра гарантирует, что будет выполняться
    /// свойство undefined_cells[i] < actual_cells[i] для любых i
    /// @param max_swap_count Максимальное число элементов, для которых может
    /// потребоваться перестановка, размеры списков ограничены данным числом.
    /// TODO: Многопоточная версия конструктора
    SwapLists(Storage& cells, int max_index, int max_swap_count) {
        undefined_cells.reserve(max_swap_count);
        for (int ic = 0; ic < max_index; ++ic) {
            if (cells[ic].is_undefined()) {
                undefined_cells.push_back(ic);
                if (undefined_cells.size() >= max_swap_count) {
                    break;
                }
            }
        }

        actual_cells.reserve(undefined_cells.size());
        for (int jc = cells.size() - 1; jc >= 0; --jc) {
            if (cells[jc].is_actual()) {
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
    void check_mapping(Storage& cells) const {
        for (int i = 0; i < cells.size(); ++i) {
            auto j = cells[i].geom().index;
            if (i != cells[j].geom().index) {
                // Перестановка не является транспозицией
                std::cout << i << " " << cells[i].geom().index << "\n";
                std::cout << j << " " << cells[j].geom().index << "\n";
                throw std::runtime_error("Only swaps are allowed");
            }
            if (i != j) {
                if (cells[i].is_actual() && cells[j].is_actual()) {
                    // Попытка обменять две актуальные ячейки
                    throw std::runtime_error("Swap two actual cells");
                }
                if (i < j && cells[i].is_actual()) {
                    // Актуальня ячейка переносится только в начало
                    throw std::runtime_error("Wrong swap #1");
                }
                if (i > j && cells[i].is_undefined()) {
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
    void set_identical_mapping_partial(Storage& cells, int from, int to) const {
        for (int ic = from; ic < to; ++ic) {
            cells[ic].geom().index = ic;
        }
    }

    /// @brief Устанавливает следующую позицию для части ячеек
    /// @param cells Ссылка на хранилище ячеек
    /// @param from, to Диапазон индексов в массивах actual_cells и undefined_cells
    void set_swap_mapping_partial(Storage& cells, int from, int to) const {
        for (int i = from; i < to; ++i) {
            int ai = actual_cells[i];
            int ui = undefined_cells[i];

            cells[ui].geom().index = ai;
            cells[ai].geom().index = ui;
        }
    }

    /// @brief Устанавливает следующие позиции ячеек в однопоточном режиме
    /// @details После выполнения операции поле element.index у ячеек указывает
    /// на следующее положение ячейки в хранилище
    void set_mapping(Storage& cells) const {
        set_identical_mapping_partial(cells, 0, cells.size());
        set_swap_mapping_partial(cells, 0, size());
#if SCRUTINY
        check_mapping(cells);
#endif
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Устанавливает следующие позиции ячеек в многопоточном режиме
    /// @details После выполнения операции поле element.index у ячеек указывает
    /// на следующее положение ячейки в хранилище
    void set_mapping(Storage& cells, ThreadPool& threads) const {
        auto num_tasks = threads.size();
        if (num_tasks < 2) {
            set_mapping(cells);
            return;
        }
        std::vector<std::future<void>> results(num_tasks);

        int bin = cells.size() / num_tasks + 1;
        int pos = 0;
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
    /// в поле element.index для части ячеек.
    /// @details Данные актуальной ячейки перемещаются на место неактуальной
    /// ячейки, индексы смежности не изменяются, они должны быть выставлены
    /// заранее.
    /// @param cells Ссылка на хранилище ячеек
    /// @param from, to Диапазон индексов в массиве undefined_cells
    void move_elements_partial(Storage &cells, int from, int to) const {
        for (int i = from; i < to; ++i) {
            int ic = undefined_cells[i];
            int jc = cells[ic].geom().index;

            cells[jc].copy_to(cells[ic]);

            cells[ic].geom().index = ic;
            cells[ic].geom().next = ic;
            cells[jc].set_undefined();
            cells[jc].geom().faces.set_undefined();
            cells[jc].geom().vertices.set_undefined();
        }
    }

    /// @brief Выполняет перестановку элементов в соответствии с индексом
    /// в поле element.index в однопоточном режиме.
    /// @details Данные актуальной ячейки перемещаются на место неактуальной
    /// ячейки, индексы смежности не изменяются, они должны быть выставлены
    /// заранее.
    void move_elements(Storage &cells) const {
        move_elements_partial(cells, 0, size());
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Выполняет перестановку элементов в соответствии с индексом
    /// в поле element.index в многопоточном режиме.
    /// @details Данные актуальной ячейки перемещаются на место неактуальной
    /// ячейки, индексы смежности не изменяются, они должны быть выставлены
    /// заранее.
    void move_elements(Storage &cells, ThreadPool& threads) const {
        auto num_tasks = threads.size();
        if (num_tasks < 2) {
            move_elements(cells);
            return;
        }
        std::vector<std::future<void>> results(num_tasks);

        std::int bin = size() / num_tasks + 1;
        std::int pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(
                    &SwapLists::move_elements_partial,
                    this, std::ref(cells),
                    pos, std::min(pos + bin, size())
            );
            pos += bin;
        }

        for (auto &result: results) result.get();
    }
#endif

    /// @brief Число ячеек для перестановки
    inline int size() const {
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
/// 2. Далее для ячеек изменяется поле element.index, которое после выполнения
/// функции set_mapping равно индексу, куда необходимо сместить ячейку.
/// 3. Выполняется функция change_adjacent, которая меняет индексы смежности
/// в соответствии с element.index
/// 4. Актуальные ячейки перемещаются на место неопределенных.
/// 5. Хранилище меняет размер, все неопределенные ячейки остаются за пределами
/// хранилища.
template<int dim>
void remove_undefined(Storage &cells, const Statistics &count) {
    static Stopwatch create_swap_timer;
    static Stopwatch set_mapping_timer;
    static Stopwatch change_adjacent_timer;
    static Stopwatch move_elements_timer;

    if (count.n_cells_large == count.n_cells_short) {
        return;
    }

    create_swap_timer.resume();
    int max_swap_count = count.n_cells_large - count.n_cells_short;
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

    cells.resize(count.n_cells_short);

#if CHECK_PERFORMANCE
    std::cout << "      Create SwapList elapsed: " << create_swap_timer.seconds() << " sec\n";
    std::cout << "      Set mapping elapsed: " << set_mapping_timer.seconds() << " sec\n";
    std::cout << "      Change adjacent elapsed: " << change_adjacent_timer.seconds() << " sec\n";
    std::cout << "      Move elements elapsed: " << move_elements_timer.seconds() << " sec\n";
#endif
}

} // namespace amr
} // namespace mesh
} // namespace zephyr