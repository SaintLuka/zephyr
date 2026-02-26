// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr::mesh::amr {

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

} // 