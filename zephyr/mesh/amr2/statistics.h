/// @file Файл содержит реализации структур PartialStatistics и Statistics, которые
/// хранят информацию о количестве ячеек с различными флагами.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh::amr2 {

using zephyr::utils::threads;
using zephyr::utils::mpi;

/// @brief Содержит статистику о числе ячеек с флагами -1, 0 и 1 в некотором
/// диапазоне ячеек. В том числе и для одной ячейки.
struct PartStatistics {
    int n_coarse = 0; ///< Число ячеек в диапазоне для огрубления
    int n_retain = 0; ///< Число ячеек в диапазоне, которые сохраняют уровень
    int n_refine = 0; ///< Число ячеек в диапазоне для разбиения

    /// @brief Конструктор по умолчанию с нулями.
    PartStatistics() = default;

    /// @brief Для суммирования статистики в тредах
    PartStatistics& operator+=(const PartStatistics& other) {
        n_coarse += other.n_coarse;
        n_retain += other.n_retain;
        n_refine += other.n_refine;
        return *this;
    }
};

inline PartStatistics cell_statistics(AmrStorage::Item& cell) {
    if (!cell.flag) {
        return {.n_retain=1};
    } else {
        if (cell.flag > 0) {
            return {.n_refine=1};
        } else {
            return {.n_coarse=1};
        }
    }
}

/// @brief Содержит статистику о числе ячеек с флагами -1, 0 и 1, число новых
/// дочерних и родительских ячеек и другие параметры
struct Statistics {
    int n_cells;     ///< Исходное количество в AmrStorage (перед адаптацией)
    int n_coarse;    ///< Число ячеек для огрубления (кратно числу дочерних ячеек CpC)
    int n_retain;    ///< Число ячеек для сохранения
    int n_refine;    ///< Число ячеек для разбиения
    int n_parents;   ///< Число новых родительских ячеек (n_coarse/CpC)
    int n_children;  ///< Число новых детей (CpC * n_refine)

    /// @brief Размер расширенного хранилища, после добавления новых ячеек,
    /// полученных разбиением или огрублением старых, всегда больше или равно
    /// исходному размеру хранилища n_cells
    int n_cells_large;

    /// @brief Итоговый размер хранилища, после удаления старых (не листовых)
    /// ячеек из хранилища
    int n_cells_short;

    /// @brief Конструктор, однопоточный сбор статистики
    /// @details В многопоточном режиме данные получаются после объединения 
    /// статистики с разных потоков
    explicit Statistics(AmrStorage &cells) {
        n_cells = cells.size();

        scrutiny_check(n_cells >= 0, "Empty AmrStorage statistics");

        PartStatistics ps = threads::sum(
                cells.begin(), cells.end(), {}, cell_statistics);

        n_coarse = ps.n_coarse;
        n_retain = ps.n_retain;
        n_refine = ps.n_refine;

        // Пересчитаем завимисые величины
        int dim = cells[0].dim;

        scrutiny_check(n_coarse % CpC(dim) == 0, "Refiner::apply() error #2")
        scrutiny_check(n_retain + n_refine + n_coarse == n_cells, "Refiner::apply() error #1");

        n_parents =  n_coarse / CpC(dim);
        n_children = n_refine * CpC(dim);

        n_cells_large = n_cells + n_parents + n_children;
        n_cells_short = n_retain + n_children + n_parents;
    }

    /// @brief Выводит статистику
    void print() const {
        mpi::for_each([this]() {
            if (!mpi::single()) {
                std::cout << "rank " << mpi::rank() << "\n";
            }
            std::cout << "Apply statistic:\n";
            std::cout << "\tn_cells:       " << n_cells << "\n";
            std::cout << "\tn_retain:      " << n_retain << "\n";
            std::cout << "\tn_coarse:      " << n_coarse << "\n";
            std::cout << "\tn_parents:     " << n_parents << "\n";
            std::cout << "\tn_refine:      " << n_refine << "\n";
            std::cout << "\tn_children:    " << n_children << "\n";
            std::cout << "\tn_cells_large: " << n_cells_large << "\n";
            std::cout << "\tn_cells_short: " << n_cells_short << "\n";
        });
    }
};

} // namespace zephyr::mesh::amr2