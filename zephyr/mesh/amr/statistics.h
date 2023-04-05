/// @file Файл содержит реализации структур PartialStatistics и Statistics, которые
/// хранят информацию о количестве ячеек с различными флагами.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/mpi.h>

namespace zephyr { namespace mesh { namespace amr {

using zephyr::utils::threads;
using zephyr::utils::mpi;

/// @struct Содержит статистику о числе ячеек с флагами -1, 0 и 1 в некотором
/// диапазоне ячеек. Вспомогательная структура для реализации многопоточности.
struct PartialStatistics {
    int n_coarse;  ///< Число ячеек в диапазоне для огрубления
    int n_retain;  ///< Число ячеек в диапазоне, которые сохраняют уровень
    int n_refine;  ///< Число ячеек в диапазоне для разбиения

    /// @brief Простой конструктор
    PartialStatistics(int n_coarse, int n_retain, int n_refine)
        : n_coarse(n_coarse), n_retain(n_retain), n_refine(n_refine) {
    }
};

/// @brief Собрать частичную статистику с части хранилища
static PartialStatistics partial_statistics(Storage& cells, int from, int to) {
    int n_coarse = 0;
    int n_retain = 0;
    int n_refine = 0;

    for (int ic = from; ic < to; ++ic) {
        switch (cells[ic].flag()) {
            case 0:
                ++n_retain;
                break;
            case 1:
                ++n_refine;
                break;
            default:
                ++n_coarse;
                break;
        }
    }
    return { n_coarse, n_retain, n_refine };
}

/// @struct Содержит статистику о числе ячеек с флагами -1, 0 и 1, число новых
/// дочерних и родительских ячеек и другие параметры
struct Statistics {
    int n_cells;     ///< Исходное количество в Storage (перед адаптацией)
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
    explicit Statistics(Storage &cells) :
        n_coarse(0), 
        n_retain(0), 
        n_refine(0) {

        n_cells = cells.size();

        if (n_cells < 1) {
            throw std::runtime_error("Empty Storage statistics");
        }

        // Составим статистику
        if (threads::is_off()) {
            // В один поток
            PartialStatistics sp = partial_statistics(cells, 0, n_cells);

            n_coarse = sp.n_coarse;
            n_retain = sp.n_retain;
            n_refine = sp.n_refine;
        }
        else {
            // В много потоков
            auto n_tasks = threads::count();

            std::vector<std::future<PartialStatistics>> results(n_tasks);
            int bin = n_cells / n_tasks + 1;
            int pos = 0;
            for (auto &res : results) {
                res = std::async(partial_statistics,
                                 std::ref(cells),
                                 pos, std::min(pos + bin, n_cells)
                );
                pos += bin;
            }

            for (auto &result: results) {
                auto count = result.get();
                n_coarse += count.n_coarse;
                n_retain += count.n_retain;
                n_refine += count.n_refine;
            }
        }

        // Пересчитаем завимисые величины

        int dim = cells[0].dim();

        if (n_retain + n_refine + n_coarse != n_cells) {
            throw std::runtime_error("Refiner::apply() error 1");
        }

        if (n_coarse % CpC(dim) != 0) {
            throw std::runtime_error("Refiner::apply() error 2");
        }

        n_parents =  n_coarse / CpC(dim);
        n_children = n_refine * CpC(dim);

        n_cells_large = n_cells + n_parents + n_children;
        n_cells_short = n_retain + n_children + n_parents;
    }

    /// @brief Выводит статистику
    void print() const {
        mpi::for_each([this]() {
            if (!mpi::is_single()) {
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

} // namespace amr
} // namespace mesh
} // namespace zephyr