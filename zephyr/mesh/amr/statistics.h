/// @file Файл содержит реализации структур StatisticsPartial и Statistics, которые
/// хранят информацию о количестве ячеек с различными флагами.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr { namespace mesh { namespace amr {

/// @struct Содержит статистику о числе ячеек с флагами -1, 0 и 1 в некотором
/// диапазоне индексов. Вспомогательная структура для реализации многопоточности.
struct StatisticsPartial {
    size_t n_coarse;  ///< Число ячеек в диапазоне для огрубления
    size_t n_retain;  ///< Число ячеек в диапазоне, которые сохраняют уровень
    size_t n_refine;  ///< Число ячеек в диапазоне для разбиения

    /// @brief Статическая функция, вызывает конструктор
    static StatisticsPartial create(Storage & cells, size_t from, size_t to) {
        return StatisticsPartial(cells, from, to);
    }

    /// @brief Базовый конструктор
    /// @param cells Хранилище ячеек
    /// @param from, to Диапазон ячеек, для которых собирается статистика
    StatisticsPartial(Storage &cells, size_t from, size_t to)
        : n_coarse(0), n_retain(0), n_refine(0) {

        for (size_t ic = from; ic < to; ++ic) {
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
    }
};

/// @struct Содержит статистику о числе ячеек с флагами -1, 0 и 1, число новых
/// дочерних и родительских ячеек и другие параметры
template<int dim>
struct Statistics {
    size_t n_cells;     ///< Исходное количество в Storage (перед адаптацией)
    size_t n_coarse;    ///< Число ячеек для огрубления (кратно числу дочерних ячеек CpC)
    size_t n_retain;    ///< Число ячеек для сохранения
    size_t n_refine;    ///< Число ячеек для разбиения
    size_t n_parents;   ///< Число новых родительских ячеек (n_coarse/CpC)
    size_t n_children;  ///< Число новых детей (CpC * n_refine)

    /// @brief Размер расширенного хранилища, после добавления новых ячеек,
    /// полученных разбиением или огрублением старых, всегда больше или равно
    /// исходному размеру хранилища n_cells
    size_t n_cells_large;

    /// @brief Итоговый размер хранилища, после удаления старых (не листовых)
    /// ячеек из хранилища
    size_t n_cells_short;

    /// @brief Конструктор, однопоточный сбор статистики
    Statistics(Storage &cells) {
        serial_constructor(cells);
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Конструктор, многопоточный сбор статистики.
    /// @details Данные получаются после объединения статистики с разных потоков
    Statistics(Storage& cells, ThreadPool& threads)
        : n_coarse(0), n_retain(0), n_refine(0) {
        n_cells = cells.size();

        auto num_tasks = threads.size();
        if (num_tasks < 2) {
            serial_constructor(cells);
            return ;
        }
        std::vector<std::future<StatisticsPartial>> results(num_tasks);
        std::size_t bin = n_cells / num_tasks + 1;
        std::size_t pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(StatisticsPartial::create,
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

        compute_dependent();
    }
#endif

    /// @brief Однопоточный конструктор класса
    void serial_constructor(Storage &cells) {
        n_cells = cells.size();
        StatisticsPartial sp(cells, 0, n_cells);

        n_coarse = sp.n_coarse;
        n_retain = sp.n_retain;
        n_refine = sp.n_refine;

        compute_dependent();
    }

    /// @brief Вычисляет значения зависимых величин, осуществляет проверки
    void compute_dependent() {
        if (n_retain + n_refine + n_coarse != n_cells) {
            throw std::runtime_error("Refiner::apply() error 1");
        }

        if (n_coarse % CpC(dim) != 0) {
            throw std::runtime_error("Refiner::apply() error 2");
        }

        n_parents = n_coarse / CpC(dim);
        n_children = n_refine * CpC(dim);

        n_cells_large = n_cells + n_parents + n_children;
        n_cells_short = n_retain + n_children + n_parents;
    }

    /// @brief Выводит статистику
    void print() const {
        std::cout << "Apply statistic:\n";
        std::cout << "\tn_cells:       " << n_cells << "\n";
        std::cout << "\tn_retain:      " << n_retain << "\n";
        std::cout << "\tn_coarse:      " << n_coarse << "\n";
        std::cout << "\tn_parents:     " << n_parents << "\n";
        std::cout << "\tn_refine:      " << n_refine << "\n";
        std::cout << "\tn_children:    " << n_children << "\n";
        std::cout << "\tn_cells_large: " << n_cells_large << "\n";
        std::cout << "\tn_cells_short: " << n_cells_short << "\n";
    }

#ifdef ZEPHYR_ENABLE_MPI
    /// @brief Выводит статистику в MPI версии
    void print(Network& network) const {
        for (int r = 0; r < network.size(); ++r) {
            if (network.rank() == r) {
                std::cout << "rank " << r << "\n";
                print();
            }
            MPI_Barrier(network.comm());
        }
    }
#endif
};

} // namespace amr
} // namespace mesh
} // namespace zephyr