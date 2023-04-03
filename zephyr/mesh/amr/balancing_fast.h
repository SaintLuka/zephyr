/// @file Файл содержит реализацию быстрой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.
/// Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_restrictions.h>

namespace zephyr { namespace mesh { namespace amr {

/// @struct Ячейки из некоторого диапазона, распределенные по уровням и флагам адаптации
/// @brief Структура содержит индексы ячеек из диапазона в хранилище, которые
/// имеют флаг адаптации -1 или 0 (огрубляются или сохраняются)
/// @details Ячейки с максимальным уровнем адаптации не хранятся, поскольку этого не
/// требует алгоритм, также не хранятся индексы ячеек с флагом адаптации 1.
/// Фактически структура осуществляет соритровку ячеек по уровням (блочная сортировка,
/// работает за линейное время).
struct CellsByLevelPartial {
    Storage &cells;   ///< Ссылка на хранилище
    size_t from, to;  ///< Диапазон ячеек в хранилище

    std::vector<size_t> n_coarse; ///< Число ячеек каждого уровня, которые хотят огрубиться
    std::vector<size_t> n_retain; ///< Число ячеек каждого уровня, которые хотят сохраниться

    /// @brief Списки ячеек по уровням, которые хотят огрубиться
    std::vector<std::vector<size_t>> coarse;

    /// @brief Списки ячеек по уровням, которые хотят сохраниться
    std::vector<std::vector<size_t>> retain;

    /// @brief Создание экземпляра класса
    static CellsByLevelPartial create(Storage& cells, unsigned int max_level, size_t from, size_t to) {
        return CellsByLevelPartial(cells, max_level, from, to);
    }

    /// @brief Конструктор класса
    explicit CellsByLevelPartial(Storage& cells, unsigned int max_level, size_t from, size_t to)
    : cells(cells), from(from), to(to) {
        set_count(max_level);
        sort_by_levels(max_level);
    }

    /// @brief Посчитать ячейки каждого типа в диапазоне
    void set_count(unsigned int max_level) {
        n_coarse.resize(max_level);
        n_retain.resize(max_level);
        for (size_t ic = from; ic < to; ++ic) {
            auto lvl = cells[ic].level();
            auto flag = cells[ic].flag();

            if (lvl < max_level && flag < 1) {
                if (flag < 0) {
                    ++n_coarse[lvl];
                } else {
                    ++n_retain[lvl];
                }
            }
        }
    }

    /// @brief Распределить индексы ячеек в списках по уровням
    void sort_by_levels(unsigned int max_level) {
        retain.resize(max_level);
        coarse.resize(max_level);
        for (unsigned int lvl = 0; lvl < max_level; ++lvl) {
            coarse[lvl].reserve(n_coarse[lvl]);
            retain[lvl].reserve(n_retain[lvl]);
        }

        for (size_t ic = from; ic < to; ++ic) {
            auto lvl = cells[ic].level();
            auto flag = cells[ic].flag();
            if (lvl < max_level && flag < 1) {
                if (flag < 0) {
                    coarse[lvl].push_back(ic);
                } else {
                    retain[lvl].push_back(ic);
                }
            }
        }
    }
};

/// @struct Ячейки из хранилища, распределенные по уровням и флагам адаптации
/// Структура по содержанию и назначению соответствует структуре CellsByLevelPartial,
/// только характеризует ячейки со всего хранилища. В многопоточном режиме строится
/// на основе множества структур CellsByLevelPartial.
struct CellsByLevel {
    Storage &cells;   ///< Ссылка на хранилище

    std::vector<size_t> n_coarse; ///< Число ячеек каждого уровня, которые хотят огрубиться
    std::vector<size_t> n_retain; ///< Число ячеек каждого уровня, которые хотят сохраниться

    /// @brief Списки ячеек по уровням, которые хотят огрубиться
    std::vector<std::vector<size_t>> coarse;
    /// @brief Списки ячеек по уровням, которые хотят сохраниться
    std::vector<std::vector<size_t>> retain;

    /// @brief Однопоточный конструктор класса
    /// @details Вызывается конструктор CellsByLevelPartial для всего диапазона
    /// ячеек хранилища, затем данные перемещаются
    CellsByLevel(Storage &cells, unsigned int max_level)
            : cells(cells) {
        serial_constructor(max_level);
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Многопоточный конструктор класса
    /// @details Каждый поток вызывает конструктор CellsByLevelPartial для
    /// части ячеек, затем полученные данные складываются в один массив
    CellsByLevel(Storage &cells, unsigned int max_level, ThreadPool& threads)
            : cells(cells) {
        if (threads.size() < 2) {
            serial_constructor(max_level);
            return ;
        }
        auto num_tasks = threads.size();
        std::vector<std::future<CellsByLevelPartial>> parts(num_tasks);

        std::size_t bin = cells.size() / num_tasks + 1;
        std::size_t pos = 0;
        for (auto &count: parts) {
            count = threads.enqueue(&CellsByLevelPartial::create,
                                    std::ref(cells), max_level,
                                    pos, std::min(pos + bin, cells.size())
            );
            pos += bin;
        }

        n_coarse.resize(max_level, 0);
        n_retain.resize(max_level, 0);
        coarse.resize(max_level);
        retain.resize(max_level);

        for (auto &result: parts) {
            auto part = result.get();
            for (unsigned int lvl = 0; lvl < max_level; ++lvl) {
                n_retain[lvl] += part.n_retain[lvl];
                n_coarse[lvl] += part.n_coarse[lvl];

                coarse[lvl].insert(coarse[lvl].end(), part.coarse[lvl].begin(), part.coarse[lvl].end());
                retain[lvl].insert(retain[lvl].end(), part.retain[lvl].begin(), part.retain[lvl].end());
            }
        }
    }
#endif

private:
    /// @brief Однопоточный конструктор класса
    void serial_constructor(unsigned int max_level) {
        CellsByLevelPartial part(cells, max_level, 0, cells.size());
        n_coarse = std::move(part.n_coarse);
        n_retain = std::move(part.n_retain);
        coarse = std::move(part.coarse);
        retain = std::move(part.retain);
    }
};

/// @brief Обновляет флаг ячейки, которая имеет флаг 0.
/// Повышает флаг адаптации, если один из соседей хочет уровень на два выше
inline void retain_update_flag(Storage &cells, Storage::Item cell) {
    scrutiny_check(cell.flag() == 0, "retain_update_flag: cell.flag != 0")

    for (const auto &face: cell.geom().faces.list()) {
        if (face.is_undefined() or face.is_boundary()) {
            continue;
        }

        size_t neib_idx = face.adjacent.index;
        scrutiny_check(neib_idx < cells.size(), "retain_update_flag: neib_idx >= cells.size()")

        auto neib = cells[neib_idx];
        int neib_wanted = neib.level() + neib.flag();
        if (neib_wanted > cell.level() + 1) {
            cell.geom().flag = 1;
            return;
        }
    }
}

/// @brief Обновляет флаг ячейки, которая имеет флаг -1.
/// Ставит флаг адаптации 0, если сосед хочет уровень на 1 выше,
/// ставит флаг адаптации 1, если сосед хочет уровень на 2 выше,
/// флаги сиблингов не рассматриваются.
inline void coarse_update_flag(Storage &cells, Storage::Item cell) {
    scrutiny_check(cell.flag() < 0, "coarse_update_flag: cell.flag != -1")

    for (auto &face: cell.geom().faces.list()) {
        if (face.is_undefined() or face.is_boundary()) {
            continue;
        }

        auto neib_idx = face.adjacent.index;
        scrutiny_check(neib_idx < cells.size(), "coarse_update_flag: neib_idx >= cells.size()")

        auto neib = cells[neib_idx];
        int neib_wanted = neib.level() + neib.flag();

        if (neib_wanted > cell.level()) {
            if (neib_wanted > cell.level() + 1) {
                cell.geom().flag = 1;
                return;
            } else {
                cell.geom().flag = 0;
            }
        }
    }
}

/// @brief Обход диапазона ячеек с флагом = 0, данные ячейки могут только повысить
/// свой уровень до 1 за счет высокоуровневых соседей, после выпонения обхода
/// уровень ячеек больше не меняется
/// @param cells Ссылка на хранилище
/// @param indices Индексы ячеек с флагом = 0 в хранилище
/// @param from, to Диапазон индексов внутри массива indices
void round_1(Storage &cells, const std::vector<size_t> &indices, size_t from, size_t to) {
    for (size_t i = from; i < to; ++i) {
        scrutiny_check(indices[i] < cells.size(), "round_1: indices[i] >= cells.size()")

        auto cell = cells[indices[i]];

        scrutiny_check(cell.flag() == 0, "round_1: flag != 0")

        retain_update_flag(cells, cell);
    }
}

/// @brief Обход диапазона ячеек с флагом = -1, данные ячейки могут повысить
/// свой уровень до 0 или 1 за счет высокоуровневых соседей, флаги сиблингов
/// не используются. После выпонения обхода флаги ячеек ещё могут измениться
/// (@see round_3, @see round_4), но только у тех ячеек, которые не изменили
/// флаги на данном обходе (сохранили флаг -1).
/// @param cells Ссылка на хранилище
/// @param indices Индексы ячеек с флагом = -1 в хранилище
/// @param from, to Диапазон индексов внутри массива indices
void round_2(Storage &cells, const std::vector<size_t> &indices, size_t from, size_t to) {
    for (size_t i = from; i < to; ++i) {
        scrutiny_check(indices[i] < cells.size(), "round_2: indices[i] >= cells.size()")

        auto cell = cells[indices[i]];

        scrutiny_check(cell.flag() < 0, "round_2: flag != -1")

        coarse_update_flag(cells, cell);
    }
}

/// @brief Обход диапазона ячеек с изначальным флагом = -1, данные ячейки могут
/// повысить свой флаг только до 0 за счет соседей своего же уровня (если при
/// выполнении round_2 у кого-то из соседей появился флаг 1). После выполнения
/// обхода флаги ячеек ещё могут измениться (@see round_4).
/// @param cells Ссылка на хранилище
/// @param indices Индексы ячеек с исходным флагом = -1 в хранилище
/// @param from, to Диапазон индексов внутри массива indices
void round_3(Storage &cells, const std::vector<size_t> &indices, size_t from, size_t to) {
    for (size_t i = from; i < to; ++i) {
        scrutiny_check(indices[i] < cells.size(), "round_3: indices[i] >= cells.size()")

        auto cell = cells[indices[i]];
        if (cell.flag() < 0) {
            coarse_update_flag(cells, cell);
        }
    }
}

/// @brief Заключительный обход диапазона ячеек, у которых после базовых
/// ограничений был установлен флаг = -1, данные ячейки могут повысить свой
/// уровень до 0, если есть сиблинги, которые не хотят огрубляться, флаги
/// соседей уже не используются, после выпонения обхода флаги ячеек больше
/// не меняются
/// @param cells Ссылка на хранилище
/// @param indices Индексы ячеек с исходным флагом = -1 в хранилище
/// @param from, to Диапазон индексов внутри массива indices
template <unsigned int dim>
void round_4(Storage &cells, const std::vector<size_t> &indices, size_t from, size_t to) {
    for (size_t i = from; i < to; ++i) {
        scrutiny_check(indices[i] < cells.size(), "round_4: indices[i] >= cells.size()")

        size_t ic = indices[i];
        auto cell = cells[ic];

        if (cell.flag() >= 0) {
            continue;
        }

        auto sibs = get_siblings<dim>(cells, ic);
        for (auto is: sibs) {
            scrutiny_check(is < cells.size(), "round_4: Sibling index out of range")

            auto sib = cells[is];
            if (sib.flag() >= 0) {
                cell.geom().flag = 0;
                break;
            }
        }
    }
}

/// @brief Быстрая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище в соответствии с общими
/// ограничениями (ячейка нижнего уровня не огрубляется, ячейка верхнего -
/// не разбивается), а также с сохранением баланса 1:2 у соседних ячеек.
/// Баланс достигается повышением флагов адаптации (-1 и 0) у части ячеек.
/// Функция выполняется в однопоточном режиме.
/// @param cells Ссылка на хранилище ячеек
/// @param max_level Максимальный уровень ячеек
/// @details Детали алгоритма. На первом этапе на флаги адаптации накладываются
/// базовые ограничения. На втором этапе ячейки сортируются по уровням (только
/// те ячейки, которые имеют флаг 0 или -1 и уровень меньше максимального).
/// На третьем этапе выполняется цикл по уровням (начиная с максимального уровня
/// до базового). На каждом уровне:
/// 1. Сначала обходятся ячейки, которые хотят сохранить свой уровень
/// (@see round_1), их флаги могут быть повышены до 1.
/// 2. Дважды выполняется обход по ячейкам, которые хотят огрубиться
/// (@see round_2, @see round_3). В round_2 флаги могут повыситься до 0 и 1,
/// в round_3 флаги могут повыситься только до 0.
/// 4. На заключительном обходе координируются сиблинги, которые хотят
/// огрубиться (@see round_4), флаги могут повыситься до 0.
template <unsigned int dim = 1234>
void balance_flags_fast(Storage& cells, unsigned int max_level) {
    static Stopwatch restriction_timer;
    static Stopwatch sorting_timer;
    static Stopwatch round_timer_1;
    static Stopwatch round_timer_2;
    static Stopwatch round_timer_3;
    static Stopwatch round_timer_4;

    restriction_timer.resume();
    base_restrictions<dim>(cells, max_level);
    restriction_timer.stop();

    sorting_timer.resume();
    CellsByLevel sorted(cells, max_level);
    sorting_timer.stop();

    for (int lvl = int(max_level) - 1; lvl >= 0; --lvl) {
        round_timer_1.resume();
        round_1(cells, sorted.retain[lvl], 0, sorted.n_retain[lvl]);
        round_timer_1.stop();

        round_timer_2.resume();
        round_2(cells, sorted.coarse[lvl], 0, sorted.n_coarse[lvl]);
        round_timer_2.stop();

        round_timer_3.resume();
        round_3(cells, sorted.coarse[lvl], 0, sorted.n_coarse[lvl]);
        round_timer_3.stop();

        round_timer_4.resume();
        round_4<dim>(cells, sorted.coarse[lvl], 0, sorted.n_coarse[lvl]);
        round_timer_4.stop();
    }

#if CHECK_PERFORMANCE
    std::cout << "    Restriction elapsed: " << restriction_timer.seconds() << " sec\n";
    std::cout << "    Sorting elapsed: " << sorting_timer.seconds() << " sec\n";
    std::cout << "    Round 1 elapsed: " << round_timer_1.seconds() << " sec\n";
    std::cout << "    Round 2 elapsed: " << round_timer_2.seconds() << " sec\n";
    std::cout << "    Round 3 elapsed: " << round_timer_3.seconds() << " sec\n";
    std::cout << "    Round 4 elapsed: " << round_timer_4.seconds() << " sec\n";
#endif
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
template <>
void balance_flags_fast<1234>(Storage& cells, unsigned int max_level) {
    if (cells.empty())
        return;

    auto dim = cells[0].dim();

    if (dim < 3) {
        amr::balance_flags_fast<2>(cells, max_level);
    }
    else {
        amr::balance_flags_fast<3>(cells, max_level);
    }

#if SCRUTINY
    throw std::runtime_error("add check #2");
    // amr::check_flags(cells, max_level);
#endif
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Быстрая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище в соответствии с общими
/// ограничениями (ячейка нижнего уровня не огрубляется, ячейка верхнего -
/// не разбивается), а также с сохранением баланса 1:2 у соседних ячеек.
/// Баланс достигается повышением флагов адаптации (-1 и 0) у части ячеек.
/// Функция выполняется в многопоточном режиме.
/// @param cells Ссылка на хранилище ячеек
/// @param max_level Максимальный уровень ячеек
/// @param threads Ссылка на пул тредов
/// @details Детали алгоритма описаны у однопоточной весрии.
/// Многопоточная фунция быстрой балансировки может уступать по
/// производительности обычной итерационной функции балансировки.
template <unsigned int dim = 1234>
void balance_flags_fast(Storage& cells, unsigned int max_level, ThreadPool& threads) {
    using zephyr::performance::timer::Stopwatch;
    static Stopwatch restriction_timer;
    static Stopwatch sorting_timer;
    static Stopwatch round_timer_1;
    static Stopwatch round_timer_2;
    static Stopwatch round_timer_3;
    static Stopwatch round_timer_4;

    std::size_t num_tasks = threads.size();
    if (num_tasks < 2) {
        // В однопоточном режиме
        balance_flags_fast<dim>(cells, max_level);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    restriction_timer.resume();
    base_restrictions<dim>(cells, max_level, threads);
    restriction_timer.stop();

    sorting_timer.resume();
    CellsByLevel sorted(cells, max_level, threads);
    sorting_timer.stop();

    for (int lvl = int(max_level) - 1; lvl >= 0; --lvl) {
        round_timer_1.resume();
        std::size_t bin = sorted.n_retain[lvl] / num_tasks + 1;
        std::size_t pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(round_1,
                                  std::ref(cells), std::ref(sorted.retain[lvl]),
                                  pos, std::min(pos + bin, sorted.n_retain[lvl])
            );
            pos += bin;
        }
        for (auto &result: results) result.get();
        round_timer_1.stop();

        round_timer_2.resume();
        bin = sorted.n_coarse[lvl] / num_tasks + 1;
        pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(round_2,
                                  std::ref(cells), std::ref(sorted.coarse[lvl]),
                                  pos, std::min(pos + bin, sorted.n_coarse[lvl])
            );
            pos += bin;
        }
        for (auto &result: results) result.get();
        round_timer_2.stop();

        round_timer_3.resume();
        pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(round_3,
                                  std::ref(cells), std::ref(sorted.coarse[lvl]),
                                  pos, std::min(pos + bin, sorted.n_coarse[lvl])
            );
            pos += bin;
        }
        round_timer_3.stop();

        round_timer_4.resume();
        pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(round_4<dim>,
                                  std::ref(cells), std::ref(sorted.coarse[lvl]),
                                  pos, std::min(pos + bin, sorted.n_coarse[lvl])
            );
            pos += bin;
        }
        for (auto &result: results) result.get();
        round_timer_4.stop();
    }

#if CHECK_PERFORMANCE
    std::cout << "    Restriction elapsed: " << restriction_timer.times().wall() << "\n";
    std::cout << "    Sorting elapsed: " << sorting_timer.times().wall() << "\n";
    std::cout << "    Round 1 elapsed: " << round_timer_1.times().wall() << "\n";
    std::cout << "    Round 2 elapsed: " << round_timer_2.times().wall() << "\n";
    std::cout << "    Round 3 elapsed: " << round_timer_3.times().wall() << "\n";
    std::cout << "    Round 4 elapsed: " << round_timer_4.times().wall() << "\n";
#endif
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
template <>
void balance_flags_fast<1234>(Storage& cells, unsigned int max_level, ThreadPool& threads) {
    if (cells.empty())
        return;

    auto dim = cells[0][element].dimension;

    if (dim < 3) {
        amr::balance_flags_fast<2>(cells, max_level, threads);
    }
    else {
        amr::balance_flags_fast<3>(cells, max_level, threads);
    }

#if SCRUTINY
    amr::check_flags(cells, max_level);
#endif
}
#endif

} // namespace amr
} // namespace mesh
} // namespace zephyr