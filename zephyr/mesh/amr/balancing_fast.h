/// @file Файл содержит реализацию быстрой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.
/// Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <numeric>
#include <iomanip>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_restrictions.h>

namespace zephyr::mesh::amr {

/// @brief Ячейки из некоторого диапазона, распределенные по уровням и флагам адаптации
/// Структура содержит индексы ячеек из диапазона в хранилище, которые имеют флаг
/// адаптации -1 или 0 (огрубляются или сохраняются)
/// @details Ячейки с максимальным уровнем адаптации не хранятся, поскольку этого не
/// требует алгоритм, также не хранятся индексы ячеек с флагом адаптации 1.
/// Фактически структура осуществляет соритровку ячеек по уровням (блочная сортировка,
/// работает за линейное время).
struct CellsByLevelPartial {
    AmrCells &cells;   ///< Ссылка на хранилище
    index_t from, to;  ///< Диапазон ячеек в хранилище

    std::vector<index_t> n_coarse; ///< Число ячеек каждого уровня, которые хотят огрубиться
    std::vector<index_t> n_retain; ///< Число ячеек каждого уровня, которые хотят сохраниться

    /// @brief Списки ячеек по уровням, которые хотят огрубиться
    std::vector<std::vector<index_t>> coarse;

    /// @brief Списки ячеек по уровням, которые хотят сохраниться
    std::vector<std::vector<index_t>> retain;

    /// @brief Создание экземпляра класса
    static CellsByLevelPartial create(AmrCells& cells, int max_level, index_t from, index_t to) {
        return CellsByLevelPartial(cells, max_level, from, to);
    }

    /// @brief Конструктор класса
    explicit CellsByLevelPartial(AmrCells& cells, int max_level, index_t from, index_t to)
    : cells(cells), from(from), to(to) {
        set_count(max_level);
        sort_by_levels(max_level);
    }

    /// @brief Посчитать ячейки каждого типа в диапазоне
    void set_count(int max_level) {
        n_coarse.resize(max_level);
        n_retain.resize(max_level);
        for (index_t ic = from; ic < to; ++ic) {
            auto lvl = cells.level[ic];
            auto flag = cells.flag[ic];

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
    void sort_by_levels(int max_level) {
        retain.resize(max_level);
        coarse.resize(max_level);
        for (int lvl = 0; lvl < max_level; ++lvl) {
            coarse[lvl].reserve(n_coarse[lvl]);
            retain[lvl].reserve(n_retain[lvl]);
        }

        for (index_t ic = from; ic < to; ++ic) {
            auto lvl = cells.level[ic];
            auto flag = cells.flag[ic];
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

/// @brief Ячейки из хранилища, распределенные по уровням и флагам адаптации
/// Структура по содержанию и назначению соответствует структуре CellsByLevelPartial,
/// только характеризует ячейки со всего хранилища. В многопоточном режиме строится
/// на основе множества структур CellsByLevelPartial.
struct CellsByLevel {
    AmrCells &cells;   ///< Ссылка на хранилище

    std::vector<index_t> n_coarse; ///< Число ячеек каждого уровня, которые хотят огрубиться
    std::vector<index_t> n_retain; ///< Число ячеек каждого уровня, которые хотят сохраниться

    /// @brief Списки ячеек по уровням, которые хотят огрубиться
    std::vector<std::vector<index_t>> coarse;

    /// @brief Списки ячеек по уровням, которые хотят сохраниться
    std::vector<std::vector<index_t>> retain;

    /// @brief Однопоточный конструктор класса
    /// @details Вызывается конструктор CellsByLevelPartial для всего диапазона
    /// ячеек хранилища, затем данные перемещаются
    CellsByLevel(AmrCells &cells, int max_level)
            : cells(cells) {
        serial_constructor(max_level);
    }

#ifdef ZEPHYR_MULTITHREADING
    /// @brief Многопоточный конструктор класса
    /// @details Каждый поток вызывает конструктор CellsByLevelPartial для
    /// части ячеек, затем полученные данные складываются в один массив
    CellsByLevel(AmrStorage &cells, int max_level, ThreadPool& threads)
            : cells(cells) {
        if (threads.size() < 2) {
            serial_constructor(max_level);
            return ;
        }
        auto num_tasks = threads.size();
        std::vector<std::future<CellsByLevelPartial>> parts(num_tasks);

        std::index_t bin = cells.size() / num_tasks + 1;
        std::index_t pos = 0;
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
            for (int lvl = 0; lvl < max_level; ++lvl) {
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
    void serial_constructor(int max_level) {
        CellsByLevelPartial part(cells, max_level, 0, cells.size());
        n_coarse = std::move(part.n_coarse);
        n_retain = std::move(part.n_retain);
        coarse = std::move(part.coarse);
        retain = std::move(part.retain);
    }
};

/// @brief Обновляет флаг ячейки под индексом index, которая имеет флаг 0.
/// Повышает флаг адаптации, если один из соседей хочет уровень на два выше
inline void retain_update_flag(int ic, AmrCells &cells) {
    scrutiny_check(ic < cells.size(), "round_1: cell_idx >= cells.size()")
    scrutiny_check(cells.flag[ic] == 0, "retain_update_flag: cell.flag != 0")

    for (auto iface: cells.faces_range(ic)) {
        if (cells.faces.is_undefined(iface) || cells.faces.is_boundary(iface)) {
            continue;
        }

        int jc = cells.faces.adjacent.index[iface];
        scrutiny_check(jc < cells.size(), "retain_update_flag: neib_idx >= cells.n_cells()")

        int neib_wanted = cells.level[jc] + cells.flag[jc];
        if (neib_wanted > cells.level[ic] + 1) {
            cells.flag[ic] = 1;
            return;
        }
    }
}

/// @brief Обновляет флаг ячейки, которая имеет флаг -1.
/// Ставит флаг адаптации 0, если сосед хочет уровень на 1 выше,
/// ставит флаг адаптации 1, если сосед хочет уровень на 2 выше,
/// флаги сиблингов не рассматриваются.
inline void coarse_update_flag(int ic, AmrCells &cells) {
    scrutiny_check(ic < cells.size(), "round_2/3: cell_idx >= cells.size()")

    if (cells.flag[ic] >= 0) {
        return;
    }

    for (auto iface: cells.faces_range(ic)) {
        if (cells.faces.is_undefined(iface) || cells.faces.is_boundary(iface)) {
            continue;
        }

        index_t jc = cells.faces.adjacent.index[iface];
        scrutiny_check(jc < cells.size(), "coarse_update_flag: neib_idx >= cells.n_cells()")

        int neib_wanted = cells.level[jc] + cells.flag[jc];

        if (neib_wanted > cells.level[ic]) {
            if (neib_wanted > cells.level[ic] + 1) {
                cells.flag[ic] = 1;
                return;
            } else {
                cells.flag[ic] = 0;
            }
        }
    }
}

/// @brief Обновляет флаг ячейки, которая имеет флаг -1.
/// Ставит флаг адаптации 0, если один из сиблингов не хочет огрубляться.
template <int dim>
inline void coarse_update_flag_by_sibs(int ic, AmrCells &cells) {
    scrutiny_check(ic < cells.size(), "round_4: cell_idx >= cells.size()")

    if (cells.flag[ic] >= 0) {
        return;
    }

    auto sibs = get_siblings<dim>(cells, ic);
    for (auto is: sibs) {
        scrutiny_check(is < cells.size(), "round_4: Sibling index out of range")

        if (cells.flag[is] >= 0) {
            cells.flag[ic] = 0;
            break;
        }
    }
}

/// @brief Обход диапазона ячеек с флагом = 0, данные ячейки могут только повысить
/// свой уровень до 1 за счет высокоуровневых соседей, после выпонения обхода
/// уровень ячеек больше не меняется
/// @param indices Индексы ячеек с флагом = 0 в хранилище
/// @param cells Ссылка на хранилище
void round_1(const std::vector<index_t> &indices, AmrCells &cells) {
    threads::for_each(
            indices.begin(), indices.end(),
            retain_update_flag, std::ref(cells));
}

/// @brief Обход диапазона ячеек с флагом = -1, данные ячейки могут повысить
/// свой уровень до 0 или 1 за счет высокоуровневых соседей, флаги сиблингов
/// не используются. После выпонения обхода флаги ячеек ещё могут измениться
/// (@see round_3, @see round_4), но только у тех ячеек, которые не изменили
/// флаги на данном обходе (сохранили флаг -1).
/// @param indices Индексы ячеек с флагом = -1 в хранилище
/// @param cells Ссылка на хранилище
void round_2(const std::vector<index_t> &indices, AmrCells &cells) {
    threads::for_each(
            indices.begin(), indices.end(),
            coarse_update_flag, std::ref(cells));
}

/// @brief Обход диапазона ячеек с изначальным флагом = -1, данные ячейки могут
/// повысить свой флаг только до 0 за счет соседей своего же уровня (если при
/// выполнении round_2 у кого-то из соседей появился флаг 1). После выполнения
/// обхода флаги ячеек ещё могут измениться (@see round_4).
/// @param indices Индексы ячеек с исходным флагом = -1 в хранилище
/// @param cells Ссылка на хранилище
void round_3(const std::vector<index_t> &indices, AmrCells &cells) {
    threads::for_each(
            indices.begin(), indices.end(),
            coarse_update_flag, std::ref(cells));
}

/// @brief Заключительный обход диапазона ячеек, у которых после базовых
/// ограничений был установлен флаг = -1, данные ячейки могут повысить свой
/// уровень до 0, если есть сиблинги, которые не хотят огрубляться, флаги
/// соседей уже не используются, после выпонения обхода флаги ячеек больше
/// не меняются
/// @param indices Индексы ячеек с исходным флагом = -1 в хранилище
/// @param cells Ссылка на хранилище
template <int dim>
void round_4(const std::vector<index_t> &indices, AmrCells &cells) {
    threads::for_each(
            indices.begin(), indices.end(),
            coarse_update_flag_by_sibs<dim>, std::ref(cells));
}

/// @brief Быстрая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище в соответствии с общими
/// ограничениями (ячейка нижнего уровня не огрубляется, ячейка верхнего -
/// не разбивается), а также с сохранением баланса 1:2 у соседних ячеек.
/// Баланс достигается повышением флагов адаптации (-1 и 0) у части ячеек.
/// @param cells Ссылка на хранилище ячеек
/// @param max_level Максимальный уровень адаптации ячеек
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
template <int dim>
void balance_flags_fast(AmrCells& cells, int max_level) {
    static Stopwatch restriction_timer;
    static Stopwatch sorting_timer;
    static Stopwatch round_timer_1;
    static Stopwatch round_timer_2;
    static Stopwatch round_timer_3;
    static Stopwatch round_timer_4;
    static int n_step = 0;
    static int n_total_retain = 0;
    static int n_total_coarse = 0;

    restriction_timer.resume();
    base_restrictions<dim>(cells, max_level);
    restriction_timer.stop();

    sorting_timer.resume();
    // TODO: Сортировка по тредам
    CellsByLevel sorted(cells, max_level);
    sorting_timer.stop();

    for (int lvl = max_level - 1; lvl >= 0; --lvl) {
        round_timer_1.resume();
        round_1(sorted.retain[lvl], cells);
        round_timer_1.stop();

        round_timer_2.resume();
        round_2(sorted.coarse[lvl], cells);
        round_timer_2.stop();

        round_timer_3.resume();
        round_3(sorted.coarse[lvl], cells);
        round_timer_3.stop();

        round_timer_4.resume();
        round_4<dim>(sorted.coarse[lvl], cells);
        round_timer_4.stop();
    }

#if CHECK_PERFORMANCE
    n_step += 1;
    n_total_coarse += std::accumulate(sorted.n_coarse.begin(), sorted.n_coarse.end(), 0);
    n_total_retain += std::accumulate(sorted.n_retain.begin(), sorted.n_retain.end(), 0);

    static index_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Base restrictions: " << std::setw(8) << restriction_timer.milliseconds() << " ms\n";
        mpi::cout << "    Sort by level:     " << std::setw(8) << sorting_timer.milliseconds() << " ms\n";
        mpi::cout << "    Round 1 elapsed:   " << std::setw(8) << round_timer_1.milliseconds() << " ms  ";
        mpi::cout << "    (avg retain: " << std::setw(10) << n_total_retain / n_step << ")\n";
        mpi::cout << "    Round 2 elapsed:   " << std::setw(8) << round_timer_2.milliseconds() << " ms  ";
        mpi::cout << "    (avg coarse: " << std::setw(10) << n_total_coarse / n_step << ")\n";
        mpi::cout << "    Round 3 elapsed:   " << std::setw(8) << round_timer_3.milliseconds() << " ms\n";
        mpi::cout << "    Round 4 elapsed:   " << std::setw(8) << round_timer_4.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

} // namespace zephyr::mesh::amr