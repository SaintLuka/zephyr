/// @file Файл содержит реализацию простой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для
/// разработчиков. Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_restrictions.h>

namespace zephyr { namespace mesh { namespace amr {

#ifdef ZEPHYR_ENABLE_MULTITHREADING
using ::zephyr::multithreading::dummy_pool;
#endif

/// @brief Простая структура, содержит ссылки на данные (уровни и флаги)
/// сосдених ячеек и сиблингов, требуется для быстрого доступа к информации
/// об обновлении соседей или сиблингов.
/// Для ячеек с флагом = 1 данные не хранятся (нет необходимости), для
/// ячеек с флагом = 0 не хранятся данные о сиблингах.
template<int dim>
struct CellsAround {
private:

    /// @brief Тип данных amr.flag
    using flag_type = int;

    /// @brief Максимальное число соседей
    static const int max_neibs_size = FpC(dim) * FpF(dim);

    /// @brief Максимальное число сиблингов
    static const int max_sibs_size = CpC(dim) - 1;

    int neib_count;                               ///< Действительное количество соседей
    int neib_levels[max_neibs_size];              ///< Уровни адаптации соседей
    const flag_type *neib_flags[max_neibs_size];  ///< Ссылки на флаги адаптации соседей

    int sibs_count;                               ///< Количество сиблингов (0 или max_sibs_size)
    const flag_type *sibs_flags[max_sibs_size];   ///< Ссылки на флаги адаптации сиблингов

public:

    int max_neib_wanted_level;                    ///< Максимальный желаемый уровень среди соседей
    int max_sibs_flag;                            ///< Максимальный флаг сиблинга


    /// @brief Конструктор по умолчанию - данные отсутствуют
    CellsAround() : neib_count(0), sibs_count(0) {};

    /// @brief Заполнение списков соседей для ячейки
    /// @param locals Ссылка на хранилище
    /// @param ic Индекс ячейки в хранилище
    void setup(Storage &locals, Storage& aliens, int ic) {
        scrutiny_check(ic < locals.size(), "setup: ic >= locals.size()")

        Cell& cell = locals[ic].geom();

        // Эти ячейки не интересуют
        if (cell.flag > 0) {
            neib_count = 0;
            sibs_count = 0;
            return;
        }

        // Поиск соседей
        neib_count = 0;
        for (auto &face: cell.faces) {
            if (face.is_undefined() || face.is_boundary()) {
                continue;
            }

            auto &adj = face.adjacent;
#if SCRUTINY
            if (adj.rank == cell.rank) {
                if (adj.index >= locals.size()) {
                    std::cerr << "rank: " << adj.rank << "\n";
                    std::cerr << "locals.size: " << locals.size() << "\n";
                    std::cerr << "adj.index: " << adj.index << "\n";
                    throw std::runtime_error("CellsAround::setup error #11");
                }
            } else {
                if (adj.ghost >= aliens.size()) {
                    throw std::runtime_error("CellsAround::setup error #2");
                }
            }
#endif
            Cell& neib = (adj.rank == cell.rank ? locals[adj.index] : aliens[adj.ghost]).geom();

            neib_levels[neib_count] = neib.level;
            neib_flags[neib_count] = &neib.flag;
            ++neib_count;
        }
        scrutiny_check(neib_count <= max_neibs_size, "nei_count > max_neibs_size")

        // Поиск сиблингов для ячеек с флагом -1
        if (cell.flag  < 0) {
            sibs_count = 0;
            if (cell.flag < 0) {
                scrutiny_check(can_coarse<dim>(locals, ic), "Can't coarse")

                // Код выполняется только если can_coarse было равно истинно,
                // поэтому можно не волноваться о начилии сиблингов
                auto sibs = get_siblings<dim>(locals, ic);
                for (auto is: sibs) {
                    scrutiny_check(is < locals.size(), "CellsAround: Sibling index out of range")

                    Cell& sib = locals[is].geom();
                    sibs_flags[sibs_count] = &sib.flag;
                    ++sibs_count;
                }

                scrutiny_check(sibs_count == max_sibs_size, "Wrong siblings count")
            }
        }
    }

    /// @brief Собрать актуальные данные о соседях и сиблингах. Функция
    /// инициализирует переменные max_neib_wanted_level и max_sibs_flag.
    void collect() {
        // Максимальный уровень, который хочет получить сосед после адаптации
        // Значение зависит от текущего состояния окружения
        max_neib_wanted_level = 0;
        for (int i = 0; i < neib_count; ++i) {
            scrutiny_check(neib_flags[i] != nullptr, "Undefined neighbor")
            max_neib_wanted_level = std::max(
                    neib_levels[i] + *neib_flags[i],
                    max_neib_wanted_level
            );
        }

        // Максимальный уровень, который хочет получить сосед после адаптации
        // Значение зависит от текущего состояния окружения
        max_sibs_flag = -1;
        for (int i = 0; i < sibs_count; ++i) {
            scrutiny_check(sibs_flags[i] != nullptr, "Undefined neighbor")
            max_sibs_flag = std::max(max_sibs_flag, int(*sibs_flags[i]));
        }
    }
};

/// @brief Список из структур типа CellsAround. Список может строится как в
/// однопоточном, так и в параллельном режиме.
/// Для ячеек с флагом = 1 данные не хранятся (нет необходимости), для
/// ячеек с флагом = 0 не хранятся данные о сиблингах.
template <int dim>
struct CellsAroundList {
    std::vector<CellsAround<dim>> m_list;

    /// @brief Однопоточный конструктор для построения окружения
    CellsAroundList(Storage& locals, Storage& aliens) {
        serial_constructor(locals, aliens);
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Многопоточный конструктор для построения окружения
    CellsAroundList(Storage& locals, Storage& aliens, ThreadPool& threads) {
        if (threads.size() < 2) {
            // Вызов последовательного кода
            serial_constructor(locals, aliens);
            return;
        }
        m_list.resize(locals.size());

        auto num_tasks = threads.size();
        std::vector<std::future<void>> results(num_tasks);

        std::int bin = locals.size() / num_tasks + 1;
        std::int pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(
                    [this](Storage& locals, Storage& aliens, int from, int to) {
                        for(int ic = from; ic < to; ++ic) {
                            m_list[ic].setup(locals, aliens, ic);
                        }
                    },
                    std::ref(locals), std::ref(aliens),
                    pos, std::min(pos + bin, locals.size())
            );
            pos += bin;
        }
        for (auto &result: results)
            result.get();
    }
#endif

    /// @brief Собрать актуальные данные о соседях и сиблингах для всех ячеек.
    /// Функция инициализирует max_neib_wanted_level и max_sibs_flag.
    void collect() {
        collect_partial(0, m_list.size());
    }

#ifdef ZEPHYR_ENABLE_MULTITHREADING
    /// @brief Собрать актуальные данные о соседях и сиблингах для всех ячеек
    /// в многопоточном режиме. Функция инициализирует переменные
    /// max_neib_wanted_level и max_sibs_flag.
    void collect(ThreadPool& threads) {
        auto num_tasks = threads.size();
        if (num_tasks < 2) {
            return collect();
        }
        std::vector<std::future<void>> results(num_tasks);

        std::int bin = m_list.size() / num_tasks + 1;
        std::int pos = 0;
        for (auto &res : results) {
            res = threads.enqueue(
                    [this](int from, int to) {
                        collect_partial(from, to);
                    },
                    pos, std::min(pos + bin, m_list.size())
            );
            pos += bin;
        }
        for (auto &result: results)
            result.get();
    }
#endif

    /// @brief Оператор доступа к элементу списка
    const CellsAround<dim>& operator[](int ic) const {
        scrutiny_check(ic < m_list.size(), "ic >= CellsAround.m_list.size()")
        return m_list[ic];
    }

private:
    /// @brief Однопоточный конструктор
    void serial_constructor(Storage& locals, Storage& aliens) {
        m_list.resize(locals.size());
        for (int ic = 0; ic < locals.size(); ++ic) {
            m_list[ic].setup(locals, aliens, ic);
        }
    }

    /// @brief Собрать актуальные данные о соседях и сиблингах для части ячеек.
    /// Функция инициализирует переменные max_neib_wanted_level и max_sibs_flag.
    void collect_partial(int from, int to) {
        for (int i = from; i < to; ++i) {
            m_list[i].collect();
        }
    }
};

/// @brief Выполняет обход части ячеек хранилища, повышает флаги адаптации для
/// ячеек, соседи которых убегают вперед (хотят уровень на 2 выше), а также для
/// ячеек, которые хотят огрубиться в отличае от своих сиблингов.
/// @param locals Ссылка на хранилище
/// @param around Ссылка на массив с окружением ячеек
/// @param from, to Диапазон индексов ячеек для обхода
/// @return true если хотя бы одна ячейка в диапазоне изменила свой флаг
template <int dim>
bool round_partial(Storage& locals, const CellsAroundList<dim>& around, int from, int to) {
    bool changed = false;

    for (int ic = from; ic < to; ++ic) {
        scrutiny_check(ic < locals.size(), "round_partial: ic >= locals.size()")

        auto cell = locals[ic];
        if (cell.flag() > 0) {
            continue;
        }

        // Максимальный желаемый уровень соседей
        auto neibs_wanted_lvl = around[ic].max_neib_wanted_level;

        // Ячейка не хочет меняться
        if (cell.flag() == 0) {
            // Один из соседей хочет слишком высокий уровень
            if (neibs_wanted_lvl > cell.level() + 1) {
                cell.geom().flag = 1;
                changed = true;
            }
        }
        else {
            // Ячейка хочет огрубиться, cell.flag() < 0
            scrutiny_check(cell.flag() < 0, "round_partial: wrong assumption")

            // Один из соседей хочет слишком высокий уровень
            if (neibs_wanted_lvl > cell.level() + 1) {
                cell.geom().flag = 1;
                changed = true;
                continue;
            }

            // Один из соседей хочет достаточно высокий уровень
            if (neibs_wanted_lvl > cell.level()) {
                cell.geom().flag = 0;
                changed = true;
                continue;
            }

            // Один из сиблингов не хочет огрубляться
            if (around[ic].max_sibs_flag > -1) {
                cell.geom().flag = 0;
                changed = true;
            }
        }
    }

    return changed;
}

/// @brief Выполняет функцию round_partial для всех ячеек в однопоточном режиме
template <int dim>
bool round(Storage& locals, const CellsAroundList<dim>& around) {
    return round_partial<dim>(locals, around, 0, locals.size());
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Выполняет функцию round_partial для всех ячеек в многопоточном режиме
/// @param threads Ссылка на пул тредов
template <int dim>
bool round(Storage& locals, const CellsAroundList<dim>& around, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        return round_partial<dim>(locals, around, 0, locals.size());
    }
    std::vector<std::future<bool>> results(num_tasks);

    std::int bin = locals.size() / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(round_partial<dim>,
                              std::ref(locals), std::ref(around),
                              pos, std::min(pos + bin, locals.size())
        );
        pos += bin;
    }
    bool changed = false;
    for (auto &result: results)
        changed = changed || result.get();
    
    return changed;
}
#endif

/// @brief Простая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище с сохранением баланса
/// 1:2 у соседних ячеек. Баланс достигается повышением флагов адаптации
/// (-1 и 0) у части ячеек.
/// @param cells Ссылка на хранилище ячеек
/// @param max_level Максимальный уровень ячеек
/// @details Детали алгоритма. На первом этапе собирается информация об
/// окружении ячеек (ссылки на флаги и уровни соседей и сиблингов), которая
/// заносится в отдельный массив. Сбор информации необходим для быстрого
/// доступа к меняющимся в течении алгоритма данным соседей.
/// На втором этапе выполняется многократный обход всех ячеек с повышением
/// их флагов адатации (если необходимо). Итерации заканчиваются, если во
/// время обхода ни одна из ячеек не изменила свой флаг. Количество итераций
/// (обходов), требуемое для достижения баланса примерно равно максимальному
/// количеству уровней. Алгоритм достаточно прост в реализации, для него легко
/// (и достаточно эффективно) реализется многопоточность, также алгоритм легко
/// обобщается на многопроцессорную систему.
template<int dim>
void balance_flags(Storage &locals, Storage& aliens, int max_level) {
    static Stopwatch restrictions_timer;
    static Stopwatch setup_around_timer;
    static Stopwatch round_timer;

    restrictions_timer.resume();
    base_restrictions<dim>(locals, max_level);
    restrictions_timer.stop();

    setup_around_timer.resume();
    
    CellsAroundList<dim> around(locals, aliens);
    setup_around_timer.stop();

    round_timer.resume();
    bool changed = true;
    while (changed) {
        around.collect();
        changed = round(locals, around);
    }
    round_timer.stop();

#if CHECK_PERFORMANCE
    std::cout << "    Restriction elapsed: " << restrictions_timer.seconds() << " sec\n";
    std::cout << "    Setup around elapsed: " << setup_around_timer.seconds() << " sec\n";
    std::cout << "    Round elapsed: " << round_timer.seconds() << " sec\n";
#endif
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
void balance_flags_slow(Storage &locals, Storage& aliens, int max_level) {
    if (locals.empty())
        return;

    auto dim = locals[0].dim();

    if (dim < 3) {
        amr::balance_flags<2>(locals, aliens, max_level);
    }
    else {
        amr::balance_flags<3>(locals, aliens, max_level);
    }

#if SCRUTINY
    amr::check_flags(locals, locals, max_level);
#endif
}

#ifdef ZEPHYR_ENABLE_MPI

/// @brief Простая версия функции балансировки флагов.
/// @details Смотреть однопоточную версию.
template<int dim>
void balance_flags(
        Decomposition &decomposition,
        int max_level
        if_multithreading(, ThreadPool& threads = dummy_pool))
{
    using zephyr::performance::timer::Stopwatch;
    static Stopwatch restrictions_timer;
    static Stopwatch setup_around_timer;
    static Stopwatch round_timer;

    Network& network = decomposition.network();
    Storage& locals = decomposition.inner_elements();
    Storage& aliens = decomposition.outer_elements();

    restrictions_timer.resume();
    base_restrictions<dim>(locals, max_level if_multithreading(, threads));
    restrictions_timer.stop();

    setup_around_timer.resume();
    CellsAroundList<dim> around(locals, aliens if_multithreading(, threads));
    setup_around_timer.stop();

    round_timer.resume();
    bool changed = true;
    while (changed) {
        decomposition.exchange_start();
        decomposition.exchange_end();

        around.collect(if_multithreading(threads));
        changed = round(locals, around if_multithreading(, threads));
        changed = network.max(int(changed));
    }
    round_timer.stop();

#if CHECK_PERFORMANCE
    std::cout << "    Restriction elapsed: " << restriction_timer.wall() << "\n";
    std::cout << "    Setup around elapsed: " << setup_around_timer.times().wall() << "\n";
    std::cout << "    Round elapsed: " << round_timer.times().wall() << "\n";
#endif
}

/// @brief Специализация для процессов без ячеек
template<>
void balance_flags<0>(
        Decomposition &decomposition,
        int max_level
        if_multithreading(, ThreadPool& threads))
{
    Network& network = decomposition.network();

    bool changed = true;
    while (changed) {
        decomposition.exchange_start();
        decomposition.exchange_end();
        changed = network.max(0);
    }
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
void balance_flags_slow(
        Decomposition &decomposition,
        int max_level
        if_multithreading(, ThreadPool& threads))
{
    Storage& cells = decomposition.inner_elements();

    if (cells.empty()) {
        amr::balance_flags<0>(decomposition, max_level if_multithreading(, threads));
    }
    else {
        auto dim = cells[0][element].dimension;

        if (dim < 3) {
            amr::balance_flags<2>(decomposition, max_level if_multithreading(, threads));
        }
        else {
            amr::balance_flags<3>(decomposition, max_level if_multithreading(, threads));
        }
    }

#if SCRUTINY
    auto& aliens = decomposition.outer_elements();
    amr::check_flags(cells, max_level, aliens);
#endif
}

#endif

} // namespace amr
} // namespace mesh
} // namespace zephyr