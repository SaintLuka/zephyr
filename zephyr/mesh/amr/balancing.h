/// @file Файл содержит реализацию простой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для
/// разработчиков. Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_restrictions.h>

namespace zephyr { namespace mesh { namespace amr {

/// @struct Vicinity - окрестность, прикольное слово, да?
/// @brief Простая структура, содержит ссылки на данные (уровни и флаги)
/// сосдених ячеек и сиблингов, требуется для быстрого доступа к информации
/// об обновлении соседей или сиблингов.
/// Для ячеек с флагом = 1 данные не хранятся (нет необходимости), для
/// ячеек с флагом = 0 не хранятся данные о сиблингах.
template<int dim>
struct Vicinity {
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
    Vicinity() : neib_count(0), sibs_count(0) {};

    /// @brief Заполнение списков соседей для ячейки
    /// @param cell Целевая ячейка, для которой определяется окрестность
    /// @param locals Ссылка на локальное хранилище
    /// @param aliens Ссылка на хранилище ячеек с других процессов
    void setup(const AmrCell& cell, AmrStorage &locals, AmrStorage& aliens) {
        scrutiny_check(cell.index < locals.size(), "setup: ic >= locals.size()")

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
                    throw std::runtime_error("Vicinity::setup error #11");
                }
            } else {
                if (adj.ghost >= aliens.size()) {
                    throw std::runtime_error("Vicinity::setup error #2");
                }
            }
#endif
            AmrCell& neib = adj.rank == cell.rank ? locals[adj.index] : aliens[adj.ghost];

            neib_levels[neib_count] = neib.level;
            neib_flags[neib_count] = &neib.flag;
            ++neib_count;
        }
        scrutiny_check(neib_count <= max_neibs_size, "nei_count > max_neibs_size")

        // Поиск сиблингов для ячеек с флагом -1
        if (cell.flag  < 0) {
            sibs_count = 0;
            if (cell.flag < 0) {
                scrutiny_check(can_coarse<dim>(cell, locals), "Can't coarse")

                // Код выполняется только если can_coarse было истинно,
                // поэтому можно не волноваться о начилии сиблингов
                auto sibs = get_siblings<dim>(cell, locals);
                for (auto is: sibs) {
                    scrutiny_check(is < locals.size(), "Vicinity: Sibling index out of range")

                    AmrCell& sib = locals[is];
                    sibs_flags[sibs_count] = &sib.flag;
                    ++sibs_count;
                }

                scrutiny_check(sibs_count == max_sibs_size, "Wrong siblings count")
            }
        }
    }

    /// @brief Собрать актуальные данные о соседях и сиблингах. Функция
    /// инициализирует переменные max_neib_wanted_level и max_sibs_flag.
    void update() {
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

        // Максимальный уровень, который хочет получить сиблинг после адаптации
        // Значение зависит от текущего состояния окружения
        max_sibs_flag = -1;
        for (int i = 0; i < sibs_count; ++i) {
            scrutiny_check(sibs_flags[i] != nullptr, "Undefined neighbor")
            max_sibs_flag = std::max(max_sibs_flag, int(*sibs_flags[i]));
        }
    }
};

/// @struct VicinityList - массив структур типа Vicinity.
/// @brief Для ячеек с флагом = 1 данные не хранятся (нет необходимости),
/// для ячеек с флагом = 0 не хранятся данные о сиблингах.
template <int dim>
struct VicinityList {
    std::vector<Vicinity<dim>> m_list;

    /// @brief Конструктор построения окружения
    /// @param locals Ссылка на локальное хранилище
    /// @param aliens Ссылка на хранилище ячеек с других процессов
    VicinityList(AmrStorage& locals, AmrStorage& aliens)
        : m_list(locals.size()) {
        threads::for_each(
                locals.begin(), locals.end(),
                [this, &locals, &aliens](const AmrStorage::Item &item) {
                    m_list[item.index].setup(item, locals, aliens);
                });
    }

    /// @brief Собрать актуальные данные о соседях и сиблингах для всех ячеек.
    /// Функция инициализирует max_neib_wanted_level и max_sibs_flag.
    void update() {
        threads::for_each(
                m_list.begin(), m_list.end(),
                [](Vicinity<dim> &ca) {
                    ca.update();
                });
    }

    /// @brief Оператор доступа к элементу списка
    const Vicinity<dim>& operator[](int ic) const {
        scrutiny_check(ic < m_list.size(), "ic >= Vicinity.m_list.size()")
        return m_list[ic];
    }
};

/// @brief Повышает флаг адаптации для ячейки, соседи которой убегают вперед
/// (хотят уровень на 2 выше), а также для ячейки, которая хочет огрубиться
/// в отличае от своих сиблингов.
/// @param cell Целевая ячейка
/// @param locals Ссылка на локальное хранилище
/// @param vicinity_list Ссылка на массив с окружением ячеек
/// @return true если ячейка изменила свой флаг
template <int dim>
bool update_flag(AmrCell& cell, AmrStorage& locals, const VicinityList<dim>& vicinity_list) {
    int ic = cell.index;
    auto vicinity = vicinity_list[ic];

    scrutiny_check(ic < locals.size(), "update_flag: ic >= locals.size()")

    if (cell.flag > 0) {
        return false;
    }

    // Максимальный желаемый уровень соседей
    auto neibs_wanted_lvl = vicinity.max_neib_wanted_level;

    // Ячейка не хочет меняться
    if (cell.flag == 0) {
        // Один из соседей хочет слишком высокий уровень
        if (neibs_wanted_lvl > cell.level + 1) {
            cell.flag = 1;
            return true;
        }
    } else {
        // Ячейка хочет огрубиться, cell.flag < 0
        scrutiny_check(cell.flag < 0, "update_flag: wrong assumption")

        // Один из соседей хочет слишком высокий уровень
        if (neibs_wanted_lvl > cell.level + 1) {
            cell.flag = 1;
            return true;
        }

        // Один из соседей хочет достаточно высокий уровень
        if (neibs_wanted_lvl > cell.level) {
            cell.flag = 0;
            return true;
        }

        // Один из сиблингов не хочет огрубляться
        if (vicinity.max_sibs_flag > -1) {
            cell.flag = 0;
            return true;
        }
    }

    return false;
}

/// @brief Выполняет функцию update_flag для всех ячеек
/// @return true если хотя бы одна ячейка изменила свой флаг
template <int dim>
bool flag_balancing_step(AmrStorage& locals, const VicinityList<dim>& vicinity_list) {
    // Функция max в данном контексте заменяет логическое "И"
    return threads::max(
            locals.begin(), locals.end(),
            update_flag<dim>, std::ref(locals), std::ref(vicinity_list)
    );
}

/// @brief Простая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище с сохранением баланса
/// 1:2 у соседних ячеек. Баланс достигается повышением флагов адаптации
/// (-1 и 0) у части ячеек.
/// @param locals Ссылка на локальное хранилище ячеек
/// @param aliens Ссылка на хранилище ячеек с других процессов
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
/// (и достаточно эффективно) реализуется многопоточность, также алгоритм легко
/// обобщается на многопроцессорную систему.
template<int dim>
void balance_flags(AmrStorage &locals, AmrStorage& aliens, int max_level) {
    static Stopwatch base_restrictions_timer;
    static Stopwatch setup_vicinity_timer;
    static Stopwatch flag_balancing_timer;

    base_restrictions_timer.resume();
    base_restrictions<dim>(locals, max_level);
    base_restrictions_timer.stop();

    setup_vicinity_timer.resume();
    VicinityList<dim> vicinity_list(locals, aliens);
    setup_vicinity_timer.stop();

    flag_balancing_timer.resume();
    bool changed = true;
    while (changed) {
        vicinity_list.update();
        changed = flag_balancing_step(locals, vicinity_list);
    }
    flag_balancing_timer.stop();

#if CHECK_PERFORMANCE
    std::cout << "    Restrictions elapsed:   " << base_restrictions_timer.milliseconds() << " sec\n";
    std::cout << "    Setup vicinity elapsed: " << setup_vicinity_timer.milliseconds() << " sec\n";
    std::cout << "    Flag balancing elapsed: " << flag_balancing_timer.milliseconds() << " sec\n";
#endif
}

/// @brief Специализация по умолчанию с автоматическим выбором размерности
/// @param locals Ссылка на локальное хранилище ячеек
/// @param aliens Ссылка на хранилище ячеек с других процессов
/// @param max_level Максимальный уровень ячеек
void balance_flags_slow(AmrStorage &locals, AmrStorage& aliens, int max_level) {
    if (locals.empty())
        return;

    auto dim = locals[0].dim;

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
    AmrStorage& locals = decomposition.inner_elements();
    AmrStorage& aliens = decomposition.outer_elements();

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
    AmrStorage& cells = decomposition.inner_elements();

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