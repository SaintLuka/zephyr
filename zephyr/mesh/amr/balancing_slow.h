/// @file Файл содержит реализацию простой функции балансировки флагов адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для
/// разработчиков. Хотя вряд ли кто-то когда-нибудь попробует в этом разобраться.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/balancing_restrictions.h>
#include <zephyr/io/vtu_file.h>

namespace zephyr::mesh::amr {

// Окрестность ячейки. Структура содержит ссылки на данные (уровни и флаги)
// соседних ячеек и сиблингов, используется для быстрого доступа к информации
// об обновлении соседей или сиблингов.
// Для ячеек с флагом = 1 данные не хранятся (нет необходимости), для
// ячеек с флагом = 0 не хранятся данные о сиблингах.
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
    /// @param ic Целевая ячейка, для которой определяется окрестность
    /// @param locals Ссылка на локальное хранилище
    /// @param aliens Ссылка на хранилище ячеек с других процессов
    void setup(index_t ic, AmrCells &locals, AmrCells& aliens) {
        scrutiny_check(ic < locals.size(), "setup: ic >= locals.size()")

        neib_count = 0;
        sibs_count = 0;

        // Эти ячейки не интересуют
        if (locals.flag[ic] > 0) {
            return;
        }

        auto &adj = locals.faces.adjacent;

        // Поиск соседей
        for (auto iface: locals.faces_range(ic)) {
            if (locals.faces.is_undefined(iface) ||
                locals.faces.is_boundary(iface)) {
                continue;
            }

#if SCRUTINY
            if (adj.rank[iface] == locals.rank[ic]) {
                if (adj.index[iface] >= locals.size()) {
                    std::cerr << "rank: " << adj.rank[iface] << "\n";
                    std::cerr << "locals.size: " << locals.size() << "\n";
                    std::cerr << "adj.index: " << adj.index[iface] << "\n";
                    throw std::runtime_error("Vicinity::setup error #11");
                }
            } else {
                if (adj.alien[iface] >= aliens.size()) {
                    throw std::runtime_error("Vicinity::setup error #2");
                }
            }
#endif
            auto [neibs, jc] = adj.get_neib(iface, locals, aliens);

            neib_levels[neib_count] = neibs.level[jc];
            neib_flags[neib_count] = &neibs.flag[jc];
            ++neib_count;
        }
        scrutiny_check(neib_count <= max_neibs_size, "nei_count > max_neibs_size")

        // Поиск сиблингов для ячеек с флагом -1
        if (locals.flag[ic] < 0) {
#if SCRUTINY
            if (!can_coarse<dim>(locals, ic)) {
                locals.print_info(ic);
            }
            scrutiny_check(can_coarse<dim>(locals, ic), "Can't coarse")
#endif

            // Код выполняется только если can_coarse было истинно,
            // поэтому можно не волноваться о наличии сиблингов
            auto sibs = get_siblings<dim>( locals, ic);
            for (auto is: sibs) {
                scrutiny_check(is < locals.size(), "Vicinity: Sibling index out of range")

                sibs_flags[sibs_count] = &locals.flag[is];
                ++sibs_count;
            }

            scrutiny_check(sibs_count == max_sibs_size, "Wrong siblings count")
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
            max_sibs_flag = std::max(max_sibs_flag, static_cast<int>(*sibs_flags[i]));
        }
    }
};

// Массив структур типа Vicinity.
// Для ячеек с флагом = 1 данные не хранятся (нет необходимости),
// для ячеек с флагом = 0 не хранятся данные о сиблингах.
template <int dim>
struct VicinityList {
    std::vector<Vicinity<dim>> m_list;

    /// @brief Конструктор только с размером
    /// Удобно для создания статического списка
    VicinityList() = default;

    /// @brief Конструктор построения окружения
    /// @param locals Ссылка на локальное хранилище
    /// @param aliens Ссылка на хранилище ячеек с других процессов
    void fill(AmrCells& locals, AmrCells& aliens) {
        m_list.resize(locals.size());
        threads::parallel_for(index_t{0}, index_t{locals.size()},
                [this, &locals, &aliens](index_t ic) {
                    m_list[ic].setup(ic, locals, aliens);
                });
    }

    /// @brief Собрать актуальные данные о соседях и сиблингах для всех ячеек.
    /// Функция инициализирует max_neib_wanted_level и max_sibs_flag.
    void update() {
        threads::for_each( m_list.begin(), m_list.end(),
            [](Vicinity<dim> &ca) { ca.update(); });
    }

    /// @brief Оператор доступа к элементу списка
    const Vicinity<dim>& operator[](int ic) const {
        scrutiny_check(ic < m_list.size(), "ic >= Vicinity.m_list.size()")
        return m_list[ic];
    }
};

/// @brief Повышает флаг адаптации для ячейки, соседи которой убегают вперед
/// (хотят уровень на 2 выше), а также для ячейки, которая хочет огрубиться
/// в отличие от своих сиблингов.
/// @param cell Целевая ячейка
/// @param locals Ссылка на локальное хранилище
/// @param vicinity_list Ссылка на массив с окружением ячеек
/// @return true если ячейка изменила свой флаг
template <int dim>
bool update_flag(index_t ic, AmrCells& locals, const VicinityList<dim>& vicinity_list) {
    scrutiny_check(ic < locals.size(), "update_flag error: ic >= locals.size()")

    if (locals.flag[ic] > 0) { return false; }

    const auto& vicinity = vicinity_list[ic];

    // Максимальный желаемый уровень соседей
    auto neibs_wanted_lvl = vicinity.max_neib_wanted_level;

    // Ячейка не хочет меняться
    if (locals.flag[ic] == 0) {
        // Один из соседей хочет слишком высокий уровень
        if (neibs_wanted_lvl > locals.level[ic] + 1) {
            locals.flag[ic] = 1;
            return true;
        }
    } else {
        // Ячейка хочет огрубиться, cell.flag < 0
        scrutiny_check(locals.flag[ic] < 0, "update_flag: wrong assumption")

        // Один из соседей хочет слишком высокий уровень
        if (neibs_wanted_lvl > locals.level[ic] + 1) {
            locals.flag[ic] = 1;
            return true;
        }

        // Один из соседей хочет достаточно высокий уровень
        if (neibs_wanted_lvl > locals.level[ic]) {
            locals.flag[ic] = 0;
            return true;
        }

        // Один из сиблингов не хочет огрубляться
        if (vicinity.max_sibs_flag > -1) {
            locals.flag[ic] = 0;
            return true;
        }
    }

    return false;
}

/// @brief Выполняет функцию update_flag для всех ячеек
/// @return true если хотя бы одна ячейка изменила свой флаг
template <int dim>
bool flag_balancing_step(AmrCells& locals, const VicinityList<dim>& vicinity_list) {
    // Функция max в данном контексте заменяет логическое "И"
    utils::range<index_t> range(0, locals.size());
    return threads::max(
            range.begin(), range.end(),
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
void balance_flags_slow(AmrCells& locals, int max_level) {
    static Stopwatch base_restrictions_timer;
    static Stopwatch setup_vicinity_timer;
    static Stopwatch flag_balancing_timer;

    static AmrCells aliens;

    base_restrictions_timer.resume();
    base_restrictions<dim>(locals, max_level);
    base_restrictions_timer.stop();

    // Делаем статическим, чтобы не выделять каждый раз память (гениально)
    static VicinityList<dim> vicinity_list;

    setup_vicinity_timer.resume();
    vicinity_list.fill(locals, aliens);
    setup_vicinity_timer.stop();

    flag_balancing_timer.resume();
    bool changed = true;
    while (changed) {
        vicinity_list.update();
        changed = flag_balancing_step(locals, vicinity_list);
    }
    flag_balancing_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "    Restrictions elapsed:   " << std::setw(10) << base_restrictions_timer.milliseconds() << " ms\n";
        std::cout << "    Setup vicinity elapsed: " << std::setw(10) << setup_vicinity_timer.milliseconds() << " ms\n";
        std::cout << "    Flag balancing elapsed: " << std::setw(10) << flag_balancing_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

#ifdef ZEPHYR_MPI

/// @brief Простая версия функции балансировки флагов.
/// @details Смотреть однопоточную версию.
template<int dim>
void balance_flags_slow(AmrCells &locals, AmrCells &aliens,
        int max_level, EuMesh& mesh) {
    using zephyr::utils::Stopwatch;
    static Stopwatch restrictions_timer;
    static Stopwatch setup_vicinity_timer;
    static Stopwatch flag_balancing_timer;

    // Делаем статическим, чтобы не выделять каждый раз память (гениально)
    static VicinityList<dim> vicinity_list;

    restrictions_timer.resume();
    base_restrictions<dim>(locals, max_level);
    restrictions_timer.stop();

    setup_vicinity_timer.resume();
    vicinity_list.fill(locals, aliens);
    setup_vicinity_timer.stop();

    flag_balancing_timer.resume();
    int changed = 1;
    while (changed) {
        // mesh.sync(); ???????
        vicinity_list.update();

        changed = flag_balancing_step(locals, vicinity_list);
        changed = mpi::max(changed);
    }
    flag_balancing_timer.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "    Restrictions elapsed:   " << std::setw(10) << base_restrictions_timer.milliseconds() << " ms\n";
        std::cout << "    Setup vicinity elapsed: " << std::setw(10) << setup_vicinity_timer.milliseconds() << " ms\n";
        std::cout << "    Flag balancing elapsed: " << std::setw(10) << flag_balancing_timer.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}

/// @brief Специализация для процессов без ячеек
template<>
void balance_flags_slow<0>(AmrCells &locals, AmrCells &aliens,
                           int max_level, EuMesh& mesh) {

    // ???
    throw std::runtime_error("IMPLEMENT");
    /*
    bool changed = true;
    while (changed) {
        decomposition.send();
        decomposition.recv();
        changed = network.max(0);
    }
     */
}

#endif

} // namespace zephyr::mesh::amr