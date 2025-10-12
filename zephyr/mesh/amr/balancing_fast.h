// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

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
    std::vector<index_t> n_coarse; ///< Число ячеек каждого уровня, которые хотят огрубиться
    std::vector<index_t> n_retain; ///< Число ячеек каждого уровня, которые хотят сохраниться

    /// @brief Списки ячеек по уровням, которые хотят огрубиться
    std::vector<std::vector<index_t>> coarse;

    /// @brief Списки ячеек по уровням, которые хотят сохраниться
    std::vector<std::vector<index_t>> retain;

    /// @brief Однопоточный конструктор класса
    /// @details Вызывается конструктор CellsByLevelPartial для всего диапазона
    /// ячеек хранилища, затем данные перемещаются
    CellsByLevel(AmrCells &cells, int max_level) {
        serial_constructor(cells, max_level);
    }

#ifdef ZEPHYR_MPI
    CellsByLevel(Tourism& tourism, AmrCells &locals, AmrCells& aliens, int max_level) {
        // Также пересылает все флаги
        coarse_to_send = std::vector(max_level, std::vector<std::vector<int>>(mpi::size()));
        retain_to_send = std::vector(max_level, std::vector<std::vector<int>>(mpi::size()));

        coarse_to_recv = std::vector(max_level, std::vector<std::vector<int>>(mpi::size()));
        retain_to_recv = std::vector(max_level, std::vector<std::vector<int>>(mpi::size()));

        const auto& router = tourism.m_cell_route;
        const auto& border_indices = tourism.border_indices();

        for (int r = 0; r < mpi::size(); ++r) {
            // Индексы на отправку
            for (index_t j: router.send_indices(r)) {
                index_t ic = border_indices[j];

                // Теперь можно не делать prepare
                tourism.m_border.flag[j] = locals.flag[ic];

                if (locals.level[ic] == max_level || locals.flag[ic] > 0) {
                    continue;
                }
                if (locals.flag[ic] == 0) {
                    retain_to_send[locals.level[ic]][r].push_back(ic);
                }
                else {
                    coarse_to_send[locals.level[ic]][r].push_back(ic);
                }
            }
        }
        // Уже подготовили border
        auto send_flag = tourism.isend<MpiTag::FLAG>();
        auto recv_flag = tourism.irecv<MpiTag::FLAG>(aliens);

        // TODO: Сортировка по тредам
        serial_constructor(locals, max_level);

        // Получили все флаги
        send_flag.wait();
        recv_flag.wait();

        // Теперь прикинем индексы
        for (int r = 0; r < mpi::size(); ++r) {
            // Индексы на получение
            for (index_t ic: router.recv_indices(r)) {
                if (aliens.level[ic] == max_level || aliens.flag[ic] > 0) {
                    continue;
                }
                if (aliens.flag[ic] == 0) {
                    retain_to_recv[aliens.level[ic]][r].push_back(ic);
                }
                else {
                    coarse_to_recv[aliens.level[ic]][r].push_back(ic);
                }
            }
        }

        // Зададим буферы
        send_flags_buf.resize(mpi::size());
        recv_flags_buf.resize(mpi::size());

        send_requests.resize(mpi::size(), MPI_REQUEST_NULL);
        recv_requests.resize(mpi::size(), MPI_REQUEST_NULL);
    }

    void send_retain(AmrCells &locals, AmrCells &aliens, int lvl, int mpi_tag) {
        send_flags_impl(mpi_tag, locals.flag, retain_to_send[lvl], aliens.flag, retain_to_recv[lvl]);
    }

    void send_coarse(AmrCells& locals, AmrCells& aliens, int lvl, int mpi_tag) {
        send_flags_impl(mpi_tag, locals.flag, coarse_to_send[lvl], aliens.flag, coarse_to_recv[lvl]);
    }

    /// @brief Отправить только часть флагов из массива local_flags в массив alien_flags
    /// @param local_flags Ссылка на массив локальных флагов адаптации
    /// @param alien_flags Ссылка на массив флагов адаптации с других процессов
    /// @param send_indices Индексы флагов из массива local_flags, которые необходимо отправить
    /// другим процессам, send_indices[rank] нужно отправить на rank.
    /// @param recv_indices Индексы флагов из массива alien_flags, которые мы получаем с других
    /// процессов, recv_indices[rank] получаем с других процессов
    void send_flags_impl(int mpi_tag,
        const std::vector<int>& local_flags, const std::vector<std::vector<index_t>>& send_indices,
              std::vector<int>& alien_flags, const std::vector<std::vector<index_t>>& recv_indices) {

        for (int r = 0; r < mpi::size(); ++r) {
            send_requests[r] = MPI_REQUEST_NULL;
            if (!send_indices[r].empty()) {
                // Запакуем в send_flags
                send_flags_buf[r].resize(send_indices[r].size());
                for (index_t i = 0; i < send_indices[r].size(); ++i) {
                    send_flags_buf[r][i] = local_flags[send_indices[r][i]];
                }
                // Асинхронная отправка
                MPI_Isend(send_flags_buf[r].data(), send_indices[r].size(),
                          mpi::type<int>(), r, mpi_tag, mpi::comm(), &send_requests[r]);
            }
        }

        for (int r = 0; r < mpi::size(); ++r) {
            recv_requests[r] = MPI_REQUEST_NULL;
            if (!recv_indices[r].empty()) {
                recv_flags_buf[r].resize(recv_indices[r].size());
                // Асинхронное получение
                MPI_Irecv(recv_flags_buf[r].data(), recv_indices[r].size(),
                          mpi::type<int>(), r, mpi_tag, mpi::comm(), &recv_requests[r]);
            }
        }

        // Ждем завершения отправки / получения
        for (int r = 0; r < mpi::size(); ++r) {
            if (!send_indices[r].empty()) {
                MPI_Wait(&send_requests[r], MPI_STATUS_IGNORE);
            }
        }
        for (int r = 0; r < mpi::size(); ++r) {
            if (!recv_indices[r].empty()) {
                MPI_Wait(&recv_requests[r], MPI_STATUS_IGNORE);
            }
        }

        // Переносим флаги в alien массив
        for (int r = 0; r < mpi::size(); ++r) {
            for (index_t i = 0; i < recv_indices[r].size(); ++i) {
                alien_flags[recv_indices[r][i]] = recv_flags_buf[r][i];
            }
        }
    }
#endif

private:
    /// @brief Однопоточный конструктор класса
    void serial_constructor(AmrCells& cells, int max_level) {
        CellsByLevelPartial part(cells, max_level, 0, cells.size());
        n_coarse = std::move(part.n_coarse);
        n_retain = std::move(part.n_retain);
        coarse = std::move(part.coarse);
        retain = std::move(part.retain);
    }

#ifdef ZEPHYR_MPI
    /// @brief Массив индексов local (border) ячеек, флаги которых retain / coarse,
    /// и уровень меньше максимального coarse_to_send[lvl][rank]
    std::vector<std::vector<std::vector<int>>> coarse_to_send;
    std::vector<std::vector<std::vector<int>>> retain_to_send;

    /// @brief Массив индексов alien ячеек, флаги которых retain / coarse,
    /// и уровень меньше максимального coarse_to_recv[lvl][rank]
    std::vector<std::vector<std::vector<int>>> coarse_to_recv;
    std::vector<std::vector<std::vector<int>>> retain_to_recv;

    // Следующие переменные дополнительные, используются при отправках и получениях,
    // Их можно было бы создать локально в функции, но я прелагаю переиспользовать

    // Буферы для отправки / получения флагов
    std::vector<std::vector<int>> send_flags_buf;
    std::vector<std::vector<int>> recv_flags_buf;

    // Реквесты отправки / получения каждому rank'у
    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;
#endif
};

/// @brief Обновляет флаг ячейки под индексом index, которая имеет флаг 0.
/// Повышает флаг адаптации, если один из соседей хочет уровень на два выше
inline void retain_update_flag(index_t ic, AmrCells &locals, AmrCells& aliens) {
    scrutiny_check(ic < locals.size(), "round_1: cell_idx >= cells.size()")
    scrutiny_check(locals.flag[ic] == 0, "retain_update_flag: cell.flag != 0")

    for (auto iface: locals.faces_range(ic)) {
        if (locals.faces.is_undefined(iface) || locals.faces.is_boundary(iface)) {
            continue;
        }

        auto [neibs, jc] = locals.faces.adjacent.get_neib(iface, locals, aliens);
        scrutiny_check(jc < neibs.size(), "retain_update_flag: neib_idx >= cells.n_cells()")

        int neib_wanted = neibs.level[jc] + neibs.flag[jc];
        if (neib_wanted > locals.level[ic] + 1) {
            locals.flag[ic] = 1;
            return;
        }
    }
}

/// @brief Обновляет флаг ячейки, которая имеет флаг -1.
/// Ставит флаг адаптации 0, если сосед хочет уровень на 1 выше,
/// ставит флаг адаптации 1, если сосед хочет уровень на 2 выше,
/// флаги сиблингов не рассматриваются.
inline void coarse_update_flag(index_t ic, AmrCells &locals, AmrCells& aliens) {
    scrutiny_check(ic < locals.size(), "round_2/3: cell_idx >= cells.size()")

    if (locals.flag[ic] >= 0) {
        return;
    }

    for (auto iface: locals.faces_range(ic)) {
        if (locals.faces.is_undefined(iface) || locals.faces.is_boundary(iface)) {
            continue;
        }

        auto [neibs, jc] = locals.faces.adjacent.get_neib(iface, locals, aliens);
        scrutiny_check(jc < neibs.size(), "coarse_update_flag: neib_idx >= cells.n_cells()")

        int neib_wanted = neibs.level[jc] + neibs.flag[jc];

        if (neib_wanted > locals.level[ic]) {
            if (neib_wanted > locals.level[ic] + 1) {
                locals.flag[ic] = 1;
                return;
            } else {
                locals.flag[ic] = 0;
            }
        }
    }
}

/// @brief Обновляет флаг ячейки, которая имеет флаг -1.
/// Ставит флаг адаптации 0, если один из сиблингов не хочет огрубляться.
template <int dim>
void coarse_update_flag_by_sibs(index_t ic, AmrCells &locals) {
    scrutiny_check(ic < locals.size(), "round_4: cell_idx >= cells.size()")

    if (locals.flag[ic] >= 0) {
        return;
    }

    auto sibs = get_siblings<dim>(locals, ic);
    for (auto is: sibs) {
        scrutiny_check(is < locals.size(), "round_4: Sibling index out of range")

        if (locals.flag[is] >= 0) {
            locals.flag[ic] = 0;
            break;
        }
    }
}

/// @brief Обход диапазона ячеек с флагом = 0, данные ячейки могут только повысить
/// свой уровень до 1 за счет высокоуровневых соседей, после выпонения обхода
/// уровень ячеек больше не меняется
/// @param indices Индексы ячеек с флагом = 0 в хранилище
/// @param locals Ссылка на хранилище
inline void round_1(const std::vector<index_t> &indices, AmrCells &locals, AmrCells &aliens) {
    threads::for_each(indices.begin(), indices.end(),
                      retain_update_flag, std::ref(locals), std::ref(aliens));
}

/// @brief Обход диапазона ячеек с флагом = -1, данные ячейки могут повысить
/// свой уровень до 0 или 1 за счет высокоуровневых соседей, флаги сиблингов
/// не используются. После выпонения обхода флаги ячеек ещё могут измениться
/// (@see round_3, @see round_4), но только у тех ячеек, которые не изменили
/// флаги на данном обходе (сохранили флаг -1).
/// @param indices Индексы ячеек с флагом = -1 в хранилище
/// @param cells Ссылка на хранилище
inline void round_2(const std::vector<index_t> &indices, AmrCells &cells, AmrCells& aliens) {
    threads::for_each(indices.begin(), indices.end(),
                      coarse_update_flag, std::ref(cells), std::ref(aliens));
}

/// @brief Обход диапазона ячеек с изначальным флагом = -1, данные ячейки могут
/// повысить свой флаг только до 0 за счет соседей своего же уровня (если при
/// выполнении round_2 у кого-то из соседей появился флаг 1). После выполнения
/// обхода флаги ячеек ещё могут измениться (@see round_4).
/// @param indices Индексы ячеек с исходным флагом = -1 в хранилище
/// @param cells Ссылка на хранилище
inline void round_3(const std::vector<index_t> &indices, AmrCells &cells, AmrCells& aliens) {
    threads::for_each(indices.begin(), indices.end(),
                      coarse_update_flag, std::ref(cells), std::ref(aliens));
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
    threads::for_each(indices.begin(), indices.end(),
                      coarse_update_flag_by_sibs<dim>, std::ref(cells));
}

/// @brief Быстрая версия функции балансировки флагов.
/// Функция меняет флаги адаптации ячеек в хранилище в соответствии с общими
/// ограничениями (ячейка нижнего уровня не огрубляется, ячейка верхнего -
/// не разбивается), а также с сохранением баланса 1:2 у соседних ячеек.
/// Баланс достигается повышением флагов адаптации (-1 и 0) у части ячеек.
/// @param locals Ссылка на хранилище ячеек
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
void balance_flags_fast(AmrCells& locals, int max_level) {
    static Stopwatch restriction_timer;
    static Stopwatch sorting_timer;
    static Stopwatch round_timer_1;
    static Stopwatch round_timer_2;
    static Stopwatch round_timer_3;
    static Stopwatch round_timer_4;
    static int n_step = 0;
    static int n_total_retain = 0;
    static int n_total_coarse = 0;

    static AmrCells aliens;

    restriction_timer.resume();
    base_restrictions<dim>(locals, max_level);
    restriction_timer.stop();

    sorting_timer.resume();
    // TODO: Сортировка по тредам
    CellsByLevel sorted(locals, max_level);
    sorting_timer.stop();

    for (int lvl = max_level - 1; lvl >= 0; --lvl) {
        round_timer_1.resume();
        round_1(sorted.retain[lvl], locals, aliens);
        round_timer_1.stop();

        round_timer_2.resume();
        round_2(sorted.coarse[lvl], locals, aliens);
        round_timer_2.stop();

        round_timer_3.resume();
        round_3(sorted.coarse[lvl], locals, aliens);
        round_timer_3.stop();

        round_timer_4.resume();
        round_4<dim>(sorted.coarse[lvl], locals);
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

#ifdef ZEPHYR_MPI
/// Отличия параллельной реализации?
/// Объяснить безумие с отправкой флагов.
template <int dim>
void balance_flags_fast(Tourism& tourism, AmrCells &locals, AmrCells &aliens, int max_level) {
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
    base_restrictions<dim>(locals, max_level);
    restriction_timer.stop();

    sorting_timer.resume();
    CellsByLevel sorted(tourism, locals, aliens, max_level);
    sorting_timer.stop();

    for (int lvl = max_level - 1; lvl >= 0; --lvl) {
        round_timer_1.resume();
        round_1(sorted.retain[lvl], locals, aliens);
        sorted.send_retain(locals, aliens, lvl, 1000 + lvl);
        round_timer_1.stop();

        round_timer_2.resume();
        round_2(sorted.coarse[lvl], locals, aliens);
        sorted.send_coarse(locals, aliens, lvl, 2000 + lvl);
        round_timer_2.stop();

        round_timer_3.resume();
        round_3(sorted.coarse[lvl], locals, aliens);
        round_timer_3.stop();

        round_timer_4.resume();
        round_4<dim>(sorted.coarse[lvl], locals);
        sorted.send_coarse(locals, aliens, lvl, 3000 + lvl);
        round_timer_4.stop();
    }

#if CHECK_PERFORMANCE
    n_step += 1;
    n_total_coarse += std::accumulate(sorted.n_coarse.begin(), sorted.n_coarse.end(), 0);
    n_total_retain += std::accumulate(sorted.n_retain.begin(), sorted.n_retain.end(), 0);

    static index_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "    Base restrictions: " << std::setw(8) << restriction_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Sort by level:     " << std::setw(8) << sorting_timer.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Round 1 elapsed:   " << std::setw(8) << round_timer_1.milliseconds_mpi() << " ms  ";
        mpi::cout << "    (avg retain: " << std::setw(10) << n_total_retain / n_step << ")\n";
        mpi::cout << "    Round 2 elapsed:   " << std::setw(8) << round_timer_2.milliseconds_mpi() << " ms  ";
        mpi::cout << "    (avg coarse: " << std::setw(10) << n_total_coarse / n_step << ")\n";
        mpi::cout << "    Round 3 elapsed:   " << std::setw(8) << round_timer_3.milliseconds_mpi() << " ms\n";
        mpi::cout << "    Round 4 elapsed:   " << std::setw(8) << round_timer_4.milliseconds_mpi() << " ms\n";
    }
    ++counter;
#endif
}
#endif

} // namespace zephyr::mesh::amr