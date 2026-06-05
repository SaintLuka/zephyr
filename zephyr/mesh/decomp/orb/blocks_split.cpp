#include <iostream>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <span>
#include <list>
#include <ranges>

#include <zephyr/mesh/decomp/orb/blocks.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh::decomp {

#ifdef ZEPHYR_MPI
using zephyr::utils::mpi;

// Сумма чисел в массиве
inline int sum(const std::vector<int>& arr) {
    return std::accumulate(arr.begin(), arr.end(), 0);
}

// Вывести массив в поток
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& arr) {
    for (const auto& v: arr) os << v << ", ";
    return os;
}

// Минимальное и максимальное значения в контейнере
template<typename Container>
auto minmax_value(const Container& container)
    -> std::pair<typename Container::value_type, typename Container::value_type>
{
    if (container.empty()) {
        throw std::invalid_argument("Container is empty");
    }

    auto [min_it, max_it] = std::minmax_element(std::begin(container), std::end(container));
    return {*min_it, *max_it};
}

// Упорядоченный массив содержит значение?
inline bool contains(const std::vector<int>& arr, int val) {
    auto it = std::ranges::lower_bound(arr, val);
    return it != arr.end() && *it == val;
}

// Объединить отсортированные массивы в один отсортированный массив
inline void merge_sorted(const std::vector<std::vector<double>>& values, std::vector<double>& result) {
    result.clear();

    if (values.empty()) return;

#ifdef ZEPHYR_ASSERTS
    for (auto& arr: values) {
        if (!std::ranges::is_sorted(arr.begin(), arr.end())) {
            throw std::runtime_error("wrong assumption merge_sorted #1");
        }
    }
#endif

    // Тривиальный случай, единственный массив
    if (values.size() == 1) {
        result = values.front();
        return;
    }

    // Подсчитываем общий размер
    size_t total_size = 0;
    for (const auto& arr: values) {
        total_size += arr.size();
    }

    result.reserve(total_size);

    // Индексы для каждого вектора
    std::vector<size_t> indices(values.size(), 0);

    // Слияние всех отсортированных векторов
    while (true) {
        double min_val = std::numeric_limits<double>::max();
        int min_idx = -1;

        for (size_t i = 0; i < values.size(); ++i) {
            if (indices[i] < values[i].size() && values[i][indices[i]] < min_val) {
                min_val = values[i][indices[i]];
                min_idx = i;
            }
        }

        if (min_idx == -1) break;

        result.push_back(min_val);
        indices[min_idx]++;
    }
}

// Распределить значения из контейнера по корзинам
template <class Container>
void bucket_sort(const Container& values, std::span<const double> barriers,
    std::span<std::list<double>> buckets) {
    if (values.empty()) { return; }

    if (barriers.size() < 2 || barriers.size() != buckets.size() + 1) {
        throw std::runtime_error("wrong assumption bucket_sort #1");
    }

#ifdef ZEPHYR_ASSERTS
    if (!std::ranges::is_sorted(barriers)) {
        throw std::runtime_error("wrong assumption bucket_sort #2");
    }
    auto [min_coord, max_coord] = minmax_value(values);
    if (min_coord < barriers.front() || min_coord > barriers.back()) {
        std::cout << "min_coord: " << min_coord << "\n";
        std::cout << "max_coord: " << max_coord << "\n";
        std::cout << "barriers: [" << barriers.front() << ", " << barriers.back() << "]\n";
        throw std::runtime_error("wrong assumption bucket_sort #3");
    }
#endif

    for (auto p: values) {
        int i = std::ranges::lower_bound(barriers, p) - barriers.begin();
        if (i < 1 || i >= barriers.size()) {
            std::cout << barriers.front() <<" -- " << p << " -- " << barriers.back() << "\n";
            throw std::runtime_error("bucket_sort error: bad range");
        }
        buckets[i - 1].push_back(p);
    }
}


// На каждом процессе MPI-группы имеется массив values, собрать все данные
// на корневом процессе группы и отсортировать, на остальных процессах
// очистить массив values.
void sort_gather(const mpi::group& group, std::vector<double>& values) {
    // Сортировка на каждом процессе
    std::ranges::sort(values);

    if (group.single()) return;

    // Собрать все массивы на корневом (rank = 0)
    auto result = group.gather(values);

    // На остальных процессах пустые массивы
    if (result.empty()) {
        values.clear();
        return;
    }

    // Склеить отсортированные массивы
    merge_sorted(result, values);
}

// Распределенная гистограмма
class Histogram {
public:
    Histogram(const mpi::group&    group,    // Группа MPI-процессов
        const std::vector<double>& values,   // Распределенный массив чисел
        const std::vector<int>&    k_stat,   // Индексы, по которым разбить
        const std::vector<double>& barriers, // Предыдущие разделители/барьеры
        const int                  n_parts   // На сколько частей разбить существующие блоки
        ) : group(group) {

        z_assert(barriers.size() > 1, "Bad barriers size");
        const int n_blocks = static_cast<int>(barriers.size()) - 1;

        // Сечения на подблоки
        sections.reserve(n_parts * n_blocks + 1);

        // Заполним подразбиение линейно
        // TODO: проверить варианты со сгущением
        sections.push_back(barriers[0]);
        for (int i = 0; i < n_blocks; ++i) {
            for (int k = 0; k < n_parts; ++k) {
                double x = barriers[i] + (barriers[i + 1] - barriers[i]) * (k + 1) / n_parts;
                sections.push_back(x);
            }
        }

        // Раскидаем точки по подблокам
        bins.resize(n_parts * n_blocks);
        bucket_sort(values, sections, bins);

        // Заполним массив с размерами
        counts.reserve(bins.size());
        for (const auto& b: bins) {
            counts.push_back(static_cast<int>(b.size()));
        }

        // Сумма совпадает с числом значений
        z_assert(values.size() == sum(counts), "Bad sum #195");

        calc_cumulative_counts_(); // Кумулятивная гистограмма
        find_target_bins_(k_stat); // Индексы интересующих столбцов
    }

    // Число элементов в итересующих столбцах гистограммы
    int remaining_values() const {
        int count = 0;
        for (int b: target_bins) {
            count += cumulative[b + 1] - cumulative[b];
        }
        return count;
    }

    // Приближенные разделители, возвращает массив на корневом процессе группы
    std::vector<double> approx_barriers(const std::vector<int>& k_stat) const {
        std::vector<double> barriers(target_bins.size());
        for (int i = 0; i < target_bins.size(); ++i) {
            int idx = target_bins[i];
            // Линейная интерполяция
            double xi = double(k_stat[i] - cumulative[idx]) / (cumulative[idx + 1] - cumulative[idx]);
            barriers[i] = sections[idx] + (sections[idx + 1] - sections[idx]) * xi;
        }
        return barriers;
    }

    // Точные разделители, возвращает массив на корневом процессе группы
    std::vector<double> exact_barriers(const std::vector<int>& k_stat) const {
        std::vector<double> barriers(target_bins.size());
        for (int i = 0; i < target_bins.size(); ++i) {
            int idx = target_bins[i];
            std::vector coords(bins[idx].begin(), bins[idx].end());

            // Объединение и сортировка на корневом процессе группы
            sort_gather(group, coords);
            z_assert(k_stat[i] >= cumulative[idx] && k_stat[i] <= cumulative[idx + 1], "Bad k-stat");
            if (group.master()) {
                int loc_idx = k_stat[i] - cumulative[idx] - 1;
                z_assert(0 <= loc_idx && loc_idx < coords.size(), "Bad loc index exact ORB");
                barriers[i] = (1.0 + 1.1652e-14) * coords[loc_idx];
            }
        }
        return barriers;
    }

    // Выполнить разбиение интересующих столбцов на n_parts частей.
    // Неинтересные столбцы удаляются.
    void specify(const std::vector<int>& k_stat, int n_parts) {
        // Запомним существующие значения
        std::vector<double> prev_sections = std::move(sections);
        std::vector<int> prev_counts = std::move(counts);
        std::vector<std::list<double>> prev_buckets = std::move(bins);

        // Делаем склейку путем заполнения по новой
        const int int_blocks = target_bins.size();
        sections.clear();
        sections.reserve(int_blocks * n_parts + int_blocks + 2);

        counts.clear();
        counts.reserve(int_blocks * n_parts + int_blocks + 1);

        bins.clear();
        bins.reserve(int_blocks * n_parts + int_blocks + 1);

        sections.push_back(prev_sections[0]);

        bool prev_interesting = true;
        for (int i = 0; i < prev_buckets.size(); ++i) {
            // Номер нового блока для добавления
            int curr = static_cast<int>(bins.size());

            if (contains(target_bins, i)) {
                // Интересный блок.
                // Добавим линейно распределенные узлы и пустые корзины
                for (int k = 0; k < n_parts; ++k) {
                    double xi = prev_sections[i] + (prev_sections[i + 1] - prev_sections[i]) * (k + 1) / n_parts;
                    sections.push_back(xi);
                    bins.emplace_back();
                }

                // распределим точки по корзинам
                bucket_sort(prev_buckets[i], std::span(&sections[curr], n_parts + 1), std::span(&bins[curr], n_parts));
                for (int k = 0; k < n_parts; ++k) {
                    counts.push_back(static_cast<int>(bins[curr + k].size()));
                }
                prev_interesting = true;

#ifdef ZEPHYR_ASSERTS
                // Новое разбиение в сумме дает старое число точек
                int prev_count = prev_counts[i];
                int new_count = std::accumulate(counts.data() + curr, counts.data() + curr + n_parts, 0);
                if (prev_count != new_count) {
                    throw std::runtime_error("Blocks::split_exact: bad counts");
                }
#endif
            }
            else {
                // неинтересный блок
                if (prev_interesting) {
                    // Предыдущий интересный, копируем данные
                    sections.push_back(prev_sections[i + 1]);
                    // Не переносим список точек, это же неинтересный блок
                    bins.emplace_back();
                    // Но число точек запоминаем
                    counts.push_back(prev_counts[i]);
                }
                else {
                    // Предыдущий неинтересный, склейка с ним
                    if (bins.empty()) {
                        throw std::runtime_error("Impossible exact ORB error");
                    }
                    sections.back() = prev_sections[i + 1];
                    counts.back() += prev_counts[i];
                }
                prev_interesting = false;
            }
        }

        calc_cumulative_counts_(); // Кумулятивная гистограмма
        find_target_bins_(k_stat); // Индексы интересующих столбцов
    }

private:
    // Кумулятивная гистограмма по всем процессам
    void calc_cumulative_counts_() {
        cumulative.resize(counts.size() + 1);
        cumulative[0] = 0;
        for (int i = 0; i < counts.size(); ++i) {
            cumulative[i + 1] = cumulative[i] + counts[i];
        }
        // Кумулятивная сумма по столбцам гистограммы со всех процессов
        cumulative = group.sum(cumulative);
#ifdef ZEPHYR_ASSERTS
        const int n_values = group.sum(sum(counts));
        if (cumulative.back() != n_values) {
            throw std::runtime_error("Blocks::split_exact: bad comm_summ");
        }
#endif
    }

    // Индексы столбцов, куда попадает k-статистика
    // k_stat -- Индексы k-ой статистики
    void find_target_bins_(const std::vector<int>& k_stat) {
        target_bins.resize(k_stat.size());
        for (int i = 0; i < k_stat.size(); ++i) {
            int m = std::ranges::lower_bound(cumulative, k_stat[i]) - cumulative.begin();
            if (m < 1 || m >= cumulative.size()) {
                throw std::runtime_error("Blocks::split_exact: bad range");
            }
            target_bins[i] = m - 1;
        }
    }


    // Группа MPI-процессов
    mpi::group group;

    // Границы столбцов гистограммы, синхронизованы на всех процессах группы,
    // полностью охватывают интервал.
    std::vector<double> sections;

    // Число значений в каждом столбце гистограммы. Считаются только локальные точки!
    std::vector<int> counts;

    // Кумулятивная гистограмма, собирается по всем процессам группы
    std::vector<int> cumulative;

    // Списки значений в каждом столбце гистограммы. Списки хранятся только для
    // важных частей гистограммы, в этом случае bins[i].size() == counts[i].
    std::vector<std::list<double>> bins;

    // Искомые столбцы гистограммы
    std::vector<int> target_bins;
};


// Поиск k-статистик "в лоб" через сортировку и сборку на одном процессе
void gather_split(
    const mpi::group&          group,    // Группа MPI-процессов
    const std::vector<double>& values,   // Распределенный массив чисел
    const std::vector<int>&    k_stat,   // Индексы, по которым разбить
    std::vector<double>&       barriers  // Предыдущие разделители/барьеры
    ) {
    std::vector<double> sorted = values;
    std::ranges::sort(sorted);
    auto res = group.gather(sorted);
    merge_sorted(res, sorted);
    if (group.master()) {
        for (int i = 0; i < k_stat.size(); ++i) {
            barriers[i + 1] = (1.0 + 1.0e-13) * sorted[k_stat[i]];
        }
    }
}

// Разделить распределенный массив координат на части.
// Обновляет разделители/барьеры на корневом процессе группы.
void split_coords(
    const mpi::group&          group,    // Группа MPI-процессов
    const std::vector<double>& coords,   // Распределенный массив чисел
    const std::vector<int>&    k_stat,   // Индексы, по которым разбить
    std::vector<double>&       barriers  // Предыдущие разделители/барьеры
    ) {

    z_assert(barriers.size() >= 2, "Wrong barriers size");
    z_assert(k_stat.size() + 2 == barriers.size(), "Wrong k_stat size");

    debug_code {
        // Выход за границы барьера
        if (!coords.empty()) {
            auto [min, max] = minmax_value(coords);
            if (min < barriers.front() || max > barriers.back()) {
                std::cout << barriers.front() <<" -- " << min << " -- " << max << " -- " << barriers.back() << "\n";
                throw std::runtime_error("Some point is out of bounds");
            }
        }

        // Число точек на группу процессов
        const int n_points = group.sum(static_cast<int>(coords.size()));

        // k_stat должен возрастать, между соседними может не быть разницы
        for (int i = 0; i < k_stat.size() - 1; ++i) {
            if (k_stat[i + 1] - k_stat[i] <= 0) {
                throw std::runtime_error("k_stat is not sorted");
            }
        }

        // Барьеры будем выбирать от points[k] до points[k + 1], поэтому
        // значения k_stat должны быть в интервале [0, n_points - 1).
        auto [k_min, k_max] = minmax_value(k_stat);
        if (k_min < 0 || k_max + 1 >= n_points) {
            throw std::runtime_error("k_stat is out of bounds");
        }
    }

    // Для небольшого числа значений: сбор на одном процессе
    int n_points = group.sum(static_cast<int>(coords.size()));
    if (n_points < 10'000) {
        gather_split(group, coords, k_stat, barriers);
        return;
    }

    // Я бы делил примерно на столько частей, тогда максимальная
    // пересылка получится порядка 10'000 значений
    const int M = 10'000 / int(group.size() * barriers.size());

    Histogram hist(group, coords, k_stat, barriers, M);
    while (hist.remaining_values() > 10'000) {
        hist.specify(k_stat, M);
    }

    //auto sections = hist.approx_barriers(k_stat);
    auto sections = hist.exact_barriers(k_stat);
    if (group.master()) {
        z_assert(sections.size() + 2 == barriers.size(), "new barriers bad size");
        for (int i = 0; i < sections.size(); ++i) {
            barriers[i + 1] = sections[i];
        }
    }
}

// Раскидать по процессам в соответствии с рангами
inline void scatter_by_rank(const mpi::group& group,
    const std::vector<int>& ranks,
    std::vector<double>& coords_x,
    std::vector<double>& coords_y,
    std::vector<double>& coords_z) {
    z_assert(coords_x.size() == ranks.size(), "scatter_by_rank: bad size #1");
    z_assert(coords_y.size() == ranks.size(), "scatter_by_rank: bad size #2");
    z_assert(coords_z.size() == ranks.size(), "scatter_by_rank: bad size #3");

    int nproc = group.size();

    // Подсчитать send_counts
    std::vector<int> send_counts(nproc, 0);
    for (int r: ranks) {
        if (r < 0 || r >= nproc) {
            throw std::runtime_error("Invalid rank in `rank` array");
        }
        ++send_counts[r];
    }

    // Отправить send_counts → получить recv_counts
    std::vector<int> recv_counts(nproc, 0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, group.comm());

    // Смещения и общий размер входящих данных
    std::vector<int> send_offsets(nproc);
    send_offsets[0] = 0;
    for (int i = 1; i < nproc; ++i) {
        send_offsets[i] = send_offsets[i-1] + send_counts[i-1];
    }

    std::vector<int> recv_offsets(nproc);
    recv_offsets[0] = 0;
    for (int i = 1; i < nproc; ++i) {
        recv_offsets[i] = recv_offsets[i-1] + recv_counts[i-1];
    }
    int recv_count = recv_offsets.back() + recv_counts.back();

    // Сгруппировать по получателям
    std::vector<double> x_send_buff(coords_x.size());
    std::vector<double> y_send_buff(coords_y.size());
    std::vector<double> z_send_buff(coords_z.size());

    std::vector<int> pos = send_offsets;
    for (int i = 0; i < ranks.size(); ++i) {
        int dest = ranks[i];
        int idx = pos[dest];
        x_send_buff[idx] = coords_x[i];
        y_send_buff[idx] = coords_y[i];
        z_send_buff[idx] = coords_z[i];
        ++pos[dest];
    }

    // Подготовим буфер для получения
    std::vector<double> x_recv_buff(recv_count);
    std::vector<double> y_recv_buff(recv_count);
    std::vector<double> z_recv_buff(recv_count);

    // Выполнить обмен
    std::vector<MPI_Request> send_reqs(3 * nproc, MPI_REQUEST_NULL);
    std::vector<MPI_Request> recv_reqs(3 * nproc, MPI_REQUEST_NULL);
    int nsend = 0, nrecv = 0;

    for (int r = 0; r < nproc; ++r) {
        if (send_counts[r] > 0) {
            MPI_Isend(x_send_buff.data() + send_offsets[r], send_counts[r], MPI_DOUBLE, r, 0, group.comm(), &send_reqs[nsend++]);
            MPI_Isend(y_send_buff.data() + send_offsets[r], send_counts[r], MPI_DOUBLE, r, 1, group.comm(), &send_reqs[nsend++]);
            MPI_Isend(z_send_buff.data() + send_offsets[r], send_counts[r], MPI_DOUBLE, r, 2, group.comm(), &send_reqs[nsend++]);
        }
        if (recv_counts[r] > 0) {
            MPI_Irecv(x_recv_buff.data() + recv_offsets[r], recv_counts[r], MPI_DOUBLE, r, 0, group.comm(), &recv_reqs[nrecv++]);
            MPI_Irecv(y_recv_buff.data() + recv_offsets[r], recv_counts[r], MPI_DOUBLE, r, 1, group.comm(), &recv_reqs[nrecv++]);
            MPI_Irecv(z_recv_buff.data() + recv_offsets[r], recv_counts[r], MPI_DOUBLE, r, 2, group.comm(), &recv_reqs[nrecv++]);
        }
    }

    // Дождёмся завершения всех операций
    MPI_Waitall(nsend, send_reqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(nrecv, recv_reqs.data(), MPI_STATUSES_IGNORE);

    // Заменить значения на полученные
    coords_x.swap(x_recv_buff);
    coords_y.swap(y_recv_buff);
    coords_z.swap(z_recv_buff);
}

// Извлечь барьеры в форме массива из простых чисел
inline std::vector<double> get_barriers(const std::vector<barrier>& arr) {
    std::vector<double> res(arr.size());
    for (int i = 0; i < arr.size(); ++i) {
        res[i] = arr[i];
    }
    return res;
}

// Установить барьеры по заданному массиву чисел
inline void set_barriers(std::vector<barrier>& arr, const std::vector<double>& barriers) {
    z_assert(arr.size() == barriers.size(), "set_barriers: bad size");
    for (int i = 0; i < arr.size(); ++i) {
        arr[i] = barriers[i];
    }
}
#endif // ZEPHYR_MPI

void Blocks::split_exact(const std::vector<Vector3d>& points) {
#ifndef ZEPHYR_MPI
    std::cerr << "Blocks::split_exact: implemented only for MPI programs\n";
    return;
#else
    if (mpi::single()) {
        if (m_size != 1) {
            throw std::runtime_error("Blocks::split_exact error: single MPI process, but not single block");
        }
        return;
    }

    // Массивы координат после отображения
    std::vector<double> coords_x(points.size());
    std::vector<double> coords_y(points.size());
    std::vector<double> coords_z(points.size());
    for (int i = 0; i < points.size(); ++i) {
        Vector3d v = m_map(points[i]);
        coords_x[i] = v.x();
        coords_y[i] = v.y();
        coords_z[i] = v.z();
    }

    // Повторно проставить ранги у точек и перераспределить точки между процессами в группе
    auto redistribute = [this, &coords_x, &coords_y, &coords_z](const mpi::group& group) {
        if (!group || group.single()) return;

        z_assert(coords_x.size() == coords_y.size(), "Blocks::redistribute: bad size #1");
        z_assert(coords_x.size() == coords_z.size(), "Blocks::redistribute: bad size #2");

        // Определить новый ранг для каждой точки
        std::vector<int> ranks(coords_x.size());
        for (int i = 0; i < ranks.size(); ++i) {
            int r = this->rank_mapped({coords_x[i], coords_y[i], coords_z[i]});
            ranks[i] = group.local_rank(r);
        }

        // Раскидать точки в соответствии с рангами
        scatter_by_rank(group, ranks, coords_x, coords_y, coords_z);

#ifdef ZEPHYR_ASSERTS
        // Проверить, что все точки на своих процессах
        for (int i = 0; i < coords_x.size(); ++i) {
            int r = this->rank_mapped({coords_x[i], coords_y[i], coords_z[i]});
            if (r != mpi::rank()) {
                throw std::runtime_error("Blocks::redistribute: rank mapping failed");
            }
        }
#endif
    };

    // ---------------------------- Разбиение по X ----------------------------

    mpi::group world = mpi::group::world();
    if (m_nx > 1) {
        // Полное число точек
        const int n_points = world.sum(static_cast<int>(coords_x.size()));

        // k-порядковая статистика (индексы точек для разбиения)
        std::vector<int> k_stat(m_nx - 1);
        int n_proc = 0;
        for (int i = 0; i < m_nx - 1; ++i) {
            n_proc += size(i);
            k_stat[i] = n_proc * n_points / m_size;
        }

        auto barriers = get_barriers(m_x_coords);
        split_coords(world, coords_x, k_stat, barriers);

        mpi::broadcast(world.root_rank(), barriers);
        set_barriers(m_x_coords, barriers);
    }

    // ---------------------------- Разбиение по Y ----------------------------

    // Коммуникаторы по "слоям"
    std::map<int, mpi::group> groups_x;
    for (int i = 0; i < m_nx; ++i) {
        if (m_ny[i] < 2) continue;

        // Собрать ранги блоков в "слоях"
        std::vector<int> rs; rs.reserve(size(i));
        for (int j = 0; j < m_ny[i]; ++j) {
            for (int k = 0; k < m_nz[i][j]; ++k) {
                rs.push_back(m_ranks[i][j][k]);
            }
        }
        z_assert(rs.size() == size(i), "wrong blocks count");

        // Нулевой коммуникатор для процессов, которые не участвуют
        groups_x[i] = mpi::group::subgroup(rs);
    }

    // Если нужно искать барьеры в подгруппах, то раскидаем координаты снова
    if (m_map.dim() > 1) {
        redistribute(world);
    }
    else {
        return; // одномерная декомпозиция, всё готово
    }

    // Выполнить поиск барьеров внутри каждой группы
    for (const auto& [i, group]: groups_x) {
        if (!group) continue;

        int n_points = group.sum(static_cast<int>(coords_y.size()));

        // k-порядковая статистика (индексы точек для разбиения)
        std::vector<int> k_stat(m_ny[i] - 1);
        int n_proc = 0;
        for (int j = 0; j < m_ny[i] - 1; ++j) {
            n_proc += size(i, j);
            k_stat[j] = n_proc * n_points / size(i);
        }

        auto barriers = get_barriers(m_y_coords[i]);
        split_coords(group, coords_y, k_stat, barriers);
        if (group.master()) {
            set_barriers(m_y_coords[i], barriers);
        }
    }

    // Правильные барьеры только у мастера каждой группы
    for (const auto& [i, group]: groups_x) {
        std::vector<double> barriers(m_y_coords[i].size());
        if (group.master()) {
            barriers = get_barriers(m_y_coords[i]);
        }
        mpi::broadcast(group.root_rank(), barriers);
        set_barriers(m_y_coords[i], barriers);
    }

    // ---------------------------- Разбиение по Z ----------------------------

    // Коммуникаторы по "столбцам"
    std::map<std::tuple<int, int>, mpi::group> groups_xy;
    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < m_ny[i]; ++j) {
            if (m_nz[i][j] < 2) continue;

            // Нулевой коммуникатор для процессов, которые не участвуют
            groups_xy[{i, j}] = mpi::group::subgroup(m_ranks[i][j]);
        }
    }

    // Если нужно искать барьеры в подгруппах, то раскидаем координаты снова
    if (m_map.dim() > 2) {
        for (auto& group: groups_x | std::views::values) {
            redistribute(group);
        }
    }
    else {
        return; // двумерная декомпозиция, всё готово
    }

    // Выполнить поиск барьеров внутри каждой группы
    for (const auto& [idx, group]: groups_xy) {
        if (!group) continue;

        int i = std::get<0>(idx);
        int j = std::get<1>(idx);

        int n_points = group.sum(static_cast<int>(coords_z.size()));

        // k-порядковая статистика (индексы точек для разбиения)
        std::vector<int> k_stat(m_nz[i][j] - 1);
        for (int k = 0; k < m_nz[i][j] - 1; ++k) {
            k_stat[k] = (k + 1) * n_points / m_nz[i][j];
        }

        auto barriers = get_barriers(m_z_coords[i][j]);
        split_coords(group, coords_z, k_stat, barriers);
        if (group.master()) {
            set_barriers(m_z_coords[i][j], barriers);
        }
    }

    // Правильные барьеры только у мастера каждой группы
    for (const auto& [idx, group]: groups_xy) {
        int i = std::get<0>(idx);
        int j = std::get<1>(idx);

        std::vector<double> barriers(m_z_coords[i][j].size());
        if (group.master()) {
            barriers = get_barriers(m_z_coords[i][j]);
        }
        mpi::broadcast(group.root_rank(), barriers);
        set_barriers(m_z_coords[i][j], barriers);
    }
#endif
}

} // namespace zephyr::mesh::decomp
