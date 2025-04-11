#pragma once

#include <vector>
#include <limits>
#include <zephyr/configuration.h>

#ifdef ZEPHYR_TBB
#include <execution>
#include <tbb/global_control.h>
#endif

#ifdef ZEPHYR_OPENMP
#include <omp.h>
#endif

#ifdef ZEPHYR_STD_THREADS
#include <zephyr/utils/thread-pool.h>
#endif

namespace zephyr::utils {

/// @brief Число задач на тред по умолчанию
static const int default_n_tasks_per_thread = 10;

/// @brief Минимальное число элементов на задачу по умолчанию
static const int default_min_elements_per_task = 500;


/// @brief Статический класс. Упрощенный интерфейс для многопоточности.
class threads {
public:

    /// @brief Рекомендуемое число тредов, используется
    /// по умолчанию при вызове on().
    static int recommended();

    /// @brief Включить многопоточность, будут использованы
    /// все доступные ядра, см. recommended()
    static void on();

    /// @brief Включить треды
    /// @param count Число тредов.
    static void on(int count);

    /// @brief Выключить треды
    static void off();

    /// @brief Многопоточность включена?
    static bool active();

    /// @brief Число тредов
    static int count();

    /// @brief Выводит информацию о параллельности
    static void info();

    /// @brief Выполнить функцию для элементов из диапазона
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., целевая
    /// функция может иметь возвращаемое значение, но оно игнорируется.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference>
    static void for_each(Iter begin, Iter end, Func &&func, Args &&... args);

    /// @brief Минимизировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param init Начальное "наибольшое" значение для операции минимизации
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, для данного
    /// типа должен быть определен оператор сравнения < "меньше".
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать значение типа Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type;

    /// @brief Минимизировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, должен быть
    /// арифметическим типом (int, float, double, и т. д.)
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать арифметический тип.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto min(Iter begin, Iter end, Func &&func, Args &&... args)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return min<n_tasks_per_thread, min_elements_per_task>(
                begin, end, std::numeric_limits<Value>::max(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Максимизировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param init Начальное "наименьшее" значение для операции максимизации
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, для данного
    /// типа должен быть определен оператор сравнения > "больше".
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать значение типа Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type;

    /// @brief Максимизировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, должен быть
    /// арифметическим типом (int, float, double, и т. д.)
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать арифметический тип.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto max(Iter begin, Iter end, Func &&func, Args &&... args)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return max<n_tasks_per_thread, min_elements_per_task>(
                begin, end, std::numeric_limits<Value>::lowest(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Суммировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param init Начальное значение для суммирование (нейтральный элемент по сложению)
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, для данного
    /// типа должен быть определен оператор добавления +=.
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать значение типа Value.
    /// @return Массив из частичных сумм (на каждом треде в отдельности), размер
    /// выходного массива будет равен n_threads * n_tasks_per_thread.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto partial_sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, std::vector<Value>>::type;

    /// @brief Суммировать результаты выполенения функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param init Начальное значение для суммирование (нейтральный элемент по сложению)
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, для данного
    /// типа должен быть определен оператор добавления +=.
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать значение типа Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type;

    /// @brief Операция свертки (обобщенное суммирование) результатов выполенения
    /// функции на диапазоне элементов
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param init Начальное для свертки (нейтральный элемент по операции свертки)
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @tparam DeRef Тип итератора после разыменования
    /// @tparam Value Тип возвращаемого значения целевой функции, для данного
    /// типа должен быть определен оператор добавления &=.
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter (DeRef) и набор аргументов Args...,
    /// целевая функция должна возвращать значение типа Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            class Iter, class Func, class... Args,
            class DeRef = typename std::iterator_traits<Iter>::reference,
            class Value = std::invoke_result_t<Func, DeRef, Args...>>
    static auto reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type;

public:
    /// @brief Число тредов
    static int n_threads;

#ifdef ZEPHYR_TBB
    static tbb::global_control m_control;
#endif

#ifdef ZEPHYR_STD_THREADS
    /// @brief Указатель на пул тредов
    static std::unique_ptr<ThreadPool> pool;
#endif
};

template<int n_tasks_per_thread, int min_elements_per_task>
int get_n_tasks(int n_threads, int size) {
    return std::max(1, std::min(n_tasks_per_thread * n_threads, size / min_elements_per_task));
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef>
void threads::for_each(Iter begin, Iter end, Func &&func, Args &&... args) {
#ifdef ZEPHYR_TBB
    std::for_each(std::execution::par,
        begin, end,
        [&func, &args...](DeRef elem) {
            func(elem, std::forward<Args>(args)...);
        });
    return;
#endif

#ifdef ZEPHYR_OPENMP
    #pragma omp parallel for
    for (Iter it = begin; it < end; ++it) {
        func(*it, std::forward<Args>(args)...);
    }
    return;
#endif

#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&func, &args...](const Iter &a, const Iter &b) {
                for (auto it = a; it < b; ++it) {
                    func(*it, std::forward<Args>(args)...);
                }
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return;

    // Выполняем последовательно
    if (n_threads < 2) {
        bin_function(begin, end);
        return;
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<void>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    for (auto &result : results)
        result.wait();
    return;
#endif
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto a, auto b) { return a < b ? a : b; },
        [&func, &args...](auto elem) -> Value {
            return func(elem, std::forward<Args>(args)...);
        });
#endif

#ifdef ZEPHYR_OPENMP
    //#pragma omp declare reduction(minValue:Value: \
    //omp_out = omp_in < omp_out ? omp_in : omp_out)

    Value min_val = init;
    #pragma omp parallel for reduction(min:min_val)
    for (Iter it = begin; it < end; ++it) {
        Value temp = func(*it, std::forward<Args>(args)...);
        if (temp < min_val) min_val = temp;
    }
    return min_val;
#endif

#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res(init);
                Value temp(init);
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp < res) {
                        res = temp;
                    }
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

    // Выполняем последовательно
    if (n_threads < 2) {
        return bin_function(begin, end);
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    Value res(init);
    Value temp(init);
    for (auto &result : results) {
        temp = result.get();
        if (temp < res) {
            res = temp;
        }
    }
    return res;
#endif
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {

#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto a, auto b) { return a > b ? a : b; },
        [&func, &args...](auto elem) -> Value {
            return func(elem, std::forward<Args>(args)...);
        });
#endif

#ifdef ZEPHYR_OPENMP
    #pragma omp declare reduction(maxValue:Value: \
    omp_out = omp_in > omp_out ? omp_in : omp_out)

    Value max_val = init;
    #pragma omp parallel for reduction(maxValue:max_val)
    for (Iter it = begin; it < end; ++it) {
        Value temp = func(*it, std::forward<Args>(args)...);
        if (temp > max_val) max_val = temp;
    }
    return max_val;
#endif

#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res(init);
                Value temp(init);
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp > res) {
                        res = temp;
                    }
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

    // Выполняем последовательно
    if (n_threads < 2) {
        return bin_function(begin, end);
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    Value res(init);
    Value temp(init);
    for (auto &result : results) {
        temp = result.get();
        if (temp > res) {
            res = temp;
        }
    }
    return res;
#endif
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::partial_sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
-> typename std::enable_if<!std::is_void<Value>::value, std::vector<Value>>::type {


#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res(init);
                for (auto it = a; it < b; ++it) {
                    res += func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return {init};

    // Выполняем последовательно
    if (n_threads < 2) {
        return {bin_function(begin, end)};
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    std::vector<Value> res;
    res.reserve(n_tasks);
    for (auto &r: results) {
        res.push_back(r.get());
    }

    return res;
#else
    throw std::runtime_error("Not implemented for other");
#endif
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto a, auto b) { auto c = a; c += b; return c; },
        [&func, &args...](DeRef elem) {
            return func(elem, std::forward<Args>(args)...);
        });
#endif

#ifdef ZEPHYR_OPENMP
    #pragma omp declare reduction(sumValue:Value:omp_out += omp_in)

    Value sum_val = init;
    #pragma omp parallel for reduction(sumValue:sum_val)
    for (Iter it = begin; it < end; ++it) {
        sum_val += func(*it, std::forward<Args>(args)...);
    }
    return sum_val;
#endif

#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res(init);
                for (auto it = a; it < b; ++it) {
                    res += func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

    // Выполняем последовательно
    if (n_threads < 2) {
        return bin_function(begin, end);
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    Value res(init);
    for (auto &result : results) {
        res += result.get();
    }
    return res;
#endif
}

template<int n_tpt, int min_ept,
        class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto a, auto b) { auto c = a; c &= b; return c; },
        [&func, &args...](auto elem) {
            return func(elem, std::forward<Args>(args)...);
        });
#endif

#ifdef ZEPHYR_OPENMP
#pragma omp declare reduction(redValue:Value:omp_out &= omp_in)

    Value reduced = init;
#pragma omp parallel for reduction(redValue:reduced)
    for (Iter it = begin; it < end; ++it) {
        reduced &= func(*it, std::forward<Args>(args)...);
    }
    return reduced;
#endif

#ifdef ZEPHYR_STD_THREADS
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res(init);
                for (auto it = a; it < b; ++it) {
                    res &= func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

    // Выполняем последовательно
    if (n_threads < 2) {
        return bin_function(begin, end);
    }

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<Value>> results;
    results.reserve(n_tasks);

    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    Value res(init);
    for (auto &result : results) {
        res &= result.get();
    }
    return res;
#endif
}

} // namespace zephyr::utils