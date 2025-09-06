#pragma once

#include <vector>
#include <limits>
#include <zephyr/configuration.h>

#ifdef ZEPHYR_TBB
#include <execution>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#else
#include <zephyr/utils/thread-pool.h>
#endif

namespace zephyr::utils {

/// @brief Число задач на тред по умолчанию
constexpr int default_n_tasks_per_thread = 10;

/// @brief Минимальное число элементов на задачу по умолчанию
constexpr int default_min_elements_per_task = 500;


/// @brief Упрощенный интерфейс для многопоточности.
class threads {
public:
    /// @brief Инициализация аргументами командной строки, варианты:
    /// --threads=on   -t on   (эквивалентно threads::on(), по умолчанию)
    /// --threads=off  -t off  (эквивалентно threads::off() )
    /// --threads=N    -t N    (эквивалентно threads::on(N) )
    static void init(int argc, char** argv);

    /// @brief Рекомендуемое число потоков, используется
    /// по умолчанию при вызове on().
    static int recommended();

    /// @brief Включить многопоточность, будут использованы
    /// все доступные ядра, см. recommended()
    static void on();

    /// @brief Включить многопоточность
    /// @param count Число потоков.
    static void on(int count);

    /// @brief Выключить многопоточность
    static void off();

    /// @brief Число тредов
    static int count() { return n_threads; }

    /// @brief Многопоточность выключена?
    static bool disabled() { return n_threads < 2; }

    /// @brief Выводит информацию о параллельности
    static void info();

    /// @brief Выполнить функцию для элементов из диапазона [begin, end)
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., целевая
    /// функция может иметь возвращаемое значение, но оно игнорируется.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
        int min_elements_per_task = default_min_elements_per_task,
        typename Iter, class Func, class... Args>
    static void for_each(Iter begin, Iter end, Func &&func, Args &&... args);

    /// @brief Выполнить функцию для чисел в диапазоне [begin, end)
    /// @param func Целевая функция, принимает аргументы (idx, Args...)
    /// @param begin, end Индексы, границы диапазона
    /// @param args Аргументы функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// индекс и набор аргументов Args..., целевая функция может иметь
    /// возвращаемое значение, но оно игнорируется.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Index, class Func, class... Args>
    static void parallel_for(Index begin, Index end, Func &&func, Args &&... args);

    /// @brief Минимизировать результаты выполнения функции на диапазоне
    /// элементов [begin, end).
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param init Начальное "наибольшее" значение
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @tparam Value Тип возвращаемого значения целевой функции; для типа
    /// должен быть определен оператор "меньше" (operator<).
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., возвращаемое
    /// значение должно иметь тип Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, typename Value, class Func, class... Args>
    static Value min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Минимизировать результаты выполнения функции на диапазоне
    /// элементов [begin, end). Функция должна возвращать арифметический тип.
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args...
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, class Func, class... Args>
    static auto min(Iter begin, Iter end, Func &&func, Args &&... args)
    -> std::invoke_result_t<Func, decltype(*begin), Args...> {
        using Value = std::invoke_result_t<Func, decltype(*begin), Args...>;
        static_assert(std::is_arithmetic_v<Value>, "Function must have arithmetic return type");
        return min<n_tasks_per_thread, min_elements_per_task>(
            begin, end, std::numeric_limits<Value>::max(),
            std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Максимизировать результаты выполнения функции на диапазоне
    /// элементов [begin, end).
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param init Начальное "наименьшее" значение
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @tparam Value Тип возвращаемого значения целевой функции; для типа
    /// должен быть определен оператор "больше" (operator>).
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., возвращаемое
    /// значение должно иметь тип Value.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, typename Value, class Func, class... Args>
    static Value max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Максимизировать результаты выполнения функции на диапазоне
    /// элементов [begin, end). Функция должна возвращать арифметический тип.
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args...
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, class Func, class... Args>
    static auto max(Iter begin, Iter end, Func &&func, Args &&... args)
    -> std::invoke_result_t<Func, decltype(*begin), Args...> {
        using Value = std::invoke_result_t<Func, decltype(*begin), Args...>;
        static_assert(std::is_arithmetic_v<Value>, "Function must have arithmetic return type");
        return max<n_tasks_per_thread, min_elements_per_task>(
            begin, end, std::numeric_limits<Value>::lowest(),
            std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Суммировать результаты выполнения функции на диапазоне
    /// элементов [begin, end).
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param init Начальное "нулевое" значение
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @tparam Value Тип возвращаемого значения целевой функции; для типа
    /// должен быть определен оператор "append" (operator+=).
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., для
    /// возвращаемого значения должен быть определен operator+=.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, typename Value, class Func, class... Args>
    static Value sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

    /// @brief Обобщенная операция свёртки для диапазона элементов [begin, end).
    /// @param begin, end Итераторы на начало и конец диапазона
    /// @param init Начальное "нулевое" значение
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргументы целевой функции
    /// @tparam Value Тип возвращаемого значения целевой функции; для типа
    /// должен быть определен оператор operator&=.
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., для
    /// возвращаемого значения должен быть определен operator&=.
    template<int n_tasks_per_thread = default_n_tasks_per_thread,
            int min_elements_per_task = default_min_elements_per_task,
            typename Iter, typename Value, class Func, class... Args>
    static Value reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args);

private:
    /// @brief Число потоков
    static int n_threads;

#ifdef ZEPHYR_TBB
    /// @brief Ограничитель числа потоков для TBB
    static std::unique_ptr<tbb::global_control> m_control;
#else
    /// @brief Указатель на пул тредов
    static std::unique_ptr<ThreadPool> pool;
#endif
};

template<int n_tasks_per_thread, int min_elements_per_task>
int get_n_tasks(int n_threads, int size) {
    return std::max(1, std::min(n_tasks_per_thread * n_threads, size / min_elements_per_task));
}

template<int n_tpt, int min_ept, typename Iter, class Func, class ...Args>
void threads::for_each(Iter begin, Iter end, Func &&func, Args &&... args) {
    if (threads::disabled()) {
        // Последовательное выполнение
        for (auto it = begin; it != end; ++it) {
            func(*it, std::forward<Args>(args)...);
        }
        return;
    }

#ifdef ZEPHYR_TBB
    std::for_each(std::execution::par,
        begin, end,
        [&func, &args...](auto&& elem) {
            func(elem, std::forward<Args>(args)...);
        });
#else
    auto bin_function =
            [&func, &args...](const Iter &a, const Iter &b) {
                for (auto it = a; it < b; ++it) {
                    func(*it, std::forward<Args>(args)...);
                }
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return;

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

template<int n_tpt, int min_ept, typename Index, class Func, class ...Args>
void threads::parallel_for(Index begin, Index end, Func &&func, Args &&... args) {
    if (threads::disabled()) {
        // Последовательное выполнение
        for (Index i = begin; i < end; ++i) {
            func(i, std::forward<Args>(args)...);
        }
        return;
    }

#ifdef ZEPHYR_TBB
    tbb::parallel_for(begin, end,
            [&func, &args...](Index idx) {
                func(idx, std::forward<Args>(args)...);
            });
#else
    auto bin_function =
            [&func, &args...](Index from, Index to) {
                for (Index idx = from; idx < to; ++idx) {
                    func(idx, std::forward<Args>(args)...);
                }
            };

    Index size = end - begin;

    // Пустой диапазон
    if (size < 1) return;

    int n_tasks = get_n_tasks<n_tpt, min_ept>(n_threads, size);
    size_t bin = size / n_tasks;
    std::vector<std::future<void>> results;
    results.reserve(n_tasks);

    Index from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));

    for (auto &result : results)
        result.wait();
#endif
}

template<int n_tpt, int min_ept, typename Iter, typename Value, class Func, class ...Args>
Value threads::min(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    if (threads::disabled()) { // Последовательное выполнение
        Value res{init};
        for (auto it = begin; it < end; ++it) {
            Value temp = func(*it, std::forward<Args>(args)...);
            if (temp < res) { res = temp; }
        }
        return res;
    }

#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto&& a, auto&& b) { return a < b ? a : b; },
        [&func, &args...](auto&& elem) -> Value {
            return func(elem, std::forward<Args>(args)...);
        });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                Value temp{init};
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp < res) { res = temp; }
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

    Value res{init};
    Value temp{init};
    for (auto &result : results) {
        temp = result.get();
        if (temp < res) { res = temp; }
    }
    return res;
#endif
}

template<int n_tpt, int min_ept, typename Iter, typename Value, class Func, class ...Args>
Value threads::max(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    if (threads::disabled()) { // Последовательное выполнение
        Value res(init);
        for (auto it = begin; it < end; ++it) {
            Value temp = func(*it, std::forward<Args>(args)...);
            if (temp > res) { res = temp; }
        }
        return res;
    }

#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto&& a, auto&& b) { return a > b ? a : b; },
        [&func, &args...](auto&& elem) -> Value {
            return func(elem, std::forward<Args>(args)...);
        });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                Value temp{init};
                for (auto it = a; it < b; ++it) {
                    temp = func(*it, std::forward<Args>(args)...);
                    if (temp > res) { res = temp; }
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

    Value res{init};
    Value temp{init};
    for (auto &result : results) {
        temp = result.get();
        if (temp > res) { res = temp; }
    }
    return res;
#endif
}

template<int n_tpt, int min_ept, typename Iter, typename Value, class Func, class ...Args>
Value threads::sum(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    if (threads::disabled()) { // Последовательное выполнение
        Value res(init);
        for (auto it = begin; it < end; ++it) {
            res += func(*it, std::forward<Args>(args)...);
        }
        return res;
    }

#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto&& a, auto&& b) { auto c = a; c += b; return c; },
        [&func, &args...](auto&& elem) {
            return func(elem, std::forward<Args>(args)...);
        });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                for (auto it = a; it < b; ++it) {
                    res += func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

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

    Value res{init};
    for (auto &result : results) {
        res += result.get();
    }
    return res;
#endif
}

template<int n_tpt, int min_ept, typename Iter, typename Value, class Func, class ...Args>
Value threads::reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args) {
    if (threads::disabled()) { // Последовательное выполнение
        Value res(init);
        for (auto it = begin; it < end; ++it) {
            res &= func(*it, std::forward<Args>(args)...);
        }
        return res;
    }

#ifdef ZEPHYR_TBB
    return std::transform_reduce(std::execution::par,
        begin, end, init,
        [](auto&& a, auto&& b) { auto c = a; c &= b; return c; },
        [&func, &args...](auto&& elem) {
            return func(elem, std::forward<Args>(args)...);
        });
#else
    auto bin_function =
            [&init, &func, &args...](const Iter &a, const Iter &b) -> Value {
                Value res{init};
                for (auto it = a; it < b; ++it) {
                    res &= func(*it, std::forward<Args>(args)...);
                }
                return res;
            };

    int size = end - begin;

    // Пустой диапазон
    if (size < 1) return init;

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

    Value res{init};
    for (auto &result : results) {
        res &= result.get();
    }
    return res;
#endif
}

} // namespace zephyr::utils