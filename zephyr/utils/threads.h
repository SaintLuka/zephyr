#pragma once

#include <iostream>
#include <zephyr/utils/thread-pool.h>

namespace zephyr::utils {


/// @brief Статический класс. Упрощенный интерфейс для многопоточности.
class threads {
public:

    /// @brief Включить многопоточность, будут использованы все
    /// доступные ядра.
    static void on();

    /// @brief Включить треды
    /// @param count Число тредов.
    static void on(int count);

    /// @brief Выключить треды
    static void off();

    /// @brief Многопоточность включена?
    static bool is_on();

    /// @brief Многопоточность выключена?
    static bool is_off();

    /// @brief Число тредов
    static int count();

    /// @brief Выполнить функцию для элементов из диапазона
    /// @param begin Итератор, указывающий на начало диапазона
    /// @param end Итератор, указывающий за последний элемент диапазона
    /// @param func Целевая функция, принимает аргументы (*Iter, Args...)
    /// @param args Аргуметры функции
    /// @details Целевая функция func в качестве аргументов принимает
    /// разыменованный итератор Iter и набор аргументов Args..., целевая
    /// функция может иметь возвращаемое значение, но оно игнорируется.
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args>
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
    static auto min(Iter begin, Iter end, Func &&func, Args &&... args)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return min<n_tasks_per_thread>(begin, end, std::numeric_limits<Value>::max(),
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
    static auto max(Iter begin, Iter end, Func &&func, Args &&... args)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return max<n_tasks_per_thread>(begin, end, std::numeric_limits<Value>::lowest(),
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
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
    template<int n_tasks_per_thread = 10, class Iter, class Func, class... Args,
            class DeRef = decltype(*std::declval<Iter&>()),
            class Value = typename std::result_of<Func(DeRef, Args...)>::type>
    static auto reduce(Iter begin, Iter end, const Value &init, Func &&func, Args &&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type;

public:
    /// @brief Число тредов
    static int n_threads;

    /// @brief Указатель на пул тредов
    static std::unique_ptr<ThreadPool> pool;
};


template<int n_tasks_per_thread, class Iter, class Func, class ...Args>
void threads::for_each(Iter begin, Iter end, Func&& func, Args&&... args) {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
    std::vector<std::future<void>> results;
    results.reserve(n_tasks);

#if 0
    // Тестировал вариант с непоследовательной загрузкой задач
    for (int j = 0; j < n_tasks_per_thread; ++j) {
        for (int i = 0; i < n_threads; ++i) {
            int a = (i * n_tasks_per_thread + j) * bin;
            int b = a + bin;
            if (i * n_tasks_per_thread + j == n_tasks - 1) {
                b = size;
            }
            results.emplace_back(pool->enqueue(bin_function, begin + a, begin + b));
        }
    }
#else
    Iter from = begin;
    for (int i = 0; i < n_tasks - 1; ++i) {
        results.emplace_back(pool->enqueue(bin_function, from, from + bin));
        from += bin;
    }
    results.emplace_back(pool->enqueue(bin_function, from, end));
#endif

    for (auto &result : results)
        result.wait();
}

template<int n_tasks_per_thread, class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::min(Iter begin, Iter end, const Value& init, Func&& func, Args&&... args)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
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
}

template<int n_tasks_per_thread, class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::max(Iter begin, Iter end, const Value& init, Func&& func, Args&&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
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
}

template<int n_tasks_per_thread, class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::partial_sum(Iter begin, Iter end, const Value& init, Func&& func, Args&&... args)
-> typename std::enable_if<!std::is_void<Value>::value, std::vector<Value>>::type {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
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
    for (auto& r: results) {
        res.push_back(r.get());
    }

    return res;
}

template<int n_tasks_per_thread, class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::sum(Iter begin, Iter end, const Value& init, Func&& func, Args&&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
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
}

template<int n_tasks_per_thread, class Iter, class Func, class ...Args, class DeRef, class Value>
auto threads::reduce(Iter begin, Iter end, const Value& init, Func&& func, Args&&... args)
-> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
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

    int n_tasks = n_tasks_per_thread * n_threads;
    int bin = size / n_tasks;
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
}

} // zephyr