#pragma once

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

namespace zephyr { namespace utils {

#define THREADS_RES_IT typename std::result_of<F(iterator, Args...)>::type
#define THREADS_RES    typename std::result_of<F(std::size_t, Args...)>::type

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
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает void)
    /// @param args Аргуметры функции
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto for_each(iterator a, iterator b, F&& f, Args&&... args)
    -> typename std::enable_if<std::is_void<THREADS_RES_IT>::value, void>::type;

    /// @brief Выполнить функцию для элементов из диапазона
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает не void)
    /// @param args Аргуметры функции
    /// @return Массив результатов выполнения функции
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto for_each(iterator a, iterator b, F&& f, Args&&... args)
    -> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
            std::vector<THREADS_RES_IT>>::type;

    /// @brief Применить свертку к результатам выполнения функции на каждом 
    /// элементе диапазона. Для возвращаемого типа должен быть определен 
    /// (перегружен) оператор &=
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает тип с оператором &=)
    /// @param args Аргуметры функции
    /// @return Результат свертки (тип возвращаемого значения функции f)
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto reduce(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
    -> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
            THREADS_RES_IT>::type;

    /// @brief Найти минимальное значение среди результатов выполнения функции 
    /// для каждого элемента диапазонана. Для возвращаемого типа должен быть 
    /// определен (перегружен) оператор <
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает тип с оператором <)
    /// @param args Аргуметры функции
    /// @return Минимальное значение (тип возвращаемого значения функции f)
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto min(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
    -> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
            THREADS_RES_IT>::type;

    /// @brief Найти максимимальное значение среди результатов выполнения функции 
    /// для каждого элемента диапазонана. Для возвращаемого типа должен быть 
    /// определен (перегружен) оператор >
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает тип с оператором >)
    /// @param args Аргуметры функции
    /// @return Максимальное значение (тип возвращаемого значения функции f)
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto max(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
    -> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
            THREADS_RES_IT>::type;

    /// @brief Найти сумму значений результатов выполнения функции для каждого 
    /// элемента диапазонана. Для возвращаемого типа должен быть определен 
    /// (перегружен) оператор +=
    /// @param a Итератор на начало
    /// @param b Итератор конца
    /// @param f Целевая функция (возвращает тип с оператором +=)
    /// @param args Аргуметры функции
    /// @return Сумма значений (тип возвращаемого значения функции f)
    template<std::size_t ntasks_per_thread = 16, class iterator, class F, class... Args>
    static auto sum(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
    -> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
            THREADS_RES_IT>::type;

protected:
    static int n_threads;
};

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::for_each(iterator a, iterator b, F&& f, Args&&... args)
-> typename std::enable_if<std::is_void<THREADS_RES_IT>::value, void>::type {
    std::size_t size = (b - a); if (size == 0) return;
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<void>> results(num_tasks);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a, iterator b){
        for (auto it = a; it != b; ++it)
            f(it, std::forward<Args>(args)...);
    };

    for (auto &result : results) {
        auto to = std::min((iterator&&)(pos + bin), end);
        result = std::move(enqueue(bin_function, pos, to));
        pos = pos + bin;
    }

    for (auto& result : results)
        result.wait();
};

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::for_each(iterator a, iterator b, F&& f, Args&&... args)
-> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
        std::vector<THREADS_RES_IT>>::type {
    std::size_t size = (b - a);
    using result_type = std::vector<THREADS_RES_IT>;

    if (size == 0) return result_type();
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<void>> results_future(num_tasks);
    result_type results(size);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a1, iterator b1) {
        for (auto it = a1; it != b1; ++it)
            results[it - a] = f(it, std::forward<Args>(args)...);
    };

    for (auto &result : results_future) {
        auto to = std::min((iterator&&)(pos + bin), end);
        result = enqueue(bin_function, pos, to);
        pos = pos + bin;
    }

    for (auto& result : results_future)
        result.wait();

    return results;
}

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::reduce(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
-> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
        THREADS_RES_IT>::type {
    std::size_t size = (b - a);
    using result_type = THREADS_RES_IT;

    if (size == 0) return result_type();
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<result_type>> results_future(num_tasks);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a1, iterator b1) -> result_type {
        result_type res = init;
        for (auto it = a1; it != b1; ++it) {
            res &= f(it, std::forward<Args>(args)...);
        }
        return res;
    };

    for (auto &result: results_future) {
        auto to = std::min((iterator &&) (pos + bin), end);
        result = enqueue(bin_function, pos, to);
        pos = pos + bin;
    }

    result_type result = init;
    for (auto &partial_result: results_future) {
        result &= partial_result.get();
    }

    return result;
}

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::min(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
-> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
        THREADS_RES_IT>::type {
    std::size_t size = (b - a);
    using result_type = THREADS_RES_IT;

    if (size == 0) return init;
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<result_type>> results_future(num_tasks);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a1, iterator b1) -> result_type {
        result_type res(init);
        for (auto it = a1; it != b1; ++it) {
            result_type res2 = f(it, std::forward<Args>(args)...);
            if (res2 < res) {
                res = res2;
            }
        }
        return res;
    };

    for (auto &result: results_future) {
        auto to = std::min((iterator &&) (pos + bin), end);
        result = enqueue(bin_function, pos, to);
        pos = pos + bin;
    }

    result_type result(init);
    for (auto &partial_result: results_future) {
        result_type res2 = partial_result.get();
        if (res2 < result) {
            result = res2;
        }
    }

    return result;
}

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::max(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
-> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
        THREADS_RES_IT>::type {
    std::size_t size = (b - a);
    using result_type = THREADS_RES_IT;

    if (size == 0) return init;
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<result_type>> results_future(num_tasks);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a1, iterator b1) -> result_type {
        result_type res(init);
        for (auto it = a1; it != b1; ++it) {
            result_type res2 = f(it, std::forward<Args>(args)...);
            if (res2 > res) {
                res = res2;
            }
        }
        return res;
    };

    for (auto &result: results_future) {
        auto to = std::min((iterator &&) (pos + bin), end);
        result = enqueue(bin_function, pos, to);
        pos = pos + bin;
    }

    result_type result(init);
    for (auto &partial_result: results_future) {
        result_type res2 = partial_result.get();
        if (res2 > result) {
            result = res2;
        }
    }

    return result;
}

template<std::size_t ntasks_per_thread, class iterator, class F, class... Args>
auto threads::sum(iterator a, iterator b, const THREADS_RES_IT& init, F&& f, Args&&... args)
-> typename std::enable_if<!std::is_void<THREADS_RES_IT>::value,
        THREADS_RES_IT>::type {
    std::size_t size = (b - a);
    using result_type = THREADS_RES_IT;

    if (size == 0) return init;
    std::size_t num_tasks = ntasks_per_thread * n_threads;
    std::size_t bin = size / num_tasks + 1;
    num_tasks = std::min(size / bin + 1, num_tasks);
    std::vector<std::future<result_type>> results_future(num_tasks);

    iterator pos = a;
    iterator end = b;

    auto bin_function = [&](iterator a1, iterator b1) -> result_type {
        result_type res(init);
        for (auto it = a1; it != b1; ++it) {
            res += f(it, std::forward<Args>(args)...);
        }
        return res;
    };

    for (auto &result: results_future) {
        auto to = std::min((iterator &&) (pos + bin), end);
        result = enqueue(bin_function, pos, to);
        pos = pos + bin;
    }

    result_type result(init);
    for (auto &partial_result: results_future) {
        result += partial_result.get();
    }

    return result;
}

} // utils
} // zephyr