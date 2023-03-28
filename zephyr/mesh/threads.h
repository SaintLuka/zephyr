#pragma once
#include <memory>
#include <vector>
#include <functional>
#include <iostream>
#include <limits>
#include <future>

#include <zephyr/configuration.h>

#ifdef ZEPHYR_ENABLE_MULTITHREADING
#include <zephyr/multithreading/thread-pool.h>
#endif

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
#include <zephyr/network/mpi/network.h>
#endif

namespace zephyr { namespace mesh { namespace medium {

#ifdef ZEPHYR_ENABLE_MULTITHREADING
using zephyr::multithreading::ThreadPool;
#endif


template<class T>
class ThreadResults {
public:

    /// @brief Конструктор по умолчанию
    ThreadResults() = default;

    /// @brief Выполнить функцию в этом треде
    /// @details Не вызываются никакие thread, async и future, всё выполняется
    /// совершенно обычным способо, а результат заносится в массив.
    template<class F, class... Args>
    void emplace_back(F &&func, Args &&... args) {
        if (!m_futures.empty() || !m_result.empty()) {
            throw std::runtime_error("ThreadResults::emplace_back error #1");
        }
        using T2 = typename std::result_of<F(Args...)>::type;
        static_assert(std::is_same<T, T2>::value, "ThreadResults: Types mismatch");

        m_result.emplace_back(func(std::forward<Args>(args)...));
    }

    /// @brief Добавить future
    void emplace_back(std::future<T> &&res) {
        if (!m_result.empty()) {
            throw std::runtime_error("ThreadResults::emplace_back error #2");
        }
        m_futures.emplace_back(std::move(res));
    }

    /// @brief Дождаться завершения всех тредов
    void wait() {
        for (auto &res: m_futures) {
            res.wait();
        }
    }

    /// @brief Функция reduce, так называемая "свертка списка". Функция высшего
    /// порядка, прозводит преобразование списка результатов выполнения задачи
    /// на тредах к атомарному значению с использованием некоторой функции.
    /// @tparam F Бинарная, коммутативная, ассоциативная операция над типом T
    template<class F>
    typename std::result_of<F(const T &, const T &)>::type
    reduce(F &&func) {
        if (!m_result.empty()) {
            return m_result[0];
        } else if (!m_futures.empty()) {
            T res(m_futures[0].get());
            for (size_t i = 1; i < m_futures.size(); ++i) {
                T res2(m_futures[i].get());
                res = func(res, res2);
            }
            return res;
        } else {
            throw std::runtime_error("ThreadResults::reduce() error: empty ThreadResults");
        }
    }


    /// @brief Минимум по тредам
    /// @details Для сравнения используется оператор < класса
    T min() {
        if (!m_result.empty()) {
            return m_result[0];
        } else if (!m_futures.empty()) {
            T res(m_futures[0].get());
            for (size_t i = 1; i < m_futures.size(); ++i) {
                T res2(m_futures[i].get());
                if (res2 < res) {
                    res = res2;
                }
            }
            return res;
        } else {
            throw std::runtime_error("ThreadResults::min() error: empty ThreadResults");
        }
    }

    /// @brief Максимум по тредам
    /// @details Для сравнения используется оператор > класса
    T max() {
        if (!m_result.empty()) {
            return m_result[0];
        } else if (!m_futures.empty()) {
            T res(m_futures[0].get());
            for (size_t i = 1; i < m_futures.size(); ++i) {
                T res2(m_futures[i].get());
                if (res2 > res) {
                    res = res2;
                }
            }
            return res;
        } else {
            throw std::runtime_error("ThreadResults::max() error: empty ThreadResults");
        }
    }

    /// @brief Сумма по тредам
    /// @details Для суммирования используется оператор += класса
    T sum() {
        if (!m_result.empty()) {
            return m_result[0];
        } else if (!m_futures.empty()) {
            T res(m_futures[0].get());
            for (size_t i = 1; i < m_futures.size(); ++i) {
                T res2(m_futures[i].get());
                res += res2;
            }
            return res;
        } else {
            throw std::runtime_error("ThreadResults::sum() error: empty ThreadResults");
        }
    }

private:
    /// @brief Непосредственный результат обычного (последовательного)
    /// выполнения некоторой функции. Размер массива не больше единицы, то есть
    /// хранится результат только для одного треда.
    std::vector<T> m_result;

    /// @brief Набор future для каждого треда. В программе без тредов данный
    /// массив должен быть пустым. В любом случае один из массивов m_result
    /// или m_futures должен оставаться пустым.
    std::vector<std::future<T>> m_futures;
};


/// @brief Специализация класса для типа void, то есть для параллельного
/// выполнения функций, которые ничего не возвращают. Данная специализация
/// не содержит reduce операций.
template<>
class ThreadResults<void> {
public:

    /// @brief Конструктор по умолчанию
    ThreadResults() = default;

    /// @brief Добавить тред
    void emplace_back(std::future<void> &&res) {
        m_futures.emplace_back(std::move(res));
    }

    template<class F, class... Args>
    void emplace_back(F &&func, Args &&... args) {
        using T2 = typename std::result_of<F(Args...)>::type;
        static_assert(std::is_same<void, T2>::value, "ThreadResults: Types mismatch");

        func(std::forward<Args>(args)...);
    }

    /// @brief Подождать завершения
    void wait() {
        for (auto &res: m_futures) {
            res.wait();
        }
    }

private:
    std::vector<std::future<void>> m_futures;
};

} // namespace medium
} // namespace mesh
} // namespace zephyr