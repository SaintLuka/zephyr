#pragma once

#include <string>
#include <sstream>
#include <chrono>

namespace zephyr { namespace utils {

/// @brief Секуномер
class Stopwatch {
public:
    using clock = std::chrono::high_resolution_clock;
    using duration = clock::duration;
    using time_point = clock::time_point;

	/// @brief Конструктор класса
	/// @param run Включить секундомер сразу
	explicit Stopwatch(bool run = false);

	/// @brief Измерить время выполнения функции, возвращающей void
	/// @param f Целевая функция
	/// @param args Аругменты функции
	template<class F, class... Args>
	typename std::enable_if<std::is_void<typename std::result_of<F(Args...)>::type>::value>::type
	measure(F&& f, Args&&... args);

    /// @brief Измерить время выполнения функции, возвращающей не void
    /// @param f Целевая функция
    /// @param args Аругменты функции
	template<class F, class... Args>
	typename std::enable_if<!std::is_void<typename std::result_of<F(Args...)>::type>::value,
	typename std::result_of<F(Args...)>::type>::type
	measure(F&& f, Args&&... args);

	/// @brief Включить секундомер (со сбросом времени)
    void start();

	/// @brief Остановить измерение времени
    void stop();

	/// @brief Продолжить измерение времени (включить без сброса)
    void resume();

	/// @brief Проверить включен ли секундомер (измеряет ли время)
    bool is_up() const;

    /// @brief Количество миллисекунд
    long milliseconds() const;

    /// @brief Количество секунд
    long seconds() const;

    /// @brief Количество минут
    long minutes() const;

    /// @brief Количество часов
    long hours() const;

    /// @brief Количество дней
    long days() const;

    /// @brief Расширенный формат времени с указанием количества дней,
    /// часов, минут и секунд
    std::string extended_time() const;

private:

    /// @brief Текущие показания
    duration elapsed() const;

    bool m_up;           ///< Включен ли секундомер?
	time_point m_start;  ///< Время последнего запуска
    duration m_elapsed;  ///< Прошлые замеры
};

template<class F, class... Args>
typename std::enable_if<std::is_void<typename std::result_of<F(Args...)>::type>::value>::type
Stopwatch::measure(F&& f, Args&&... args) {
    start();
    std::forward<F>(f)(std::forward<Args>(args)...);
    stop();
}

template<class F, class... Args>
typename std::enable_if<!std::is_void<typename std::result_of<F(Args...)>::type>::value,
typename std::result_of<F(Args...)>::type>::type
Stopwatch::measure(F&& f, Args&&... args) {
    start();
    auto res = std::forward<F>(f)(std::forward<Args>(args)...);
    stop();
    return res;
}
	
} // utils
} // zephyr