#pragma once

#include <string>
#include <chrono>
#include <vector>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/numpy.h>

namespace zephyr::utils {

/// @brief Секундомер
class Stopwatch {
public:
    using clock = std::chrono::high_resolution_clock;
    using duration = clock::duration;
    using time_point = clock::time_point;

	/// @brief Используется в MPI-расчетах, собирает измерения со всех процессов
	template<typename T>
	struct multiple {
		std::vector<T> m_times; ///< Набор измерений со всех процессов

		/// @brief Конструктуор собирает со всех
		explicit multiple(T t) : m_times(mpi::all_gather(t)) { }

		/// @brief Минимальное время
		T min() const { return np::min(m_times); }

		/// @brief Максимальное время
		T max() const { return np::max(m_times); }

		/// @brief Среднее время
		T avg() const { return np::mean(m_times); }
	};

	/// @brief Конструктор класса
	/// @param run Включить секундомер сразу
	explicit Stopwatch(bool run = false);

	/// @brief Измерить время выполнения функции, возвращающей void
	/// @param f Целевая функция
	/// @param args Аргументы функции
	template<class F, class... Args>
	std::enable_if_t<std::is_void_v<std::result_of_t<F(Args...)>>, void>
	measure(F&& f, Args&&... args);

    /// @brief Измерить время выполнения функции, возвращающей не void
    /// @param f Целевая функция
    /// @param args Аргументы функции
	template<class F, class... Args>
	std::enable_if_t<!std::is_void_v<std::result_of_t<F(Args...)>>, std::result_of_t<F(Args...)>>
	measure(F&& f, Args&&... args);

	/// @brief Включить секундомер (со сбросом времени)
    void start();

    /// @brief Выключить секундомер с обнулением
    void reset();

	/// @brief Остановить измерение времени
    void stop();

	/// @brief Продолжить измерение времени (включить без сброса)
    void resume();

	/// @brief Проверить включен ли секундомер (измеряет ли время)
    bool is_up() const;

    /// @brief Количество миллисекунд
	long milliseconds() const;

	/// @brief Количество миллисекунд
	multiple<long> milliseconds_mpi() const;

    /// @brief Количество секунд
	long seconds() const;

	/// @brief Количество секунд в формате float
	float fseconds() const;

	/// @brief Количество секунд в формате float
	multiple<float> fseconds_mpi() const;

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
std::enable_if_t<std::is_void_v<std::result_of_t<F(Args...)>>, void>
Stopwatch::measure(F&& f, Args&&... args) {
    start();
    std::forward<F>(f)(std::forward<Args>(args)...);
    stop();
}

template<class F, class... Args>
std::enable_if_t<!std::is_void_v<std::result_of_t<F(Args...)>>, std::result_of_t<F(Args...)>>
Stopwatch::measure(F&& f, Args&&... args) {
    start();
    auto res = std::forward<F>(f)(std::forward<Args>(args)...);
    stop();
    return res;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Stopwatch::multiple<T>& m ) {
	os << m.max();
	if (m.m_times.size() > 1) {
		os << " [" << m.min() << ", " << m.avg() << "]";
	}
	return os;
}

} // namespace zephyr::utils