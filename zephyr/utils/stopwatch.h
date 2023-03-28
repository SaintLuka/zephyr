#pragma once
#include <utility>
#include <iomanip>
#include <sstream>
#include <zephyr/performance/timer/times.h>

namespace zephyr { namespace performance { namespace timer {


/** \brief
	\~russian Класс секундомера. Измеряет время. 
 	\~english Stopwatch class. Measures time.
	\~
*/
class Stopwatch {
public:
	/** \brief
			\~russian Конструктор класса.
			\~english Class constructor.
			\~
		\details
			\~russian Создает объект секундомера.
				\param run Включить секундомер сразу
			\~english Creates a \p Stopwatch object.
			    \param run Turn on Stopwatch right away.
			\~
	*/
	explicit Stopwatch(bool run = false)
	    : start_time(), elapsed(), up(false) {
	    if (run) {
	        start();
	    }
	}

	/** \brief
		\~russian Измерить время выполнения функции, возвращающей \p void.
		\~english Measure time spent on function execution for a function 
			returning \p void.
		\~
		\param[in] f
			\~russian функция.
			\~english a function.
			\~
		\param[in] args
			\~russian аргументы функции.
			\~english arguments of a function.
			\~
	*/
	template<class F, class... Args>
	typename std::enable_if<std::is_void<typename std::result_of<F(Args...)>::type>::value>::type
	measure(F&& f, Args&&... args) {
    	start(); 
		std::forward<F>(f)(std::forward<Args>(args)...);
		stop();
	}

	/** \brief
		\~russian Измерить время выполнения функции, возвращающей не \p void.
		\~english Measure time spent on function execution for a function 
			returning non-\p void.
		\~
		\param[in] f
			\~russian функция.
			\~english a function.
			\~
		\param[in] args
			\~russian аргументы функции.
			\~english arguments of a function.
			\~
	*/
	template<class F, class... Args>
	typename std::enable_if<!std::is_void<typename std::result_of<F(Args...)>::type>::value, typename std::result_of<F(Args...)>::type>::type
	measure(F&& f, Args&&... args) {
		start(); 
		auto res = std::forward<F>(f)(std::forward<Args>(args)...);
		stop();
		return res;
	}

	/** \brief
		\~russian Начать измерение времени (включить секундомер).
		\~english Start measuring time (turn stopwatch on).
		\~
	*/	
    void start() { 
    	up = true;
    	elapsed = elapsed_time_t();
    	start_time = now();
    }

	/** \brief
		\~russian Остановить измерение времени (выключить секундомер).
		\~english Stop measuring time (turn stopwatch off).
		\~
	*/	
    void stop() {
    	if (up) {
	    	elapsed += now_proc_and_wall() - start_time;
	    	up = false;
	    }
    }
	/** \brief
		\~russian Продолжить измерение времени (включить секундомер без сброса).
		\~english Resume measuring time (turn stopwatch on without 
			starting over).
		\~
	*/	
    void resume () {
    	if (!up) {
		    start_time = now();
	    	up = true;
	    }
    }

	/** \brief
		\~russian Проверить, включен ли секундомер (измеряет ли он время).
		\~english Check if a stopwatch is on, i.e. measures time at the moment.
		\~
	*/	
    bool is_up() const { return up; }

	/** \brief
		\~russian Получить показания секундомера.
		\~english Get times measured by a stopwatch.
		\~
	*/	
	template<class T = fseconds>
    elapsed_time_t times() const { 
    	if (up) 
    		return elapsed + (now_proc_and_wall() - start_time);
    	else    
    		return elapsed;
    }

    /// @brief Расширенный формат времени с указанием количества дней,
    /// часов, минут и секунд
    std::string expended_time_format() const {
        long full = static_cast<long>(times().wall());

        long minutes = full / 60   % 60;
        long hours   = full / 3600 % 24;
        long days    = full / 86400;
        long seconds = full % 60;

        std::stringstream ss;
        if (days > 0) {
            ss << std::setw(2) << days << " d ";
        }
        else {
            ss << "     ";
        }
        if (hours > 0 || days > 0) {
            ss << std::setw(2) << hours << " h ";
        }
        else {
            ss << "     ";
        }
        ss << std::setw(2) << minutes << " m ";
        ss << std::setw(2) << seconds << " s";

        return ss.str();
    }

    /// @brief Количество миллисекунд приведенное к целому числу
    size_t milliseconds() const {
        return static_cast<size_t>(1000.0 * times().wall());
    }

private:
	elapsed_time_t elapsed;
	/**<\~russian Прошедшее время.
		\~english Elapsed time.
		\~
	*/	

	time_point_t start_time;
	/**<\~russian Время запуска секундомера.
		\~english Start time.
		\~
	*/	

    bool up;
    /**<\~russian Переменная-флаг, показывает, включен ли секундомер.
		\~english Variable that indicates if a stopwatch is on.
		\~
	*/	
};

static Stopwatch sw_dummy;
    /**<\~russian Секундомер по умолчанию.
		\~english A default dummy stopwatch.
		\~
	*/	
	
} /* timer */ } /* performance */ } /* zephyr */