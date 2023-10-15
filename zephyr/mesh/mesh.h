#pragma once

#include <zephyr/geom/generator/generator.h>
#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/distributor.h>

#include <zephyr/mesh/range.h>
#include <zephyr/utils/threads.h>



namespace zephyr { namespace mesh {

#define EXEC_RESULT typename std::result_of<F(ICell)>::type

using zephyr::utils::threads;
using namespace zephyr::geom;

class Mesh {
public:

    template<class T>
    Mesh(const T &val, Generator *gen)
            : m_locals(val), m_aliens(val) {
        initialize(gen->create());
    }

    template <class T>
    Mesh(const T&val, const Grid& grid)
            : m_locals(val), m_aliens(val) {
        initialize(grid);
    }


    ICell begin() { return {m_locals, m_aliens, 0}; }

    ICell end() { return {m_locals, m_aliens, m_locals.size()}; }

    operator Storage &() { return m_locals; }

    Storage &locals() { return m_locals; }

    Storage &aliens() { return m_aliens; }


    /// @brief Проверить базовую сетку после создания.
    /// @return -1, если сетка не подходит для адаптации.
    int check_base();

    /// @brief Проверить сетку после адаптации (для дебага).
    /// @return -1, если сетка имеет неверную структуру.
    int check_refined();

    /// @brief Максимальный допустимый уровень адаптации (>= 0)
    int max_level() const;

    /// @brief Установить максимальный допустимый уровень адаптации (>= 0)
    void set_max_level(int max_level);

    /// @brief Установить распределитель данных при адаптации
    void set_distributor(Distributor distr);

    /// @brief Основная функция адаптации, меняет хранилище в соответствии
    /// с флагами amr.flag ячеек.
    void refine();

    template<int n_tasks_per_thread = 1, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    template<int n_tasks_per_thread = 1, class Func,
            class Value = typename std::result_of<Func(ICell)>::type>
    auto min(const Value &init, Func &&func)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    template<int n_tasks_per_thread = 1, class Func,
            class Value = typename std::result_of<Func(ICell)>::type>
    auto min(Func &&func)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), std::forward<Func>(func));
    }

private:

    void initialize(const Grid& gen);

    /// @brief Осуществляет инициализацию хранилища перед использованием
    /// функций адаптации, выполняется один раз после создания хранилища.
    void init_amr();

    int m_max_level = 0;
    Distributor distributor;

    Storage m_locals;
    Storage m_aliens;
};


} // mesh
} // zephyr