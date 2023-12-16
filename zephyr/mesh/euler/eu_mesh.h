#pragma once

#include <zephyr/utils/threads.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bnodes.h>

#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/eu_face.h>
#include <zephyr/mesh/euler/eu_cell.h>

namespace zephyr::mesh {

using zephyr::utils::threads;
using namespace zephyr::geom;

class EuMesh {
public:

    template<class T>
    EuMesh(const T &val, Generator *gen)
            : m_locals(val), m_aliens(val) {
        initialize(gen->make());
    }

    template <class T>
    EuMesh(const T&val, const Grid& grid)
            : m_locals(val), m_aliens(val) {
        initialize(grid);
    }

    inline int n_cells() { return m_locals.size(); }

    EuCell begin() { return {m_locals, m_aliens, 0}; }

    EuCell end() { return {m_locals, m_aliens, m_locals.size()}; }

    operator AmrStorage &() { return m_locals; }

    AmrStorage &locals() { return m_locals; }

    AmrStorage &aliens() { return m_aliens; }


    /// @brief Проверить базовую сетку после создания.
    /// @return -1, если сетка не подходит для адаптации.
    int check_base();

    /// @brief Проверить сетку после адаптации (для дебага).
    /// @return -1, если сетка имеет неверную структуру.
    int check_refined();

    /// @brief Адаптивная сетка?
    bool is_adaptive() const;

    /// @brief Максимальный допустимый уровень адаптации (>= 0)
    int max_level() const;

    /// @brief Установить максимальный допустимый уровень адаптации (>= 0)
    void set_max_level(int max_level);

    /// @brief Установить распределитель данных при адаптации
    void set_distributor(Distributor distr);

    /// @brief Основная функция адаптации, меняет хранилище в соответствии
    /// с флагами amr.flag ячеек.
    void refine();

    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    template<int n_tasks_per_thread = 10, class Func,
            class Value = typename std::result_of<Func(EuCell&)>::type>
    auto min(const Value &init, Func &&func)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    template<int n_tasks_per_thread = 10, class Func,
            class Value = typename std::result_of<Func(EuCell&)>::type>
    auto min(Func &&func)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), std::forward<Func>(func));
    }

    /// @brief Найти описывающий параллелепипед (Bounding box)
    geom::Box bbox();

    /// @brief Заполнен ли массив с уникальными узлами?
    bool has_nodes() const;

    /// @brief Собрать массивы уникальных узлов
    void collect_nodes();

    /// @brief Ссылка на массив уникальных узлов
    const std::vector<Vector3d>& nodes() const {
        return m_nodes;
    }

private:
    /// @brief Создать эйлерову сетку из сетки общего вида
    void initialize(const Grid& gen);

    /// @brief Очистить массив узлов, функция вызывается при любом
    /// перестроении сетки (вроде адаптации)
    void break_nodes();

    /// @brief Осуществляет инициализацию хранилища перед использованием
    /// функций адаптации, выполняется один раз после создания хранилища.
    void init_amr();

    int m_max_level = 0;
    Distributor distributor;

    AmrStorage m_locals;
    AmrStorage m_aliens;

    /// @brief Массив уникальных узлов. Используется в редких алгоритмах,
    /// по умолчанию пустой, заполняются при вызове функции collect_nodes().
    std::vector<Vector3d> m_nodes;
};


} // namespace zephyr::mesh