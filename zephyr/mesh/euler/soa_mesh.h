#pragma once

#include <vector>

#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/soa_prim.h>

#include <zephyr/utils/threads.h>

namespace zephyr::geom {
class Line;
class Polygon;
class Polyhedron;
class Generator;
}

namespace zephyr::mesh {

class SoaMesh;
class AmrCells;
class QCell;

using utils::threads;
using utils::Storable;

class SoaMesh {
public:
    AmrCells m_locals;
    AmrCells m_aliens;

    int m_max_level = 0;
    Distributor distributor;


    /// @brief Сетка как набор полигонов
    SoaMesh(int dim, bool adaptive, bool axial = false)
        : m_locals(dim, adaptive, axial),
          m_aliens(dim, adaptive, axial) {
    }

    /// @brief Создание сеточным генератором
    SoaMesh(geom::Generator& gen);


    /// @brief Установить максимальный допустимый уровень адаптации (>= 0)
    void set_max_level(int max_level);

    int max_level() const;

    bool is_adaptive() const;

    /// @brief Установить распределитель данных при адаптации,
    /// допустимые значения: "empty", "simple".
    void set_distributor(const std::string& name);

    /// @brief Установить распределитель данных при адаптации
    void set_distributor(Distributor distr);

    void init_amr();

    int dim() const { return m_locals.dim(); }

    bool axial() const { return m_locals.axial(); }

    void balance_flags();

    void apply_flags();

    void refine();

    int check_base() const;

    int check_refined() const;

    CellIt begin();

    CellIt end();

    QCell operator[](index_t cell_idx) {
        return {&m_locals, cell_idx, &m_aliens};
    }


    template <typename T>
    Storable<T> add(const std::string& name) {
        auto res1 = m_locals.data.add<T>(std::string(name));
        auto res2 = m_aliens.data.add<T>(std::string(name));
        if (res1.idx != res2.idx) {
            throw std::runtime_error("EuMesh error: bad add<T>");
        }
        return res1;
    }

    template <typename T, typename U>
    Storable<T> add(const char*& name) {
        return add<T>(std::string(name));
    }

    // Вспомогательная функция для применения add<T> к каждому элементу
    template <typename T, typename... Args>
    auto add_multiply(Args... args) {
        return std::make_tuple(std::invoke(&SoaMesh::add<T>, *this, args)...);
    }

    // Основная функция, которая принимает произвольное количество строк и возвращает кортеж double
    template <typename T, typename... Args>
    auto append(Args... args) {
        return add_multiply<T>(args...);
    }

    template <typename T>
    void swap(const Storable<T>& val1, const Storable<T>& val2) {
        m_locals.data.swap<T>(val1, val2);
        m_aliens.data.swap<T>(val1, val2);
    }

    void push_back(const geom::Line& line);

    void push_back(const geom::Polygon& poly);

    void push_back(const geom::Polyhedron& poly);


    AmrCells& locals() { return m_locals; }

    const AmrCells& locals() const { return m_locals; }



    /// @brief Получить ячейку по нескольким индексам подразумевая,
    /// что сетка является структурированной. Индексы периодически
    /// замкнуты (допускаются отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    QCell operator()(int i, int j);

    /// @brief Получить ячейку по нескольким индексам подразумевая,
    /// что сетка является структурированной. Индексы периодически
    /// замкнуты (допускаются отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    QCell operator()(int i, int j, int k);




    geom::Box bbox() const;


    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, QCell>>
    auto min(Func &&func, const Value &init)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @brief Число ячеек по оси x для структурированных сеток,
    /// число всех ячеек для сеток общего вида
    int nx() const { return m_nx; };

    /// @brief Число ячеек по оси y для структурированных сеток,
    /// единица для сеток общего вида
    int ny() const { return m_ny; };

    /// @brief Число ячеек по оси z для структурированных сеток,
    /// единица для сеток общего вида
    int nz() const { return m_nz; };

    /// @brief Структура сетки, если предполагается, что сетка декартова.
    bool structured = false;
    int m_nx = 1, m_ny = 1, m_nz = 1;
};

}