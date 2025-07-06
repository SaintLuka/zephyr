#pragma once

#include <vector>

#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/eu_prim.h>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/mpi.h>

#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/mesh/euler/migration.h>

#include <zephyr/mesh/decomp/ORB.h>

namespace zephyr::geom {
class Line;
class Polygon;
class Polyhedron;
class Generator;
}

namespace zephyr::mesh {

class EuMesh;
class AmrCells;
class EuCell;

using utils::threads;
using utils::Storable;

using zephyr::utils::threads;
using zephyr::utils::mpi;

using zephyr::mesh::decomp::Decomposition;
using zephyr::mesh::decomp::ORB;

class EuMesh {
public:

    /// @brief Создание сетки сеточным генератором
    explicit EuMesh(geom::Generator& gen);

    /// @brief Конструируемая сетка, можно добавлять полигоны, но связи
    /// с соседями не восстанавливаются, используется для визуализации.
    EuMesh(int dim, bool adaptive, bool axial = false)
        : m_locals(dim, adaptive, axial) {
    }

    /// @brief Установить максимальный допустимый уровень адаптации (>= 0)
    void set_max_level(int max_level);

    /// @brief Допустимый уровень адаптации
    int max_level() const;

    /// @brief Адаптивная сетка?
    bool adaptive() const;

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

    EuCellIt begin();

    EuCellIt end();

    EuCell operator[](index_t cell_idx) {
        return {&m_locals, cell_idx, &m_aliens};
    }


    template <typename T>
    Storable<T> add(const std::string& name) {
        auto res1 = m_locals.data.add<T>(name);
        auto res2 = m_aliens.data.add<T>(name);
        if (res1.idx != res2.idx) {
            throw std::runtime_error("EuMesh error: bad add<T> #1");
        }
#ifdef ZEPHYR_MPI
        auto res3 = m_tourists.m_border.data.add<T>(name);
        if (res1.idx != res3.idx) {
            throw std::runtime_error("EuMesh error: bad add<T> #2");
        }
        auto res4 = m_migrants.migrants.data.add<T>(name);
        if (res1.idx != res4.idx) {
            throw std::runtime_error("EuMesh error: bad add<T> #3");
        }
#endif
        return res1;
    }

    template <typename T, typename U>
    Storable<T> add(const char*& name) {
        return add<T>(std::string(name));
    }

    // Основная функция, которая принимает произвольное количество строк и возвращает кортеж double
    template <typename T, typename... Args>
    auto append(Args... args) {
        return std::make_tuple(std::invoke(&EuMesh::add<T>, *this, args)...);
    }

    template <typename T>
    void swap(Storable<T> val1, Storable<T> val2) {
        m_locals.data.swap<T>(val1, val2);
        m_aliens.data.swap<T>(val1, val2);
#ifdef ZEPHYR_MPI
        m_tourists.m_border.data.swap<T>(val1, val2);
        m_migrants.migrants.data.swap<T>(val1, val2);
#endif
    }

    void push_back(const geom::Line& line);

    void push_back(const geom::Polygon& poly);

    void push_back(const geom::Polyhedron& poly);


    AmrCells& locals() { return m_locals; }

    const AmrCells& locals() const { return m_locals; }


    /// @brief Выставяет новые ранги ячейкам
    void setup_ranks();

    /// @brief Перераспределить ячейки между процессами, в соответствии
    /// со значениями, которые выдает m_decomp
    template <typename T>
    void migrate(Storable<T> var) {
#ifdef ZEPHYR_MPI
        if (mpi::single()) return;

        setup_ranks();
        m_migrants.migrate(m_tourists, m_locals, m_aliens, var);
#endif
    }

    /// @brief Собрать обменные слои (aliens и сопутствующие члены),
    /// принимаются актуальные данные с border слоя.
    void build_aliens();

    /// @brief Получить ячейку по нескольким индексам подразумевая,
    /// что сетка является структурированной. Индексы периодически
    /// замкнуты (допускаются отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    EuCell operator()(int i, int j);

    /// @brief Получить ячейку по нескольким индексам подразумевая,
    /// что сетка является структурированной. Индексы периодически
    /// замкнуты (допускаются отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    EuCell operator()(int i, int j, int k);


    /// @brief Ссылка на декомпозицию
    const Decomposition& decomp() const { return *m_decomp; }

    /// @brief Добавить декомпозицию сетки, основная функция
    /// @param decmp Умный указатель на декомпозицию
    /// @param update Сразу перераспределить ячейки
    void set_decomposition(Decomposition::Ref decmp, bool update=true);

    /// @brief Добавить ORB декомпозицию сетки
    /// @param orb Ссылка на ORB декомпозицию, внутри функции заменяется
    /// @param update Сразу перераспределить ячейки
    void set_decomposition(ORB& orb, bool update=true);

    /// @brief Добавить ORB декомпозицию сетки, используется простейший
    /// конструктор ORB декомпозиции, ячейки сразу перераспределяются.
    /// @param type Тип ORB декомпозиции
    void set_decomposition(const std::string& type);


    /// @brief Перераспределить ячейки между процессами в соответствии с рангом,
    /// который выдает функция m_decomp::rank().
    /// До вызова redistribute распределенная сетка должна быть согласована
    /// и после вызова остается согласованной (массивы locals и aliens
    /// корректно связаны).
    template <typename T>
    void redistribute(Storable<T> var) {
#ifdef ZEPHYR_MPI
        if (mpi::single()) return;

        migrate(var);
        build_aliens();
#endif
    }

    /// @brief Обмен данными между процессами, в массивы aliens записываются
    /// данные с других процессов. Последовательное выволнение send и recv
    template <typename T>
    void sync(Storable<T> var) {
#ifdef ZEPHYR_MPI
        if (mpi::single()) return;

        mpi::for_each([]() {
            std::cout << "START SYNC " << mpi::rank() << "\n";
        });
        m_tourists.sync(m_locals, m_aliens, var);
        mpi::for_each([]() {
            std::cout << "SUCCESS SYNC " << mpi::rank() << "\n";
        });
#endif
    }

    /// @brief Предбалансировка по числу ячеек
    /// @param n_iters Количество итераций балансировки
    void prebalancing(int n_iters);

    /// @brief Балансирует нагрузку по числу ячеек
    void balancing();

    /// @brief Балансирует нагрузку согласно decomposition
    void balancing(double load);

    /// @brief Дисбаланс нагрузки
    double get_imbalance(const std::vector<double>& ws) const {
        if(m_decomp)
            return m_decomp->imbalance(ws);
        return -1;
    }


    geom::Box bbox() const;


    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, EuCell>>
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


public:
    int m_max_level = 0;
    Distributor distributor;

    AmrCells m_locals;
    AmrCells m_aliens;

    /// @brief Метод декомпозиции
    Decomposition::Ptr m_decomp = nullptr;

#ifdef ZEPHYR_MPI
    Tourism   m_tourists;  ///< Построение обменных слоев и обмены
    Migration m_migrants;  ///< Пересылка ячеек при изменении декомпозиции
#endif

    /// @brief Структура сетки, если предполагается, что сетка декартова.
    bool structured = false;
    int m_nx = 1, m_ny = 1, m_nz = 1;
};

}