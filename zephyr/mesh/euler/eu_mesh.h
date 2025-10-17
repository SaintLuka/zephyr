#pragma once

#include <zephyr/configuration.h>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/mpi.h>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/tourism.h>
#include <zephyr/mesh/euler/migration.h>
#include <zephyr/mesh/decomp/ORB.h>

// forward declaration
namespace zephyr::geom {
class Line;
class Polygon;
class Polyhedron;
class Generator;
}

namespace zephyr::mesh {

/// @brief Эйлерова сетка.
/// @ingroup euler-mesh
///
/// Поддерживается три типа сеток:
///   1. Двумерная AMR сетка, по 8 граней, по 9 вершин на ячейку.
///   2. Трехмерная AMR сетка, по 24 грани, по 27 вершин на ячейку.
///   3. Неструктурированная/произвольная сетка. Произвольное число граней
///      и вершин на ячейку, но вершины не уникальны.
class EuMesh {
    if_mpi(using mpi = utils::mpi;)
    using threads = utils::threads;
    using ORB = decomp::ORB;
    using Decomposition = decomp::Decomposition;

public:
    /// @{ @name Создание сетки

    /// @brief Создание сетки с помощью сеточного генератора
    explicit EuMesh(geom::Generator& gen);

    /// @brief Конструируемая сетка, можно добавлять полигоны, но связи
    /// с соседями не восстанавливаются, используется для визуализации.
    EuMesh(int dim, bool adaptive, bool axial = false);

    /// @brief Инициализация сетки из json-конфига
    explicit EuMesh(const utils::Json& config);


    /// @brief Добавить на неструктурированную сетку ячейку в виде отрезка
    /// (сплюснутая четырехугольная ячейка)
    void push_back(const geom::Line& line);

    /// @brief Добавить на неструктурированную сетку ячейку в виде
    /// произвольного многоугольника
    void push_back(const geom::Polygon& poly);

    /// @brief Добавить на неструктурированную сетку ячейку в виде
    /// произвольного многогранника
    void push_back(const geom::Polyhedron& poly);

    /// @brief Добавить треугольный маркер (ячейка)
    void add_marker(const geom::Vector3d& pos, double size);

    /// @}

    /// @{ @name Общие свойства

    /// @brief Размерность сетки (= 2 для сеток с осевой симметрией)
    int dim() const { return m_locals.dim(); }

    /// @brief Сетка с осевой симметрией?
    bool axial() const { return m_locals.axial(); }

    /// @brief Отсутствуют ячейки на данном процессе?
    bool empty() const { return n_cells() == 0; }

    /// @brief Число ячеек на данном процессе
    index_t n_cells() const { return m_locals.size(); }

    /// @brief Ограничивающий прямоугольник (кубоид) области
    geom::Box bbox() const;

    /// @}

    /// @{ @name Массивы данных

    /// @brief Добавить несколько массивов данных в хранилище
    /// @param names Имена массивов данных
    /// @return Если передан один аргумент, то возвращает единственный Storable<T>,
    /// при наличии нескольких аргументов возвращает кортеж Storable<T>.
    template<typename T, typename... Names, typename = std::enable_if_t<
            (sizeof...(Names) > 0) && (std::is_convertible_v<Names, std::string> && ...)>>
    auto add(Names&&... names) {
        if constexpr (sizeof...(Names) == 1) {
            return add_one<T>(std::string(std::forward<Names>(names))...);
        } else {
            // Короче жесть, тут сложно было добиться соблюдения порядка
            // с круглыми скобками std::tuple() или std::make_tuple() не работают.
            return std::tuple{add_one<T>(std::string(std::forward<Names>(names)))...};
        }
    }

    /// @brief Поменять местами два массива данных
    template <typename T>
    void swap(Storable<T> var1, Storable<T> var2);

    /// @}

    /// @{ @name Выбор ячеек

    /// @brief Итератор, указывающий на первую ячейку
    EuCell_Iter begin();

    /// @brief Итератор, указывающий на ячейку за последней
    EuCell_Iter end();

    /// @brief Локальная ячейка по индексу
    EuCell operator[](index_t idx);

    /// @}

    /// @{ @name Многопоточное выполнение функций

    /// @brief Выполнить функцию для всех ячеек сетки
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, typename Value, class Func>
    Value min(Func &&func, const Value &init) {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, class Func>
    auto min(Func &&func) {
        return threads::min<n_tasks_per_thread>(begin(), end(), std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать максимум
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, typename Value, class Func>
    Value max(Func &&func, const Value &init) {
        return threads::max<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать максимум
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, class Func>
    auto max(Func &&func) {
        return threads::max<n_tasks_per_thread>(begin(), end(), std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать сумму
    template<int n_tasks_per_thread = utils::default_n_tasks_per_thread, class Func, typename Value>
    auto sum(Func &&func, const Value& init) {
        return threads::sum<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @}

    /// @{ @name Адаптация сетки

    /// @brief Сетка считается адаптивной, если она состоит из AMR-ячеек
    /// и выставлен максимальный уровень адаптации больше нуля.
    bool adaptive() const;

    /// @brief Максимальный уровень адаптации сетки
    int max_level() const;

    /// @brief Установить допустимый уровень адаптации (>= 0).
    /// Не сработает для сеток общего вида без возможности адаптации.
    void set_max_level(int max_level);

    /// @brief Установить распределитель данных при адаптации,
    /// допустимые значения: "empty", "simple".
    void set_distributor(const std::string& name);

    /// @brief Установить распределитель данных при адаптации
    void set_distributor(Distributor distr);

    /// @brief Выполнить адаптацию сетки
    void refine();

    /// @brief Адаптировать сетку целиком до некоторого уровня
    /// По умолчанию до максимального уровня адаптации.
    void refine_full(int level = -1);

    /// @brief Поправить референсную сетку
    void check_reference(bool fix = false);

    /// @}

    /// @{ @name Части распределенной сетки

    /// @brief Неявное преобразование в AmrCells
    operator AmrCells&() { return locals(); }

    /// @brief Локальные ячейки (принадлежат данному процессу)
    AmrCells& locals() { return m_locals; }

    /// @brief Локальные ячейки (принадлежат данному процессу)
    const AmrCells& locals() const { return m_locals; }

    /// @brief Слой обменных ячеек (с других процессов)
    AmrCells& aliens() { return m_aliens; }

    /// @brief Слой обменных ячеек (с других процессов)
    const AmrCells& aliens() const { return m_aliens; }

    /// @}

    /// @{ @name Работа с распределенной сеткой

    /// @brief Обмен данными между процессами, в массивы aliens пересылаются
    /// данные с других процессов. Последовательное выполнение send и recv.
    /// @param vars Положительное количество параметров типа Storable<T>, для
    /// которых осуществляется обмен.
    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) > 0)> >
    void sync(Args&&... vars);

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
    /// конструктор ORB декомпозиции.
    /// @param type Тип ORB декомпозиции
    /// @param update Сразу перераспределить ячейки
    void set_decomposition(const std::string& type, bool update=true);

    /// @brief Балансирует нагрузку по числу ячеек
    void balancing();

    /// @brief Балансирует нагрузку по полученной нагрузке
    void balancing(double load);

    /// @brief Предбалансировка по числу ячеек
    /// @param n_iters Количество итераций балансировки
    void prebalancing(int n_iters);

    /// @brief Перераспределить ячейки между процессами в соответствии с
    /// рангом, который возвращает функция m_decomp::rank().
    ///
    /// До вызова redistribute распределенная сетка должна быть согласована
    /// и после вызова остается согласованной (массивы locals и aliens
    /// корректно связаны). В качестве аргументов передаются параметры,
    /// которые необходимо сохранить и перенести при декомпозиции.
    template <typename... Args>
    void redistribute(Args&&... vars);

    /// @}

    /// @{ @name Функции структурированной сетки

    /// @brief Структурированная сетка (распределенная = false)
    bool structured() const { return m_structured; }

    /// @brief Число ячеек по оси x для структурированных сеток,
    /// число всех ячеек для сеток общего вида
    int nx() const { return m_nx; };

    /// @brief Число ячеек по оси y для структурированных сеток,
    /// единица для сеток общего вида
    int ny() const { return m_ny; };

    /// @brief Число ячеек по оси z для структурированных сеток,
    /// единица для сеток общего вида
    int nz() const { return m_nz; };

    /// @brief Получить ячейку по нескольким индексам подразумевая, что сетка
    /// является структурированной. Индексы периодически замкнуты (допускаются
    /// отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    EuCell operator()(index_t i, index_t j);

    /// @brief Получить ячейку по нескольким индексам подразумевая, что сетка
    /// является структурированной. Индексы периодически замкнуты (допускаются
    /// отрицательные индексы и индексы сверх нормы)
    /// @details Не актуально для распределенных сеток
    EuCell operator()(index_t i, index_t j, index_t k);

    /// @}

    /// @{ @name debug functions

    /// @brief Проверить базовую сетку
    int check_base() const;

    /// @brief Проверить сетку после адаптации
    int check_refined() const;

    /// @}

private:
    /// @brief Реальный конструктор сетки
    void build(geom::Generator& gen);

    /// @brief Добавить массив данных в хранилище
    /// @param name Имя массива данных
    template <typename T>
    Storable<T> add_one(const std::string& name);

    /// @brief Инициализация параметров AMR-ячеек
    void init_amr();

    /// @brief Балансировка флагов адаптации
    void balance_flags();

    /// @brief Применение флагов адаптации (после балансировки)
    void apply_flags();

    /// @brief Выставить новые ранги ячеек по декомпозиции
    void setup_ranks();

    /// @brief Собрать обменные слои (aliens и сопутствующие члены),
    /// принимаются актуальные данные с border слоя.
    void build_aliens();


    /// @brief Максимальный уровень адаптации для адаптивной сетки
    int m_max_level = 0;

    /// @brief Процедуры слияния и огрубления данных при адаптации
    Distributor m_distributor;

    AmrCells m_locals;  ///< Ячейки, которые принадлежат данному процессу
    AmrCells m_aliens;  ///< Ячейки, получаемые с других процессов

    /// @brief Метод декомпозиции
    Decomposition::Ptr m_decomp = nullptr;

#ifdef ZEPHYR_MPI
    Tourism   m_tourists;  ///< Построение обменных слоев и обмены
    Migration m_migrants;  ///< Пересылка ячеек при изменении декомпозиции
#endif

    /// @brief Структурированная сетка? (только для однопроцессорных)
    bool m_structured = false;
    int m_nx = 1, m_ny = 1, m_nz = 1;  ///< Размеры структурированной сетки
};


// ============================================================================
//                    Реализации шаблонных функций
// ============================================================================

template <typename T>
Storable<T> EuMesh::add_one(const std::string& name) {
    auto res1 = m_locals.data.add<T>(name);
    auto res2 = m_aliens.data.add<T>(name);
    if (res1 != res2) {
        throw std::runtime_error("EuMesh error: bad add<T> #1");
    }
#ifdef ZEPHYR_MPI
    auto res3 = m_tourists.add<T>(name);
    if (res1 != res3) {
        throw std::runtime_error("EuMesh error: bad add<T> #2");
    }
#endif
    return res1;
}

template <typename T>
void EuMesh::swap(Storable<T> var1, Storable<T> var2) {
    m_locals.data.swap<T>(var1, var2);
    m_aliens.data.swap<T>(var1, var2);
#ifdef ZEPHYR_MPI
    m_tourists.swap<T>(var1, var2);
#endif
}

template <typename... Args, typename >
void EuMesh::sync(Args&&... vars) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;
    m_tourists.sync(m_locals, m_aliens, std::forward<Args>(vars)...);
#endif
}

template <typename... Args>
void EuMesh::redistribute(Args&&... vars) {
#ifdef ZEPHYR_MPI
    if (mpi::single()) return;
    setup_ranks();
    m_migrants.migrate(
            m_tourists, m_locals, m_aliens,
            std::forward<Args>(vars)...);
    build_aliens();
#endif
}

} // namespace zephyr::mesh