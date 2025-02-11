#pragma once

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bnodes.h>

#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/eu_face.h>
#include <zephyr/mesh/euler/eu_cell.h>

#include <zephyr/mesh/euler/tourist_agency.h>
#include <zephyr/mesh/euler/migration_service.h>

#include <zephyr/mesh/decomp/ORB.h>

namespace zephyr::mesh {

using zephyr::utils::threads;
using zephyr::utils::mpi;
using namespace zephyr::geom;

using zephyr::mesh::decomp::Decomposition;
using zephyr::mesh::decomp::ORB;

/// @class Класс для хранения распределенной эйлеровой сетки.
/// Для сеток из четырехугольников и кубоидов допускается адаптация.
// TODO: Стркуктурировать поля и методы класса
class EuMesh {
public:
    /// @brief Конструктор сетки из генератора
    /// @tparam T Тип данных для хранения в ячейках
    /// @param gen Сеточный генератор
    template<class T>
    EuMesh(const T &val, Generator *gen)
            : m_locals(val, 0), m_aliens(val, 0) {
        if (mpi::master()) {
            initialize(gen->make());
        }
    }

    /// @brief Конструктор сетки из генератора
    /// @tparam T Тип данных для хранения в ячейках
    /// @param gen Сеточный генератор
    template<class T>
    EuMesh(const T &val, const std::shared_ptr<Generator>& gen)
            : m_locals(val, 0), m_aliens(val, 0) {
        if (mpi::master()) {
            initialize(gen->make());
        }
    }

    /// @brief Конструктор сетки из генератора
    /// @tparam T Тип данных для хранения в ячейках
    /// @param gen Сеточный генератор
    template<class T>
    EuMesh(const T &val, Generator &gen)
            : m_locals(val, 0), m_aliens(val, 0) {
        if (mpi::master()) {
            initialize(gen.make());
        }
    }

    /// @brief Можно сказать основной конструктор.
    /// Преобразует сетку Grid в специализированную EuMesh
    /// @tparam T Тип данных для хранения в ячейках
    /// @param grid Сетка общего вида
    template <class T>
    EuMesh(const T&val, const Grid& grid)
            : m_locals(val, 0), m_aliens(val, 0) {
        if (mpi::master()) {
            initialize(grid);
        }
    }

    /// @brief Локальное количество ячеек
    inline int n_cells() { return m_locals.size(); }

    /// @brief Итератор на начало локального массива ячеек
    EuCell begin() { return {m_locals, m_aliens, 0}; }

    /// @brief Итератор на конец локального массива ячеек
    EuCell end() { return {m_locals, m_aliens, m_locals.size()}; }

    /// @brief Доступ к локальной ячейке по индексу (idx < n_cells())
    EuCell operator[](int idx) {
        return {m_locals, m_aliens, idx};
    }

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

    /// @brief Неявное преобразование в хранилище. Позволяет передавать
    /// локальные ячейки в функции, которые работают с хранилищем.
    operator AmrStorage &() { return m_locals; }

    /// @brief Локальные ячейки
    AmrStorage &locals() { return m_locals; }

    /// @brief Ячейки с других процессов
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

    /// @brief Параллелизация по тредам
    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            class Value = typename std::result_of<Func(EuCell&)>::type>
    auto min(const Value &init, Func &&func)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            class Value = typename std::result_of<Func(EuCell&)>::type>
    auto min(Func &&func)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), std::forward<Func>(func));
    }

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
    void redistribute();

    /// @brief Обмен данными между процессами, в массивы aliens записываются
    /// данные с других процессов
    void exchange();
    void exchange_start();
    void exchange_end();

    /// @brief Балансирует нагрузку согласно decomposition
    void balancing(double load);

    /// @brief Дисбаланс нагрузки
    double get_imbalance(const std::vector<double>& ws) const {
        if(m_decomp)
            return m_decomp->imbalance(ws);
        return -1;
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

    /// @brief Число ячеек по оси x для структурированных сеток,
    /// число всех ячеек для сеток общего вида
    inline int nx() const { return m_nx; };

    /// @brief Число ячеек по оси y для структурированных сеток,
    /// единица для сеток общего вида
    inline int ny() const { return m_ny; };

    /// @brief Число ячеек по оси z для структурированных сеток,
    /// единица для сеток общего вида
    inline int nz() const { return m_nz; };

private:
    /// @brief Создать эйлерову сетку из сетки общего вида
    void initialize(const Grid& gen);

    /// @brief Очистить массив узлов, функция вызывается при любом
    /// перестроении сетки (вроде адаптации)
    void break_nodes();


public:
    // Сделать только балансировку флагов
    void balance_flags();

    // Применить флаги адаптации (после балансировки)
    void apply_flags();

    /// @brief Осуществляет инициализацию хранилища перед использованием
    /// функций адаптации, выполняется один раз после создания хранилища.
    void init_amr();

    /// @brief Перераспределить ячейки между процессами, в соответствии
    /// со значениями, которые выдает m_decomp
    void migrate();

public:
    /// @brief Собрать обменные слои (aliens и сопутствующие члены),
    /// не обязательно с пересылкой актуальных данных
    void build_aliens();


    int m_max_level = 0;
    Distributor distributor;

    AmrStorage m_locals;  ///< Локальные ячейки
    AmrStorage m_aliens;  ///< Ячейки с других процессов

    /// @brief Все данные, связанные с обменными слоями
    TouristAgency    m_tourism;

    /// @brief Метод декомпозиции
    Decomposition::Ptr m_decomp = nullptr;

    /// @brief Все данные, связанные с пересылкой ячеек
    MigrationService m_migration;

    /// @brief Массив уникальных узлов. Используется в редких алгоритмах,
    /// по умолчанию пустой, заполняются при вызове функции collect_nodes().
    std::vector<Vector3d> m_nodes;

    /// @brief Структура сетки, если предполагается, что сетка декартова.
    bool structured = false;
    int m_nx = 1, m_ny = 1, m_nz = 1;
};


} // namespace zephyr::mesh