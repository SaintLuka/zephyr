#pragma once

#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

class EuCell; // forward declaration

/// @defgroup euler-mesh Эйлерова сетка
/// @brief Обертки над сырыми массивами данных граней и ячеек.

class EuFace_Iter;

/// @brief Грань эйлеровой ячейки
/// @ingroup euler-mesh
class EuFace final {
    friend class EuFace_Iter;

    using Vector3d = geom::Vector3d;
    using Boundary = geom::Boundary;

private:
    AmrCells* m_cells;     //< Указатель на сетку (обычно locals)
    index_t   m_face_idx;  //< Индекс первой грани

    /// @brief Нулевое значение допускается, но не проверяется в целях
    /// оптимизации, поэтому ловите segfaults.
    AmrCells* m_aliens = nullptr;

public:
    /// @brief Основной конструктор
    EuFace(AmrCells* cells, index_t face_idx, AmrCells* aliens = nullptr)
        : m_cells(cells), m_face_idx(face_idx), m_aliens(aliens) { }

    /// @{ @name Тип грани

    /// @brief Флаг граничных условий
    Boundary flag() const;

    /// @brief Является ли грань граничной?
    bool is_boundary() const;

    /// @brief Является ли грань актуальной?
    bool is_actual() const;

    /// @return 'true', если грань не актуальна
    bool is_undefined() const;

    /// @brief Установить неопределенную грань
    void set_undefined() const;

    /// @}

    /// @{ @name Геометрия грани

    /// @brief Внешняя нормаль
    const Vector3d& normal() const;

    /// @brief Барицентр грани
    const Vector3d& center() const;

    /// @brief Координата x центра грани
    double x() const { return center().x(); }

    /// @brief Координата y центра грани
    double y() const { return center().y(); }

    /// @brief Координата z центра грани
    double z() const { return center().z(); }

    /// @brief Сторона, по которой расположена грань
    template <int dim = 3>
    Side<dim> side() const;

    /// @brief Площадь/длина обычной грани
    double area() const;

    /// @brief Площадь/длина на внешнюю нормаль
    Vector3d area_n() const;

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double area(bool axial) const;

    /// @brief Площадь грани осесимметричной ячейки
    double area_as() const;

    /// @brief Число вершин грани
    int n_vertices() const;

    /// @brief Локальный индекс вершины (в ячейке)
    index_t vertex_index(int idx) const;

    /// @brief Глобальный индекс вершины
    index_t node_index(int idx) const;

    /// @brief Получить вершину грани
    Vector3d vs(int idx) const;

    /// @brief Точка симметричная относительно грани
    Vector3d symm_point(const Vector3d& p) const;

    /// @}

    /// @{ @name Свойства соседней ячейки

    /// @brief Ранг соседней ячейки
    int adj_rank() const;

    /// @brief Индекс соседней ячейки в массиве locals
    index_t adj_index() const;

    /// @brief Индекс соседней ячейки в массиве aliens
    index_t adj_alien() const;

    /// @brief Индекс родительской ячейки в массиве locals
    index_t adj_basic() const;

    /// @brief Проверяет, является сосед через грань локальным
    bool local_neib() const;

    /// @brief Соседняя ячейка по внешней нормали, на границе сетки
    /// гарантированно возвращается сама ячейка
    EuCell neib() const;

    /// @brief Получить ссылку на данные соседа
    template <typename T>
    const T& neib(Storable<T> type) const;

    /// @brief Флаг адаптации соседней ячейки
    int neib_flag() const;

    /// @brief Центр соседней ячейки
    Vector3d neib_center() const;

    /// @brief Объем соседней ячейки
    double neib_volume() const;

    /// @brief Объем соседней ячейки
    double neib_volume(bool axial) const;

    /// @}
};

/// @brief Итератор по граням ячейки
class EuFace_Iter final {
    EuFace    m_eu_face;   //< Действительная грань
    index_t   m_face_end;  //< Индекс за последней гранью
    Direction m_dir;       //< Выбранное направление граней

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    EuFace_Iter(AmrCells* cells, index_t face_idx, index_t face_end,
                AmrCells* aliens, Direction dir = Direction::ANY);

    /// @brief Ссылка на грань при разыменовании
    EuFace &operator*() { return m_eu_face; }

    /// @brief Ссылка на грань при разыменовании
    const EuFace &operator*() const { return m_eu_face; }

    /// @brief Перейти к следующей определенной грани
    EuFace_Iter &operator++();

    /// @brief Сравнение итераторов
    bool operator!=(const EuFace_Iter &face) const;

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределенна или не соответствует направлению
    bool to_skip(Direction dir) const;
};


/// @brief Интерфейс для итераций по граням ячейки
class EuFaces final {
    EuFace_Iter m_begin;
    EuFace_Iter m_end;

public:
    EuFaces(AmrCells *cells, index_t cell_idx,
            AmrCells *aliens = nullptr,
            Direction dir = Direction::ANY);

    EuFace_Iter begin() const { return m_begin; }

    EuFace_Iter end() const { return m_end; }
};


class EuCell_Iter;

/// @brief Эйлерова ячейка
/// @ingroup euler-mesh
class EuCell final {
    friend class EuCell_Iter;

    using Vector3d = geom::Vector3d;

    /// @brief Характеристическая функция (функция-индикатор)
    using InFunction = std::function<bool(const Vector3d &)>;

    /// @brief Пространственная функция
    using SpFunction = std::function<double(const Vector3d &)>;

private:
    AmrCells* m_cells;  //< Указатель на сетку (обычно locals)
    index_t   m_index;  //< Индекс ячейки

    /// @brief Нулевое значение допускается, но не проверяется в целях
    /// оптимизации, поэтому ловите segfaults.
    AmrCells* m_aliens = nullptr;

public:
    EuCell(AmrCells* cells, index_t index, AmrCells* aliens = nullptr)
        : m_cells(cells), m_index(index), m_aliens(aliens) { }

    /// @{ @name Характеристики ячейки

    /// @brief Размерность ячейки
    int dim() const;

    /// @brief Адаптивная (AMR) ячейка?
    bool adaptive() const;

    /// @brief Ранг, которому принадлежит ячейка
    int rank() const;

    /// @brief Индекс базовой родительской ячейки
    index_t b_idx() const;

    /// @brief Индекс ячейки на z-кривой
    index_t z_idx() const;

    /// @brief Текущий индекс ячейки в массиве
    index_t index() const;

    /// @brief Индекс новой ячейки (в алгоритмах)
    index_t next() const;

    /// @brief Флаг адаптации ячейки
    index_t flag() const;

    /// @brief Уровень адаптации ячейки
    index_t level() const;

    /// @brief Установить ранг ячейки
    void set_rank(int rank) const;

    /// @brief Установить флаг адаптации
    ///   flag = -1: огрубление/слияние;
    ///   flag =  0: ничего не делать;
    ///   flag = +1: разбиение ячейки.
    void set_flag(int flag) const;

    /// @}

    /// @{ @name Геометрия ячейки

    /// @brief Центр ячейки
    const Vector3d& center() const;
    double x() const { return m_cells->center[m_index].x(); }
    double y() const { return m_cells->center[m_index].y(); }
    double z() const { return m_cells->center[m_index].z(); }

    /// @brief Объем ячейки (площадь в двумерном случае)
    double volume() const;

    /// @brief Объем ячейки
    double volume(bool axial) const;

    /// @brief Объем осесимметричной ячейки
    double volume_as() const;

    /// @brief Линейный размер ячейки по оси x (от левой до правой грани)
    double hx() const;

    /// @brief Линейный размер ячейки по оси y (от нижней до верхней грани)
    double hy() const;

    /// @brief Линейный размер ячейки по оси z (от задней до передней грани)
    double hz() const;

    /// @brief Линейный размер ячейки
    double linear_size() const;

    /// @brief Диаметр вписанной окружности
    double incircle_diameter() const;

    /// @brief Bounding box ячейки
    geom::Box bbox() const;

    /// @brief Создать полигон из ячейки (для 3D ячеек -- UB)
    geom::Polygon polygon() const;

    /// @brief Создать полигон из ячейки (для 2D ячеек -- UB)
    geom::Polyhedron polyhedron() const;

    /// @}

    /// @{ @name Данные ячейки

    /// @brief Ссылка на данные ячейки
    template <typename T>
    T& operator()(Storable<T> var);
    template <typename T>
    T& operator[](Storable<T> var);

    /// @brief Константная ссылка на данные ячейки
    template <typename T>
    const T& operator()(Storable<T> var) const;
    template <typename T>
    const T& operator[](Storable<T> var) const;

    /// @brief Ссылка на данные ячейки, вызов актуален, только если сетка
    /// содержит единственный массив данных такого типа
    template <typename T>
    T& operator()(const T& var);
    template <typename T>
    T& operator[](const T& var);

    /// @brief Константная ссылка на данные ячейки, вызов актуален, только
    /// если сетка содержит единственный массив данных такого типа
    template <typename T>
    const T& operator()(const T& var) const;
    template <typename T>
    const T& operator[](const T& var) const;

    /// @brief Скопировать данные в другую ячейку
    void copy_data_to(EuCell& dst_cell) const;

    /// @}

    /// @{ @name Грани ячейки

    /// @brief Число актуальных граней ячейки
    int face_count() const;

    /// @brief Получить грань по индексу в ячейке
    EuFace face(int idx) const;

    /// @brief Получить грань двумерной ячейки
    EuFace face(Side2D s) const;

    /// @brief Получить грань трёхмерной ячейки
    EuFace face(Side3D s) const;

    /// @brief Простая грань на выбранной стороне?
    bool simple_face(Side2D s) const;

    /// @brief Простая грань на выбранной стороне?
    bool simple_face(Side3D s) const;

    /// @brief Сложная грань на выбранной стороне?
    bool complex_face(Side2D s) const;

    /// @brief Сложная грань на выбранной стороне?
    bool complex_face(Side3D s) const;

    /// @brief Итератор по граням ячейки
    EuFaces faces(Direction dir = Direction::ANY) const;

    /// @}

    /// @{ @name Вершины ячейки

    /// @brief Число вершин, для AMR-ячеек полное число вершин (9 или 27).
    int node_count() const;

    /// @brief Указатель на первую вершину
    const Vector3d* vertices_data() const;

    /// @brief Вершины как набор узлов квадратичного отображения
    template <int dim>
    const SqMap<dim>& mapping() const { return m_cells->mapping<dim>(m_index); }

    /// @}

    /// @{ @name Соседние ячейки

    /// @brief Помещает на место текущей ячейки соседнюю ячейку, которая
    /// находится со стороны loc_face.
    void replace(int loc_face);

    /// @brief Получить соседнюю ячейку на двумерной сетке, заданы смещения
    /// относительно ячейки по осям. Функция работает, если вокруг ячейки можно
    /// построить структурированный шаблон, который также целиком расположен
    /// на одном процессе. В остальных случаях неопределенное поведение.
    /// @code
    ///     for (auto cell: mesh) {
    ///         auto neib_L = cell.neib(-1, 0);  // сосед слева
    ///         auto neib_R = cell.neib(+1, 0);  // сосед справа
    ///     }
    /// @endcode
    EuCell neib(index_t i, index_t j) const;

    /// @brief Получить соседнюю ячейку на трёхмерной сетке, заданы смещения
    /// относительно ячейки по осям. Функция работает, если вокруг ячейки можно
    /// построить структурированный шаблон, который также целиком расположен
    /// на одном процессе. В остальных случаях неопределенное поведение.
    /// @code
    ///     for (auto cell: mesh) {
    ///         auto neib_L = cell.neib(-1, 0, 0);  // сосед слева
    ///         auto neib_R = cell.neib(+1, 0, 0);  // сосед справа
    ///     }
    /// @endcode
    EuCell neib(index_t i, index_t j, index_t k) const;

    /// @}

    /// @{ @name Сечения и интегрирование по ячейке

    /// @brief Оценка объемной доли, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @details Относительно быстрая функция, проверяет функцию inside только
    /// на узлах ячейки, позволяет быстро выяснить, содержит ли ячейка
    /// границу двух областей. Если ячейка внутри тела, то возвращает строго
    /// единицу 1.0, если снаружи -- строго ноль 0.0.
    double approx_vol_fraction(const SpFunction& inside) const;

    /// @brief Объемная доля, которая отсекается от ячейки некоторым телом.
    /// @param inside Характеристическая функция области, возвращает true для
    /// точек, которые располагаются внутри области.
    /// @param n_points Число тестовых точек, для которых проверяется функция
    /// inside, погрешность определения объемной доли ~ 1/N.
    double volume_fraction(const SpFunction& inside, int n_points) const;

    /// @brief Функция func является константой на ячейке?
    /// @details Проверяется значение функции в узлах и в центре ячейки,
    /// если все значения совпадают, то считается, что функция принимает
    /// постоянное значение в пределах ячейки.
    bool const_function(const SpFunction& func) const;

    /// @brief Интеграл скалярной функции по ячейке
    /// @param n_points Разбиение по сторонам
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(const SpFunction& func, int n_points) const;

    /// @}
};

/// @brief Итератор по ячейкам из EuMesh или AmrCells
class EuCell_Iter final {
private:
    EuCell m_eu_cell;  ///< Реальная ячейка

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = EuCell;
    using pointer    = EuCell *;
    using reference  = EuCell &;

    /// @brief Конструктор как у ячейки
    EuCell_Iter(AmrCells *cells, index_t index, AmrCells *aliens = nullptr)
            : m_eu_cell{cells, index, aliens} { }

    /// @brief Ссылка на ячейку при разыменовании
    EuCell &operator*() { return m_eu_cell; }

    /// @brief Ссылка на ячейку при разыменовании
    const EuCell &operator*() const { return m_eu_cell; }

    /// @brief Инкремент
    EuCell_Iter &operator++() {
        ++m_eu_cell.m_index;
        return *this;
    }

    /// @brief Декремент
    EuCell_Iter &operator--() {
        --m_eu_cell.m_index;
        return *this;
    }

    /// @brief Итератор через step
    EuCell_Iter &operator+=(index_t step) {
        m_eu_cell.m_index += step;
        return *this;
    }

    /// @brief Итератор через step
    EuCell_Iter operator+(index_t step) const {
        return {m_eu_cell.m_cells,
                m_eu_cell.m_index + step,
                m_eu_cell.m_aliens};
    }

    /// @brief Оператор доступа как для указателя (random access iterator)
    EuCell operator[](index_t offset) const {
        return {m_eu_cell.m_cells,
                m_eu_cell.m_index + offset,
                m_eu_cell.m_aliens};
    }

    /// @brief Расстояние между двумя ячейками
    index_t operator-(const EuCell_Iter &cell) const {
        return m_eu_cell.m_index - cell.m_eu_cell.m_index;
    }

    /// @brief Оператор сравнения
    bool operator<(const EuCell_Iter &cell) const {
        return m_eu_cell.m_index < cell.m_eu_cell.m_index;
    }

    /// @brief Оператор сравнения
    bool operator!=(const EuCell_Iter &cell) const {
        return m_eu_cell.m_index != cell.m_eu_cell.m_index;
    }

    /// @brief Оператор сравнения
    bool operator==(const EuCell_Iter &cell) const {
        return m_eu_cell.m_index == cell.m_eu_cell.m_index;
    }
};

/// @brief Для итераций по AmrCells, aliens = nullptr,
/// поэтому проход по соседям не всегда возможен
inline EuCell_Iter begin(AmrCells& cells) {
    return {&cells, 0, nullptr};
}

/// @brief Для итераций по AmrCells, aliens = nullptr,
/// поэтому проход по соседям не всегда возможен
inline EuCell_Iter end(AmrCells& cells) {
    return {&cells, cells.n_cells(), nullptr};
}

/// @brief Набор дочерних ячеек (сиблингов).
/// @ingroup euler-mesh
///
/// Предполагается, что дочерние ячейки располагаются в локальном хранилище.
/// Во время операций split (refine) и merge (coarse) сетка может находиться
/// в не совместном состоянии, поэтому переход по соседям запрещен.
class Children final {
public:
    /// @brief Индексы дочерних ячеек
    std::array<index_t, 8> index = {-1, -1, -1, -1, -1, -1, -1, -1};

    /// @brief Обязательная инициализация локального хранилища
    explicit Children(AmrCells* locals) : m_locals(locals) { }

    /// @brief Размерность определяется по числу дочерних ячеек
    int dim() const { return index[4] < 0 ? 2 : 3; }

    /// @brief Число дочерних ячеек (4 для 2D и 8 для 3D)
    int count() const { return index[4] < 0 ? 4 : 8; }

    /// @brief Получить дочернюю ячейку по индексу
    EuCell operator[](int idx) const { return {m_locals, index[idx], nullptr}; };

    /// @brief Итератор по дочерним ячейкам
    struct iterator {
        iterator(const Children& children, int idx)
            : m_children(children), m_idx(idx) { }

        EuCell operator*() const { return m_children[m_idx]; }

        void operator++() { ++m_idx; }

        bool operator!=(const iterator &it) const {
            return m_idx != it.m_idx;
        }
    private:
        const Children& m_children;
        int m_idx;
    };

    /// @brief Первая дочерняя ячейка
    iterator begin() const { return {*this, 0}; }

    /// @brief За последней дочерней ячейкой
    iterator end() const { return {*this, count()}; }

private:
    AmrCells* m_locals;  ///< Локальное хранилище с ячейками
};

// ================================================================================================
//                                     inline функции EuFace
// ================================================================================================

inline geom::Boundary EuFace::flag() const { return m_cells->faces.boundary[m_face_idx]; }

inline bool EuFace::is_boundary() const { return m_cells->faces.is_boundary(m_face_idx); }

inline bool EuFace::is_actual() const { return m_cells->faces.is_actual(m_face_idx); }

inline bool EuFace::is_undefined() const { return m_cells->faces.is_undefined(m_face_idx); }

inline void EuFace::set_undefined() const { m_cells->faces.set_undefined(m_face_idx); }

inline const geom::Vector3d &EuFace::normal() const { return m_cells->faces.normal[m_face_idx]; }

inline const geom::Vector3d &EuFace::center() const { return m_cells->faces.center[m_face_idx]; }

template <int dim >
Side<dim> EuFace::side() const {
    index_t cell_idx = m_cells->faces.adjacent.basic[m_face_idx];
    return m_face_idx - m_cells->face_begin[cell_idx];
}

inline double EuFace::area() const { return m_cells->faces.area[m_face_idx]; }

inline geom::Vector3d EuFace::area_n() const { return m_cells->faces.area[m_face_idx] * m_cells->faces.normal[m_face_idx]; }

inline double EuFace::area(bool axial) const { return m_cells->faces.get_area(m_face_idx, axial); }

inline double EuFace::area_as() const { return m_cells->faces.area_alt[m_face_idx]; }

inline int EuFace::n_vertices() const { return m_cells->faces.n_vertices(m_face_idx); }

inline index_t EuFace::vertex_index(int idx) const { return m_cells->faces.vertices[m_face_idx][idx]; }

inline index_t EuFace::node_index(int idx) const {
    index_t cell_idx = m_cells->faces.adjacent.basic[m_face_idx];
    return m_cells->node_begin[cell_idx] + static_cast<int>(m_cells->faces.vertices[m_face_idx][idx]);
}

inline geom::Vector3d EuFace::vs(int idx) const { return m_cells->verts[node_index(idx)]; }

inline geom::Vector3d EuFace::symm_point(const Vector3d &p) const {
    return m_cells->faces.symm_point(m_face_idx, p);
}

inline int EuFace::adj_rank() const { return m_cells->faces.adjacent.rank[m_face_idx]; }

inline index_t EuFace::adj_index() const { return m_cells->faces.adjacent.index[m_face_idx]; }

inline index_t EuFace::adj_alien() const { return m_cells->faces.adjacent.alien[m_face_idx]; }

inline index_t EuFace::adj_basic() const { return m_cells->faces.adjacent.basic[m_face_idx]; }

inline bool EuFace::local_neib() const {
    return m_cells->faces.adjacent.is_local(m_face_idx);
}

inline EuCell EuFace::neib() const {
    if (utils::mpi::single() || local_neib()) {
        return {m_cells, adj_index(), m_aliens};
    }
    return {m_aliens, adj_alien(), m_aliens};
}

template <typename T>
const T& EuFace::neib(Storable<T> type) const {
    if (utils::mpi::single() || local_neib()) {
        return m_cells->data.get_val(type, adj_index());
    }
    return m_aliens->data.get_val(type, adj_alien());
}

inline int EuFace::neib_flag() const {
    if (utils::mpi::single() || local_neib()) {
        return m_cells->flag[adj_index()];
    }
    return m_aliens->flag[adj_alien()];
}

inline geom::Vector3d EuFace::neib_center() const {
    if (utils::mpi::single() || local_neib()) {
        return m_cells->center[adj_index()];
    }
    return m_aliens->center[adj_alien()];
}

inline double EuFace::neib_volume() const {
    if (utils::mpi::single() || local_neib()) {
        return m_cells->volume[adj_index()];
    }
    return m_aliens->volume[adj_alien()];
}

inline double EuFace::neib_volume(bool axial) const {
    if (axial) {
        if (utils::mpi::single() || local_neib()) {
            return m_cells->volume_alt[adj_index()];
        }
        return m_aliens->volume_alt[adj_alien()];
    }
    return neib_volume();
}

// ================================================================================================
//                                     inline функции EuCell
// ================================================================================================

inline int EuCell::dim() const { return m_cells->dim(); }

inline bool EuCell::adaptive() const { return m_cells->adaptive(); }

inline int EuCell::rank() const { return m_cells->rank[m_index]; }

inline int EuCell::flag() const { return m_cells->flag[m_index]; }

inline int EuCell::level() const { return m_cells->level[m_index]; }

inline index_t EuCell::b_idx() const { return m_cells->b_idx[m_index]; }

inline index_t EuCell::z_idx() const { return m_cells->z_idx[m_index]; }

inline index_t EuCell::index() const { return m_cells->index[m_index]; }

inline index_t EuCell::next() const { return m_cells->next[m_index]; }

inline void EuCell::set_rank(int rank) const { m_cells->rank[m_index] = rank; }

inline void EuCell::set_flag(int flag) const { m_cells->flag[m_index] = flag; }

inline const geom::Vector3d& EuCell::center() const { return m_cells->center[m_index]; }

inline double EuCell::volume() const { return m_cells->volume[m_index]; }

inline double EuCell::volume(bool axial) const { return m_cells->get_volume(m_index, axial); }

inline double EuCell::volume_as() const { return m_cells->volume_alt[m_index]; }

inline double EuCell::hx() const { return m_cells->hx(m_index); }

inline double EuCell::hy() const { return m_cells->hy(m_index); }

inline double EuCell::hz() const { return m_cells->hz(m_index); }

inline double EuCell::linear_size() const { return m_cells->linear_size(m_index); }

inline double EuCell::incircle_diameter() const { return m_cells->incircle_diameter(m_index); }

template <typename T>
T& EuCell::operator()(Storable<T> var) { return m_cells->data.get_val(var, m_index); }
template <typename T>
T& EuCell::operator[](Storable<T> var) { return m_cells->data.get_val(var, m_index); }

template <typename T>
const T& EuCell::operator()(Storable<T> var) const { return m_cells->data.get_val(var, m_index); }
template <typename T>
 const T& EuCell::operator[](Storable<T> var) const { return m_cells->data.get_val(var, m_index); }

template <typename T>
T& EuCell::operator()(const T& var) { return m_cells->data.get_val<T>(m_index); }
template <typename T>
T& EuCell::operator[](const T& var) { return m_cells->data.get_val<T>(m_index); }

template <typename T>
const T& EuCell::operator()(const T& var) const { return m_cells->data.get_val<T>(m_index); }
template <typename T>
const T& EuCell::operator[](const T& var) const { return m_cells->data.get_val<T>(m_index); }

inline void EuCell::copy_data_to(EuCell &dst_cell) const {
    m_cells->copy_data(m_index, dst_cell.m_cells, dst_cell.m_index);
}

inline int EuCell::face_count() const { return m_cells->face_count(m_index); }

inline EuFace EuCell::face(int idx) const {
    return {m_cells, m_cells->face_begin[m_index] + idx, m_aliens};
}

inline EuFace EuCell::face(Side2D s) const {
    return {m_cells, m_cells->face_begin[m_index] + s, m_aliens};
}

inline EuFace EuCell::face(Side3D s) const {
    return {m_cells, m_cells->face_begin[m_index] + s, m_aliens};
}

inline bool EuCell::simple_face(Side2D s) const { return m_cells->simple_face(m_index, s); }

inline bool EuCell::simple_face(Side3D s) const { return m_cells->simple_face(m_index, s); }

inline bool EuCell::complex_face(Side2D s) const { return m_cells->complex_face(m_index, s); }

inline bool EuCell::complex_face(Side3D s) const { return m_cells->complex_face(m_index, s); }

inline EuFaces EuCell::faces(Direction dir) const { return {m_cells, m_index, m_aliens, dir}; }

inline int EuCell::node_count() const { return m_cells->node_count(m_index); }

inline const geom::Vector3d* EuCell::vertices_data() const { return m_cells->vertices_data(m_index); }

inline double EuCell::approx_vol_fraction(const SpFunction& inside) const {
    return m_cells->approx_vol_fraction(m_index, inside);
}

inline double EuCell::volume_fraction(const SpFunction& inside, int n_points) const {
    return m_cells->volume_fraction(m_index, inside, n_points);
}

inline bool EuCell::const_function(const SpFunction& func) const {
    return m_cells->const_function(m_index, func);
}

inline double EuCell::integrate_low(const SpFunction& func, int n_points) const {
    return m_cells->integrate_low(m_index, func, n_points);
}

} // namespace zephyr::mesh