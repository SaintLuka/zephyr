#pragma once

#include <zephyr/mesh/raw/raw_cells.h>
#include <zephyr/mesh/raw/raw_nodes.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

class Cell; // forward declaration

/// @defgroup raw-mesh Эйлерова сетка
/// @brief Обертки над сырыми массивами данных граней и ячеек.

class Face_Iter;

/// @brief Грань эйлеровой ячейки
/// @ingroup raw-mesh
class Face final {
    friend class Face_Iter;

    using Vector3d = geom::Vector3d;
    using Boundary = geom::Boundary;

private:
    RawCells* cells_;     //< Указатель на сетку (обычно locals)
    index_t   face_idx_;  //< Индекс первой грани

    /// @brief Нулевое значение допускается, но не проверяется в целях
    /// оптимизации, поэтому ловите segfaults.
    RawCells* ghosts_ = nullptr;

public:
    /// @brief Основной конструктор
    Face(RawCells* cells, index_t face_idx, RawCells* ghosts = nullptr)
        : cells_(cells), face_idx_(face_idx), ghosts_(ghosts) { }

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

    /// @brief Установить флаг граничных условий
    void set_boundary(Boundary flag) const;

    /// @brief Поворот соседней ячейки
    int rotation() const;

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
    geom::Side<dim> side() const;

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

    /// @brief Индекс соседней ячейки в массиве ghosts
    index_t adj_ghost() const;

    /// @brief Индекс родительской ячейки в массиве locals
    index_t adj_basic() const;

    /// @brief Проверяет, является сосед через грань локальным
    bool local_neib() const;

    /// @brief Соседняя ячейка по внешней нормали, на границе сетки
    /// гарантированно возвращается сама ячейка
    Cell neib() const;

    /// @brief Получить ссылку на данные соседа
    template <typename T>
    const T& neib(Storable<T> type) const;

    template <typename T>
    std::span<const T> neib(Storable<T[]> type) const;

    template <typename T, size_t N>
    std::span<const T, N> neib(Storable<T[N]> type) const;

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
class Face_Iter final {
    Face    face_;      //< Текущая грань
    index_t   face_end_;  //< Индекс за последней гранью
    Direction dir_;       //< Выбранное направление граней

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    Face_Iter(RawCells* cells, index_t face_idx, index_t face_end,
                RawCells* ghosts, Direction dir = Direction::ANY);

    /// @brief Ссылка на грань при разыменовании
    Face &operator*() { return face_; }

    /// @brief Ссылка на грань при разыменовании
    const Face &operator*() const { return face_; }

    /// @brief Перейти к следующей определенной грани
    Face_Iter &operator++();

    /// @brief Сравнение итераторов
    bool operator!=(const Face_Iter &face) const;

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределенна или не соответствует направлению
    bool to_skip(Direction dir) const;
};


/// @brief Интерфейс для итераций по граням ячейки
class FacesRange final {
    Face_Iter begin_;
    Face_Iter end_;

public:
    FacesRange(RawCells *cells, index_t cell_idx,
            RawCells *ghosts = nullptr,
            Direction dir = Direction::ANY);

    Face_Iter begin() const { return begin_; }

    Face_Iter end() const { return end_; }
};


class Cell_Iter;

/// @brief Эйлерова ячейка
/// @ingroup raw-mesh
class Cell final {
    friend class Cell_Iter;

    using Vector3d = geom::Vector3d;

    /// @brief Характеристическая функция (функция-индикатор)
    using InFunction = std::function<bool(const Vector3d &)>;

    /// @brief Пространственная функция
    using SpFunction = std::function<double(const Vector3d &)>;

private:
    RawCells* cells_{nullptr};  //< Указатель на сетку (обычно locals)
    index_t   index_{-1};       //< Индекс ячейки

    /// @brief Нулевое значение допускается, но не проверяется в целях
    /// оптимизации, поэтому ловите segfaults.
    RawCells* ghosts_{nullptr};

public:
    /// @brief Конструктор по умолчанию (иногда требует TBB)
    Cell() = default;

    Cell(RawCells* cells, index_t index, RawCells* ghosts = nullptr)
        : cells_(cells), index_(index), ghosts_(ghosts) { }

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

    /// @brief Реальное положение в массиве
    index_t id() const;

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

    /// @brief Хранилище, которому принадлежит ячейка
    const RawCells& cells() const { return *cells_; }

    /// @}

    /// @{ @name Геометрия ячейки

    /// @brief Центр ячейки
    const Vector3d& center() const;
    double x() const { return cells_->center[index_].x(); }
    double y() const { return cells_->center[index_].y(); }
    double z() const { return cells_->center[index_].z(); }

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
    T& operator[](Storable<T> var);

    /// @brief Ссылка на данные ячейки
    template <typename T>
    std::span<T> operator[](Storable<T[]> var);

    /// @brief Ссылка на данные ячейки
    template <typename T, size_t N>
    std::span<T, N> operator[](Storable<T[N]> var);

    /// @brief Константная ссылка на данные ячейки
    template <typename T>
    const T& operator[](Storable<T> var) const;

    /// @brief Константная ссылка на данные ячейки
    template <typename T>
    std::span<const T> operator[](Storable<T[]> var) const;

    /// @brief Константная ссылка на данные ячейки
    template <typename T, size_t N>
    std::span<const T, N> operator[](Storable<T[N]> var) const;

    /// @brief Скопировать данные в другую ячейку
    void copy_data_to(Cell& dst_cell) const;

    /// @}

    /// @{ @name Грани ячейки

    /// @brief Число актуальных граней ячейки
    int face_count() const;

    /// @brief Получить грань по индексу в ячейке
    Face face(int idx) const;

    /// @brief Получить грань двумерной ячейки
    Face face(Side2D s) const;

    /// @brief Получить грань трёхмерной ячейки
    Face face(Side3D s) const;

    /// @brief Простая грань на выбранной стороне?
    bool simple_face(Side2D s) const;

    /// @brief Простая грань на выбранной стороне?
    bool simple_face(Side3D s) const;

    /// @brief Сложная грань на выбранной стороне?
    bool complex_face(Side2D s) const;

    /// @brief Сложная грань на выбранной стороне?
    bool complex_face(Side3D s) const;

    /// @brief Итератор по граням ячейки
    FacesRange faces(Direction dir = Direction::ANY) const;

    /// @}

    /// @{ @name Вершины ячейки

    /// @brief Число вершин, для AMR-ячеек полное число вершин (9 или 27).
    int node_count() const;

    /// @brief Указатель на первую вершину
    const Vector3d* vertices_data() const;

    /// @brief Вершины как набор узлов квадратичного отображения
    template <int dim>
    const SqMap<dim>& mapping() const { return cells_->verts.mapping<dim>(index_); }

    int node_rank(int iv) const;
    int node_index(int iv) const;
    int node_ghost(int iv) const;

    /// @}

    /// @{ @name Соседние ячейки

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
    Cell neib(index_t i, index_t j) const;

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
    Cell neib(index_t i, index_t j, index_t k) const;

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

/// @brief Итератор по ячейкам из Mesh или RawCells
class Cell_Iter final {
private:
    Cell cell_{};  ///< Реальная ячейка

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = Cell;
    using pointer    = Cell *;
    using reference  = Cell &;

    /// @brief Конструктор по умолчанию (иногда требует TBB)
    Cell_Iter() = default;

    /// @brief Конструктор как у ячейки
    Cell_Iter(RawCells *cells, index_t index, RawCells *ghosts = nullptr)
            : cell_{cells, index, ghosts} { }

    /// @brief Ссылка на ячейку при разыменовании
    Cell &operator*() { return cell_; }

    /// @brief Ссылка на ячейку при разыменовании
    const Cell &operator*() const { return cell_; }

    /// @brief Инкремент
    Cell_Iter &operator++() {
        ++cell_.index_;
        return *this;
    }

    /// @brief Декремент
    Cell_Iter &operator--() {
        --cell_.index_;
        return *this;
    }

    /// @brief Итератор через step
    Cell_Iter &operator+=(index_t step) {
        cell_.index_ += step;
        return *this;
    }

    /// @brief Итератор через step
    Cell_Iter operator+(index_t step) const {
        return {cell_.cells_,
                cell_.index_ + step,
                cell_.ghosts_};
    }

    /// @brief Оператор доступа как для указателя (random access iterator)
    Cell operator[](index_t offset) const {
        return {cell_.cells_,
                cell_.index_ + offset,
                cell_.ghosts_};
    }

    /// @brief Расстояние между двумя ячейками
    index_t operator-(const Cell_Iter &cell) const {
        return cell_.index_ - cell.cell_.index_;
    }

    /// @brief Оператор сравнения
    bool operator<(const Cell_Iter &cell) const {
        return cell_.index_ < cell.cell_.index_;
    }

    /// @brief Оператор сравнения
    bool operator!=(const Cell_Iter &cell) const {
        return cell_.index_ != cell.cell_.index_;
    }

    /// @brief Оператор сравнения
    bool operator==(const Cell_Iter &cell) const {
        return cell_.index_ == cell.cell_.index_;
    }
};

/// @brief Для итераций по RawCells, ghosts = nullptr,
/// поэтому проход по соседям не всегда возможен
inline Cell_Iter begin(RawCells& cells) {
    return {&cells, 0, nullptr};
}

/// @brief Для итераций по RawCells, ghosts = nullptr,
/// поэтому проход по соседям не всегда возможен
inline Cell_Iter end(RawCells& cells) {
    return {&cells, cells.n_cells(), nullptr};
}

/// @brief Набор дочерних ячеек (сиблингов).
/// @ingroup raw-mesh
///
/// Предполагается, что дочерние ячейки располагаются в локальном хранилище.
/// Во время операций split (refine) и merge (coarse) сетка может находиться
/// в не совместном состоянии, поэтому переход по соседям запрещен.
class Children final {
public:
    /// @brief Индексы дочерних ячеек
    std::array<index_t, 8> index = {-1, -1, -1, -1, -1, -1, -1, -1};

    /// @brief Обязательная инициализация локального хранилища
    explicit Children(RawCells* locals) : locals_(locals) { }

    /// @brief Размерность определяется по числу дочерних ячеек
    int dim() const { return index[4] < 0 ? 2 : 3; }

    /// @brief Число дочерних ячеек (4 для 2D и 8 для 3D)
    int count() const { return index[4] < 0 ? 4 : 8; }

    /// @brief Получить дочернюю ячейку по индексу
    Cell operator[](int idx) const { return {locals_, index[idx], nullptr}; };

    /// @brief Итератор по дочерним ячейкам
    struct iterator {
        iterator(const Children& children, int idx)
            : m_children(children), m_idx(idx) { }

        Cell operator*() const { return m_children[m_idx]; }

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
    RawCells* locals_;  ///< Локальное хранилище с ячейками
};

// ================================================================================================
//                                     inline функции Face
// ================================================================================================

inline geom::Boundary Face::flag() const { return cells_->faces.boundary[face_idx_]; }

inline bool Face::is_boundary() const { return cells_->faces.is_boundary(face_idx_); }

inline bool Face::is_actual() const { return cells_->faces.is_actual(face_idx_); }

inline bool Face::is_undefined() const { return cells_->faces.is_undefined(face_idx_); }

inline void Face::set_undefined() const { cells_->faces.set_undefined(face_idx_); }

inline void Face::set_boundary(Boundary flag) const {
    if (flag == Boundary::INNER || flag == Boundary::PERIODIC) {
        std::cerr << "You can't just set Ordinary or Periodic boundary flag\n";
    }
    else {
        index_t basic = cells_->faces.adjacent.basic[face_idx_];
        cells_->faces.boundary[face_idx_] = flag;
        cells_->faces.adjacent.rank[face_idx_] = cells_->rank[basic];
        cells_->faces.adjacent.index[face_idx_] = basic;
        cells_->faces.adjacent.ghost[face_idx_] = -1;
    }
}

inline int Face::rotation() const { return cells_->faces.adjacent.rotation[face_idx_]; }

inline const geom::Vector3d &Face::normal() const { return cells_->faces.normal[face_idx_]; }

inline const geom::Vector3d &Face::center() const { return cells_->faces.center[face_idx_]; }

template <int dim >
Side<dim> Face::side() const {
    index_t cell_idx = cells_->faces.adjacent.basic[face_idx_];
    return face_idx_ - cells_->faces.offsets[cell_idx];
}

inline double Face::area() const { return cells_->faces.area[face_idx_]; }

inline geom::Vector3d Face::area_n() const { return cells_->faces.area[face_idx_] * cells_->faces.normal[face_idx_]; }

inline double Face::area(bool axial) const { return cells_->faces.get_area(face_idx_, axial); }

inline double Face::area_as() const { return cells_->faces.area_alt[face_idx_]; }

inline int Face::n_vertices() const { return cells_->faces.n_vertices(face_idx_); }

inline index_t Face::vertex_index(int idx) const { return cells_->faces.vertices[face_idx_][idx]; }

inline index_t Face::node_index(int idx) const {
    index_t cell_idx = cells_->faces.adjacent.basic[face_idx_];
    return cells_->verts.offsets[cell_idx] + static_cast<int>(cells_->faces.vertices[face_idx_][idx]);
}

inline geom::Vector3d Face::vs(int idx) const { return cells_->verts[node_index(idx)]; }

inline geom::Vector3d Face::symm_point(const Vector3d &p) const {
    return cells_->faces.symm_point(face_idx_, p);
}

inline int Face::adj_rank() const { return cells_->faces.adjacent.rank[face_idx_]; }

inline index_t Face::adj_index() const { return cells_->faces.adjacent.index[face_idx_]; }

inline index_t Face::adj_ghost() const { return cells_->faces.adjacent.ghost[face_idx_]; }

inline index_t Face::adj_basic() const { return cells_->faces.adjacent.basic[face_idx_]; }

inline bool Face::local_neib() const {
    return cells_->faces.adjacent.is_local(face_idx_);
}

inline Cell Face::neib() const {
    if (utils::mpi::single() || local_neib()) {
        return {cells_, adj_index(), ghosts_};
    }
    return {ghosts_, adj_ghost(), ghosts_};
}

template <typename T>
const T& Face::neib(Storable<T> type) const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->data.get_val(type, adj_index());
    }
    return ghosts_->data.get_val(type, adj_ghost());
}

template <typename T>
std::span<const T> Face::neib(Storable<T[]> type) const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->data.get_val(type, adj_index());
    }
    return ghosts_->data.get_val(type, adj_ghost());
}

template <typename T, size_t N>
std::span<const T, N> Face::neib(Storable<T[N]> type) const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->data.get_val(type, adj_index());
    }
    return ghosts_->data.get_val(type, adj_ghost());
}

inline int Face::neib_flag() const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->flag[adj_index()];
    }
    return ghosts_->flag[adj_ghost()];
}

inline geom::Vector3d Face::neib_center() const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->center[adj_index()];
    }
    return ghosts_->center[adj_ghost()];
}

inline double Face::neib_volume() const {
    if (utils::mpi::single() || local_neib()) {
        return cells_->volume[adj_index()];
    }
    return ghosts_->volume[adj_ghost()];
}

inline double Face::neib_volume(bool axial) const {
    if (axial) {
        if (utils::mpi::single() || local_neib()) {
            return cells_->volume_alt[adj_index()];
        }
        return ghosts_->volume_alt[adj_ghost()];
    }
    return neib_volume();
}

// ================================================================================================
//                                     inline функции Cell
// ================================================================================================

inline int Cell::dim() const { return cells_->dim(); }

inline bool Cell::adaptive() const { return cells_->adaptive(); }

inline int Cell::rank() const { return cells_->rank[index_]; }

inline int Cell::flag() const { return cells_->flag[index_]; }

inline int Cell::level() const { return cells_->level[index_]; }

inline index_t Cell::b_idx() const { return cells_->b_idx[index_]; }

inline index_t Cell::z_idx() const { return cells_->z_idx[index_]; }

inline index_t Cell::id() const { return index_; }

inline index_t Cell::index() const { return cells_->index[index_]; }

inline index_t Cell::next() const { return cells_->next[index_]; }

inline void Cell::set_rank(int rank) const { cells_->rank[index_] = rank; }

inline void Cell::set_flag(int flag) const { cells_->flag[index_] = flag; }

inline const geom::Vector3d& Cell::center() const { return cells_->center[index_]; }

inline double Cell::volume() const { return cells_->volume[index_]; }

inline double Cell::volume(bool axial) const { return cells_->get_volume(index_, axial); }

inline double Cell::volume_as() const { return cells_->volume_alt[index_]; }

inline double Cell::hx() const { return cells_->hx(index_); }

inline double Cell::hy() const { return cells_->hy(index_); }

inline double Cell::hz() const { return cells_->hz(index_); }

inline double Cell::linear_size() const { return cells_->linear_size(index_); }

inline double Cell::incircle_diameter() const { return cells_->incircle_diameter(index_); }

template <typename T>
T& Cell::operator[](Storable<T> var) { return cells_->data.get_val<T>(var, index_); }

template <typename T>
std::span<T> Cell::operator[](Storable<T[]> var) { return cells_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<T, N> Cell::operator[](Storable<T[N]> var) { return cells_->data.get_val<T>(var, index_); }

template <typename T>
const T& Cell::operator[](Storable<T> var) const { return cells_->data.get_val<T>(var, index_); }

template <typename T>
std::span<const T> Cell::operator[](Storable<T[]> var) const { return cells_->data.get_val<T>(var, index_); }

template <typename T, size_t N>
std::span<const T, N> Cell::operator[](Storable<T[N]> var) const { return cells_->data.get_val<T>(var, index_); }

inline void Cell::copy_data_to(Cell &dst_cell) const {
    cells_->copy_data(index_, dst_cell.cells_, dst_cell.index_);
}

inline int Cell::face_count() const { return cells_->face_count(index_); }

inline Face Cell::face(int idx) const {
    return {cells_, cells_->faces.offsets[index_] + idx, ghosts_};
}

inline Face Cell::face(Side2D s) const {
    return {cells_, cells_->faces.offsets[index_] + s, ghosts_};
}

inline Face Cell::face(Side3D s) const {
    return {cells_, cells_->faces.offsets[index_] + s, ghosts_};
}

inline bool Cell::simple_face(Side2D s) const { return cells_->faces.is_simple(index_, s); }

inline bool Cell::simple_face(Side3D s) const { return cells_->faces.is_simple(index_, s); }

inline bool Cell::complex_face(Side2D s) const { return cells_->faces.is_complex(index_, s); }

inline bool Cell::complex_face(Side3D s) const { return cells_->faces.is_complex(index_, s); }

inline FacesRange Cell::faces(Direction dir) const { return {cells_, index_, ghosts_, dir}; }

inline int Cell::node_count() const { return cells_->verts.count(index_); }

inline int Cell::node_rank(int iv) const {
    if (cells_->verts.has_nodes()) {
        return cells_->verts.rank[cells_->verts.offsets[index_] + iv];
    }
    return utils::mpi::rank();
}

inline int Cell::node_index(int iv) const {
    if (cells_->verts.has_nodes()) {
        return cells_->verts.index[cells_->verts.offsets[index_] + iv];
    }
    return -13;
}

inline int Cell::node_ghost(int iv) const {
    if (cells_->verts.has_nodes()) {
        return cells_->verts.ghost[cells_->verts.offsets[index_] + iv];
    }
    return -1;
}

inline const geom::Vector3d* Cell::vertices_data() const { return cells_->verts.coords_data(index_); }

inline double Cell::approx_vol_fraction(const SpFunction& inside) const {
    return cells_->approx_vol_fraction(index_, inside);
}

inline double Cell::volume_fraction(const SpFunction& inside, int n_points) const {
    return cells_->volume_fraction(index_, inside, n_points);
}

inline bool Cell::const_function(const SpFunction& func) const {
    return cells_->const_function(index_, func);
}

inline double Cell::integrate_low(const SpFunction& func, int n_points) const {
    return cells_->integrate_low(index_, func, n_points);
}

} // namespace zephyr::mesh