#pragma once

#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh {

class EuCell;

/// @brief Классная грань
class EuFace final {
    using Vector3d = geom::Vector3d;
    using Boundary = geom::Boundary;

    AmrCells* m_cells;     //< Указатель на сетку (обычно locals)
    index_t   m_face_idx;  //< Индекс первой грани

    /// @brief Нулевое значение допускается, но в целях оптимизации проверок
    /// нигде нет, так что ловите segfaults.
    AmrCells* m_aliens = nullptr;

public:
    /// @brief Основной конструктор
    EuFace(AmrCells* cells, index_t face_idx, AmrCells* aliens = nullptr)
        : m_cells(cells), m_face_idx(face_idx), m_aliens(aliens) { }

    /// @brief Флаг граничных условий
    Boundary flag() const;

    /// @brief Является ли грань граничной?
    bool is_boundary() const;

    /// @brief Является ли грань актуальной?
    bool is_actual() const;

    /// @return 'true', если грань не актуальна
    bool is_undefined() const;

    /// @brief Установить неопределенную грань
    void set_undefined();


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

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double area() const;

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double area(bool axial) const;

    /// @brief Точка симметричная относительно грани
    Vector3d symm_point(const Vector3d& p) const;

    /// @brief Получить вершину грани
    Vector3d vs(int idx) const;



    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    EuCell neib() const;

    template <typename T>
    inline const T& neib(utils::Storable<T> type) const;

    int neib_rank() const;

    /// @brief Индекс соседа, мне для soa
    index_t adj_index() const;

    index_t adj_alien() const;

    Vector3d neib_center() const;
};

/// @brief Обертка для типа geom::BFace, реализует интерфейс грани,
/// необходимый для работы. Также содержит несколько новых функций.
class EuFaceIt {
private:
    AmrCells* m_cells;     //< Указатель на сетку (обычно locals)
    index_t   m_face_idx;  //< Индекс первой грани

    /// @brief Нулевое значение допускается, но в целях оптимизации проверок
    /// нигде нет, так что ловите segfaults.
    AmrCells* m_aliens = nullptr;

    index_t m_face_end;    //< Индекс за последней гранью
    Direction m_dir;       //< Выбранное направление граней

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    EuFaceIt(AmrCells* cells, index_t face_idx, index_t face_end, AmrCells* aliens, Direction dir = Direction::ANY)
            : m_cells(cells), m_face_idx(face_idx), m_aliens(aliens), m_face_end(face_end), m_dir(dir) {
        while (m_face_idx < m_face_end && to_skip(m_dir)) {
            m_face_idx += 1;
        }
    }

    EuFaceIt &operator++() {
        do {
            // Доработать, не забыть (забыл)
            m_face_idx += 1;
        } while (m_face_idx < m_face_end && to_skip(m_dir));
        return *this;
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределенна или
    /// не совпадает с выбранным направлением
    bool to_skip(Direction dir) const {
        return m_cells->faces.to_skip(m_face_idx, dir);
    }

    bool operator!=(const EuFaceIt &face) const {
        return m_face_idx != face.m_face_idx;
    }

    template <int dim>
    EuFace operator[](Side<dim> s) const {
        return EuFace(m_cells, m_face_idx + s);
    }

    // Лайфхак, начало структур совпадает
    EuFace &operator*() { return *reinterpret_cast<EuFace*>(this); }

    // Лайфхак, начало структур совпадает
    const EuFace &operator*() const { return *reinterpret_cast<const EuFace*>(this); }
};

/// @brief Интерфейс для итераций по граням ячейки
class EuFaces {
    EuFaceIt m_begin;
    EuFaceIt m_end;
public:
    EuFaces(
        AmrCells* cells,
        index_t cell_idx,
        AmrCells* aliens = nullptr,
        Direction dir = Direction::ANY)
    :
    m_begin(cells,
            cells->face_begin[cell_idx],
            cells->face_begin[cell_idx + 1],
            aliens, dir),
    m_end(  cells,
            cells->face_begin[cell_idx + 1],
            cells->face_begin[cell_idx + 1],
            aliens, dir) { }

    EuFaceIt begin() const { return m_begin; }

    EuFaceIt end() const { return m_end; }
};

/// @brief Классная ячейка
class EuCell final {
public:
    AmrCells* m_cells;  //< Указатель на сетку (обычно locals)
    index_t   m_index;  //< Индекс ячейки

    /// @brief Нулевое значение допускается, но в целях оптимизации проверок 
    /// нигде нет, так что ловите segfaults.
    AmrCells* m_aliens = nullptr;

public:
    EuCell(AmrCells* cells, index_t index, AmrCells* aliens = nullptr)
        : m_cells(cells), m_index(index), m_aliens(aliens) { }

    /// @brief Размерность ячейки
    inline int dim() const;

    /// @brief Адаптивная (AMR) ячейка?
    inline bool adaptive() const;

    /// @brief Ранг, которому принадлежит ячейка
    inline int rank() const;

    inline index_t b_idx() const;

    inline index_t z_idx() const;

    /// @brief Индекс ячейки на z-кривой
    inline index_t index() const;

    /// @brief Индекс новой ячейки (в алгоритмах)
    inline index_t next() const;

    inline index_t flag() const;

    inline index_t level() const;

    template <typename T>
    T& operator()(utils::Storable<T> type);

    template <typename T>
    const T& operator()(utils::Storable<T> type) const;

    inline EuFace face(int idx) const;

    inline EuFace face(Side2D s) const;

    inline EuFace face(Side3D s) const;

    template <int dim>
    bool complex_face(Side<dim> s) const {
        return m_cells->complex_face(m_index, s);
    }

    int node_count() const { return m_cells->node_count(m_index); }

    const geom::Vector3d* vertices_data() const {
        return m_cells->vertices_data(m_index);
    }


    inline const geom::Vector3d& center() const;

    inline double volume() const;

    inline double volume(bool axial) const;

    inline double linear_size() const;

    double incircle_diameter() const;

    inline void set_flag(int flag);

    template <int dim>
    SqMap<dim> mapping() const { return m_cells->mapping<dim>(m_index); }

    int max_face_count() const;

    EuFaces faces(Direction dir = Direction::ANY) const;

    void copy_data_to(EuCell& dst) const;

    geom::Polygon polygon() const;

    double approx_vol_fraction(const std::function<double(const geom::Vector3d &)> &inside) const {
        return m_cells->approx_vol_fraction(m_index, inside);
    }

    double volume_fraction(const std::function<double(const geom::Vector3d &)> &inside, int n_points) const {
        return m_cells->volume_fraction(m_index, inside, n_points);
    }

    bool const_function(const std::function<double(const geom::Vector3d&)>& func) const {
        return m_cells->const_function(m_index, func);
    }

    double integrate_low(const std::function<double(const geom::Vector3d&)>& func, int n_points) const {
        return m_cells->integrate_low(m_index, func, n_points);
    }
};

/// @brief Итератор по ячейкам из EuMesh или AmrCells
class EuCellIt final {
private:
    AmrCells* m_cells;  //< Указатель на сетку
    index_t   m_index;  //< Индекс ячейки

    /// @brief Нулевое значение допускается, но в целях оптимизации проверок 
    /// нигде нет, так что ловите segfaults.
    AmrCells* m_aliens = nullptr;

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = EuCell;
    using pointer    = EuCell*;
    using reference  = EuCell&;

    // Конструктор
    EuCellIt(AmrCells* cells, index_t index, AmrCells* aliens = nullptr)
        : m_cells(cells), m_index(index), m_aliens(aliens) { }

    /// @brief Разыменование итератора (для цикла for), лайфхак
    EuCell& operator*() { return *reinterpret_cast<EuCell*>(this); }

    /// @brief Инкремент
    EuCellIt &operator++() { ++m_index; return *this; }

    /// @brief Декремент
    EuCellIt &operator--() { --m_index; return *this; }

    /// @brief Итератор через step
    EuCellIt &operator+=(index_t step) { m_index += step; return *this; }

    /// @brief Итератор через step
    EuCellIt operator+(index_t step) const {
        return EuCellIt(m_cells, m_index + step);
    }

    /// @brief Оператор доступа как для указателя (требует random access iterator)
    EuCell operator[](index_t offset) const {
        return EuCell(m_cells, m_index + offset, m_aliens);
    }

    /// @brief Расстояние между двумя ячейками
    index_t operator-(const EuCellIt& cell) const { return m_index - cell.m_index; }

    /// @brief Оператор сравнения
    bool operator<(const EuCellIt& cell) const { return m_index < cell.m_index; }

    /// @brief Оператор сравнения
    bool operator!=(const EuCellIt& cell) const { return m_index != cell.m_index; }

    /// @brief Оператор сравнения
    bool operator==(const EuCellIt& cell) const { return m_index == cell.m_index; }
};

/// @brief Для итераций по AmrCells, aliens = nullptr,
/// поэтому проход по соседям не возможен
inline EuCellIt begin(AmrCells& cells) {
    return EuCellIt(&cells, 0);
}

/// @brief Для итераций по AmrCells, aliens = nullptr,
/// поэтому проход по соседям не возможен
inline EuCellIt end(AmrCells& cells) {
    return EuCellIt(&cells, cells.size());
}

/// @brief Простой интерфейс для обхода дочерних ячеек.
/// Предполагается, что дочерние ячейки располагаются в локальном хранилище.
/// Во время операций split (refine) и merge (coarse) сетка может находиться
/// в не совместном состоянии, поэтому переход по соседям запрещен.
class Children {
public:

    /// @brief Индексы дочерних ячеек
    std::array<index_t, 8> index = {-1, -1, -1, -1, -1, -1, -1, -1};


    /// @brief Обязательная инициализация локального хранилища
    explicit Children(AmrCells* locals) : m_locals(locals) { }

    /// @brief Размерность определяется по числу дочерних ячеек
    int dim() const { return index[4] < 0 ? 2 : 3; }

    /// @brief Число дочерних ячеек (4 для 2D и 8 для 3D)
    int count() const { return index[4] < 0 ? 4 : 8; }

    /// @brief Получить дочернюю ячейку
    EuCell operator[](int idx) const { return EuCell(m_locals, index[idx]); };

    /// @brief Простенький итератор по дочерним ячейкам
    struct iterator {
        iterator(Children& children, int idx)
                : m_children(children), m_idx(idx) { }

        EuCell operator*() const { return m_children[m_idx]; }

        void operator++() { ++m_idx; }

        bool operator!=(const iterator &it) const { return m_idx != it.m_idx; }

    private:
        int m_idx;
        Children& m_children;
    };

    iterator begin() { return iterator(*this, 0); }

    iterator end() { return iterator(*this, count()); }

protected:
    /// @brief Дочерние ячейки расположены в локальном хранилище
    AmrCells* m_locals;
};



// ============================================================================
//                         inline функции EuFace
// ============================================================================

inline geom::Boundary EuFace::flag() const {
    return m_cells->faces.boundary[m_face_idx];
}

inline bool EuFace::is_boundary() const {
    return m_cells->faces.is_boundary(m_face_idx);
}

inline bool EuFace::is_actual() const {
    return m_cells->faces.is_actual(m_face_idx);
}

inline bool EuFace::is_undefined() const {
    return m_cells->faces.is_undefined(m_face_idx);
}

inline void EuFace::set_undefined() {
    m_cells->faces.set_undefined(m_face_idx);
}

inline const geom::Vector3d &EuFace::normal() const {
    return m_cells->faces.normal[m_face_idx];
}

inline const geom::Vector3d &EuFace::center() const {
    return m_cells->faces.center[m_face_idx];
}

inline double EuFace::area() const {
    return m_cells->faces.area[m_face_idx];
}

inline double EuFace::area(bool axial) const {
    return m_cells->faces.get_area(m_face_idx, axial);
}

inline geom::Vector3d EuFace::symm_point(const Vector3d &p) const {
    return m_cells->faces.symm_point(m_face_idx, p);
}

inline geom::Vector3d EuFace::vs(int idx) const {
    index_t cell_idx = m_cells->faces.adjacent.basic[m_face_idx];
    return m_cells->verts[m_cells->node_begin[cell_idx] +
                          m_cells->faces.vertices[m_face_idx][idx]];
}

inline index_t EuFace::adj_index() const {
    return m_cells->faces.adjacent.index[m_face_idx];
}

inline index_t EuFace::adj_alien() const {
    return m_cells->faces.adjacent.alien[m_face_idx];
}

inline EuCell EuFace::neib() const {
    if (m_cells->faces.adjacent.is_local(m_face_idx)) {
        return EuCell(m_cells, adj_index(), m_aliens);
    }
    return EuCell(m_aliens, adj_alien(), m_aliens);
}

template <typename T>
const T& EuFace::neib(utils::Storable<T> type) const {
    if (m_cells->faces.adjacent.is_local(m_face_idx)) {
        return m_cells->data(type)[adj_index()];
    }
    return m_aliens->data(type)[adj_alien()];
}

inline int EuFace::neib_rank() const {
    return m_cells->faces.adjacent.rank[m_face_idx];
}

inline geom::Vector3d EuFace::neib_center() const {
    if (m_cells->faces.adjacent.is_local(m_face_idx)) {
        return m_cells->center[adj_index()];
    }
    return m_aliens->center[adj_alien()];
}

// ============================================================================
//                         inline функции EuCell
// ============================================================================

inline int EuCell::dim() const { return m_cells->dim(); }

inline bool EuCell::adaptive() const { return m_cells->adaptive(); }

inline int EuCell::rank() const { return m_cells->rank[m_index]; }

inline int EuCell::flag() const { return m_cells->flag[m_index]; }

inline int EuCell::level() const { return m_cells->level[m_index]; }

inline index_t EuCell::b_idx() const { return m_cells->b_idx[m_index]; }

inline index_t EuCell::z_idx() const { return m_cells->z_idx[m_index]; }

inline index_t EuCell::index() const { return m_cells->index[m_index]; }

inline index_t EuCell::next() const { return m_cells->next[m_index]; }

inline void EuCell::set_flag(int flag) { m_cells->flag[m_index] = flag; }


template <typename T>
T& EuCell::operator()(utils::Storable<T> type) {
    return m_cells->data(type)[m_index];
}

template <typename T>
const T& EuCell::operator()(utils::Storable<T> type) const {
    return m_cells->data(type)[m_index];
}

inline const geom::Vector3d& EuCell::center() const {
    return m_cells->center[m_index];
}

inline double EuCell::volume() const {
    return m_cells->volume[m_index];
}

inline double EuCell::volume(bool axial) const {
    return m_cells->get_volume(m_index, axial);
}

inline double EuCell::linear_size() const {
    return m_cells->linear_size(m_index);
}

inline EuFace EuCell::face(int idx) const {
    return EuFace(m_cells, m_cells->face_begin[m_index] + idx, m_aliens);
}

inline EuFace EuCell::face(Side2D s) const {
    return EuFace(m_cells, m_cells->face_begin[m_index] + s, m_aliens);
}

inline EuFace EuCell::face(Side3D s) const {
    return EuFace(m_cells, m_cells->face_begin[m_index] + s, m_aliens);
}


inline void EuCell::copy_data_to(EuCell &dst) const {
    m_cells->copy_data(m_index, dst.m_cells, dst.m_index);
}

inline double EuCell::incircle_diameter() const {
    return m_cells->incircle_diameter(m_index);
}

inline int EuCell::max_face_count() const {
    return m_cells->max_face_count(m_index);
}

inline EuFaces EuCell::faces(Direction dir) const {
    return EuFaces(m_cells, m_index, m_aliens, dir);
}

} // namespace zephyr::mesh