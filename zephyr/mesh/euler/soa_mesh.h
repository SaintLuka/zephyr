#pragma once

#include <vector>
#include <boost/next_prior.hpp>
#include <boost/mpl/next_prior.hpp>

#include <zephyr/utils/storage.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/utils/range.h>
#include <zephyr/mesh/primitives/bface.h>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh {
class SoaMesh;
class AmrCells;
class QCell;

using geom::generator::Strip;
using geom::generator::Rectangle;
using geom::generator::Cuboid;

/// @brief Классная грань
class QFace final {
public:
    using Vector3d = zephyr::geom::Vector3d;
    using Boundary = zephyr::geom::Boundary;

    AmrCells* m_cells;
    index_t iface;

public:
    QFace(AmrCells* cell, index_t iface)
        : m_cells(cell), iface(iface) { }

    /// @brief Является ли грань граничной?
    bool is_boundary() const;

    /// @brief Является ли грань актуальной?
    bool is_actual() const;

    /// @return 'true', если грань не актуальна
    bool is_undefined() const;

    /// @brief Установить неопределенную грань
    void set_undefined();

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    QCell neib() const;

    template <typename T>
    inline const T& neib(Storable<T> type) const;

    int neib_rank() const;

    /// @brief Индекс соседа, мне для soa
    index_t adj_index() const;

    index_t adj_alien() const;

    Vector3d neib_center() const;

    /// @brief Флаг граничных условий
    Boundary flag() const;

    /// @brief Внешняя нормаль
    const Vector3d& normal() const;

    /// @brief Барицентр грани
    const Vector3d& center() const;

    /// @brief Координата x центра грани
    double x() const;

    /// @brief Координата y центра грани
    double y() const;

    /// @brief Координата z центра грани
    double z() const;

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double area() const;

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    double area(bool axial) const;

    /// @brief Точка симметричная относительно грани
    Vector3d symm_point(const Vector3d& p) const;
};

/// @brief Обертка для типа geom::BFace, реализует интерфейс грани,
/// необходимый для работы. Также содержит несколько новых функций.
class FaceIt {
public:
    AmrCells* m_cells;
    index_t iface;
    index_t face_end;
    Direction m_dir;

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    FaceIt(AmrCells* cells, index_t iface, index_t face_end, Direction dir = Direction::ANY)
            : m_cells(cells), iface(iface), face_end(face_end), m_dir(dir) {
        while (iface < face_end && to_skip(m_dir)) {
            iface += 1;
        }
    }

    FaceIt &operator++() {
        do {
            // Доработать, не забыть
            iface += 1;
        } while (iface < face_end && to_skip(m_dir));
        return *this;
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает с выбраным направлением
    bool to_skip(Direction dir) const;

    bool operator!=(const FaceIt &face) const {
        return iface != face.iface;
    }

    template <int dim>
    QFace operator[](Side<dim> s) const {
        return QFace(m_cells, iface + s);
    }

    // Лайфках, начало структур совпадает
    QFace &operator*() { return *reinterpret_cast<QFace*>(this); }
    const QFace &operator*() const { return *reinterpret_cast<const QFace*>(this); }
};

/// @brief Интерфейс для итерирования по граням ячейки
struct FacesIts {
public:
    AmrCells* m_cells;
    index_t iface_beg, iface_end;
    Direction m_dir;

    FacesIts(AmrCells* cells,
             index_t iface_beg,
             index_t iface_end,
             Direction dir = Direction::ANY)
            : m_cells(cells), iface_beg(iface_beg),
              iface_end(iface_end), m_dir(dir) { }

    FaceIt begin() const { return FaceIt(m_cells, iface_beg, iface_end, m_dir); }

    FaceIt end() const { return FaceIt(m_cells, iface_end, iface_end, m_dir); }

    QFace operator[](int idx) const { return QFace(m_cells, iface_beg + idx); }

    template <int dim>
    QFace operator[](Side<dim> s) const { return QFace(m_cells, iface_beg + s); }

    template <int dim>
    bool complex(Side<dim> s) const;


};

/// @brief Классная ячейка
class QCell final {
public:
    AmrCells* m_cells;  //< Указатель на сетку
    index_t  m_index;  //< Индекс ячейки

    // Нулевые значения валидны, при таких значениях нельзя перейти к соседней ячейке
    AmrCells* m_locals = nullptr;
    AmrCells* m_aliens = nullptr;

public:
    QCell(AmrCells* cells, index_t index,
          AmrCells* locals = nullptr,
          AmrCells* aliens = nullptr)
        : m_cells(cells), m_index(index),
          m_locals(locals), m_aliens(aliens) { }

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
    inline T& operator()(Storable<T> type);

    inline QFace face(int idx) const;

    inline QFace face(Side2D s) const;
    inline QFace face(Side3D s) const;

    inline bool complex_face(Side2D s) const;

    inline bool complex_face(Side3D s) const;



    inline const Vector3d& center() const;

    inline double volume() const;

    inline double volume(bool axial) const;

    inline double linear_size() const;

    double diameter() const;

    inline void set_flag(int flag);

    template <int dim>
    inline typename std::conditional<dim < 3, const SqQuad&, const SqCube&>::type
    get_vertices() const;

    FacesIts faces(Direction dir = Direction::ANY);

    void copy_data_to(QCell& dst);
};

/// @brief Итератор по ячейкам SoaMesh
struct CellIt final {
private:
    AmrCells* m_cells;   //< Указатель на сетку
    index_t  m_index;  //< Индекс ячейки

    // Нулевые значения валидны, при таких значениях нельзя перейти к соседней ячейке
    AmrCells* m_locals = nullptr;
    AmrCells* m_aliens = nullptr;

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = QCell;
    using pointer    = QCell*;
    using reference  = QCell&;

    // Конструктор
    CellIt(AmrCells* cells, index_t index,
           AmrCells* locals = nullptr,
           AmrCells* aliens = nullptr)
        : m_cells(cells), m_index(index),
          m_locals(locals), m_aliens(aliens) { }

    /// @brief Разыменование итератора (для цикла for), лайфхак
    inline QCell& operator*() { return *reinterpret_cast<QCell*>(this); }

    /// @brief Инкремент
    inline CellIt &operator++() { ++m_index; return *this; }

    /// @brief Декремент
    inline CellIt &operator--() { --m_index; return *this; }

    /// @brief Итератор через step
    inline CellIt &operator+=(index_t step) { m_index += step; return *this; }

    /// @brief Итератор через step
    inline CellIt operator+(index_t step) const { return CellIt(m_cells, m_index + step); }

    /// @brief Оператор доступа как для указателя
    /// (интерфейс random access iterator)
    inline QCell operator[](index_t offset) const {
        return QCell(m_cells, m_index + offset, m_locals, m_aliens);
    }

    /// @brief Расстояние между двумя ячейками
    inline index_t operator-(const CellIt& cell) const { return m_index - cell.m_index; }

    /// @brief Оператор сравнения
    inline bool operator<(const CellIt& cell) const { return m_index < cell.m_index; }

    /// @brief Оператор сравнения
    inline bool operator!=(const CellIt& cell) const { return m_index != cell.m_index; }

    /// @brief Оператор сравнения
    inline bool operator==(const CellIt& cell) const { return m_index == cell.m_index; }
};



class SoaMesh {
public:
    AmrCells m_locals;
    AmrCells m_aliens;

    int m_max_level = 0;
    Distributor distributor;

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

    void balance_flags();

    void apply_flags();

    void refine();

    int check_base() const;

    int check_refined() const;

    CellIt begin() { return {&m_locals, 0, &m_locals, &m_aliens}; }

    CellIt end() { return {&m_locals, m_locals.size(), &m_locals, &m_aliens}; }

    QCell operator[](index_t cell_idx) {
        return {&m_locals, cell_idx, &m_locals, &m_aliens};
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

    // Вспомогательная функция для применения stringToDouble к каждому элементу
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




    SoaMesh(const Strip& gen);

    SoaMesh(const Rectangle& gen);

    SoaMesh(const Cuboid& gen);

    Box bbox() const;


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
};

// Определения inline функций


inline int QCell::dim() const { return m_cells->dim; }

inline bool QCell::adaptive() const { return m_cells->adaptive; }

inline int QCell::rank() const { return m_cells->rank[m_index]; }

inline int QCell::flag() const { return m_cells->flag[m_index]; }

inline int QCell::level() const { return m_cells->level[m_index]; }

inline index_t QCell::b_idx() const { return m_cells->b_idx[m_index]; }

inline index_t QCell::z_idx() const { return m_cells->z_idx[m_index]; }

inline index_t QCell::index() const { return m_cells->index[m_index]; }

inline index_t QCell::next() const { return m_cells->next[m_index]; }

inline void QCell::set_flag(int flag) { m_cells->flag[m_index] = flag; }

template <typename T>
inline const T& QFace::neib(Storable<T> type) const {
    return m_cells->data(type)[m_cells->faces.adjacent.index[iface]];
};

template <typename T>
inline T& QCell::operator()(Storable<T> type) {
    return m_cells->data(type)[m_index];
}

inline const Vector3d& QCell::center() const {
    return m_cells->center[m_index];
}

inline double QCell::volume() const {
    return m_cells->volume[m_index];
}

inline double QCell::volume(bool axial) const {
    return m_cells->get_volume(m_index, axial);
}

inline double QCell::linear_size() const {
    return m_cells->linear_size(m_index);
}

template <int dim>
inline
typename std::conditional<dim < 3, const SqQuad&, const SqCube&>::type
QCell::get_vertices() const {
    return m_cells->get_vertices<dim>(m_index);
}

template <int dim>
bool FacesIts::complex(Side<dim> s) const {
    return m_cells->faces.is_actual(iface_beg + s[1]);
}

inline QFace QCell::face(int idx) const {
    return QFace(m_cells, m_cells->face_begin[m_index] + idx);
};


inline QFace QCell::face(Side2D s) const {
    return QFace(m_cells, m_cells->face_begin[m_index] + s);
}

inline QFace QCell::face(Side3D s) const {
    return QFace(m_cells, m_cells->face_begin[m_index] + s);
}

inline bool QCell::complex_face(Side2D s) const {
    return m_cells->faces.is_actual(m_cells->face_begin[m_index] + s[1]);
}

inline bool QCell::complex_face(Side3D s) const {
    return m_cells->faces.is_actual(m_cells->face_begin[m_index] + s[1]);
}


class SoaChildren {
public:

    /// @brief Индексы дочерних ячеек
    std::array<index_t, 8> index = {-1, -1, -1, -1, -1, -1, -1, -1};

    /// @brief Обязательная инициализация локального хранилища
    explicit SoaChildren(AmrCells& locals) : m_locals(locals) { }


    /// @brief Размерность определяется по числу дочерних ячеек
    inline int dim() const { return index[4] < 0 ? 2 : 3; }

    /// @brief Число дочерних ячеек (4 для 2D и 8 для 3D)
    inline int count() const { return index[4] < 0 ? 4 : 8; }


    /// @brief Получить дочернюю ячейку
    inline QCell operator[](int idx) {
        return QCell(&m_locals, index[idx]);
    };

    struct iterator {
        iterator(SoaChildren& children, int idx)
            : m_children(children), m_idx(idx) { }

        QCell operator*() { return m_children[m_idx]; }

        void operator++() { ++m_idx; }

        bool operator!=(const iterator &it) const { return m_idx != it.m_idx; }

    private:
        int m_idx;
        SoaChildren& m_children;
    };

    iterator begin() {
        return iterator(*this, 0);
    }

    iterator end() {
        return iterator(*this, count());
    }

protected:
    /// @brief Дочерние ячейки расположены в локальном хранилище
    AmrCells& m_locals;
};

}