#pragma once

#include <vector>
#include <boost/next_prior.hpp>
#include <boost/mpl/next_prior.hpp>

#include <zephyr/utils/storage.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/utils/range.h>
#include <zephyr/mesh/primitives/bface.h>

namespace zephyr::mesh {
class SoaMesh;
class SoaCell;
class QCell;

using index_t = int;  // Для индексации примитивов

/// @brief Классная грань
class QFace final {
public:
    using Vector3d = zephyr::geom::Vector3d;
    using Boundary = zephyr::geom::Boundary;

    SoaCell* m_cells;
    index_t face_idx;

public:
    QFace(SoaCell* cell, index_t face_idx)
        : m_cells(cell), face_idx(face_idx) { }

    /// @brief Является ли грань граничной?
    bool is_boundary() const;

    /// @brief Является ли грань актуальной?
    bool is_actual() const;

    /// @return 'true', если грань не актуальна
    bool is_undefined() const;

    /// @brief Установить неопределенную грань
    void set_undefined();

    const Adjacent &adjacent() const;

    /// @brief Ячейка по внешней нормали, на границе сетки возвращается
    /// сама ячейка
    QCell neib() const;

    /// @brief Индекс соседа, мне для soa
    index_t neib_index() const;

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
    SoaCell* m_cells;
    index_t face_idx;
    index_t face_end;
    Direction m_dir;

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    FaceIt(SoaCell* cells, index_t face_idx, index_t face_end, Direction dir = Direction::ANY)
            : m_cells(cells), face_idx(face_idx), face_end(face_end), m_dir(dir) {
        while (face_idx < face_end && to_skip(m_dir)) {
            face_idx += 1;
        }
    }

    FaceIt &operator++() {
        do {
            // Доработать, не забыть
            face_idx += 1;
        } while (face_idx < face_end && to_skip(m_dir));
        return *this;
    }

    /// @brief Пропустить грань?
    /// @return 'true' если грань неопределена или
    /// не совпадает с выбраным направлением
    bool to_skip(Direction dir) const;

    bool operator!=(const FaceIt &iface) const {
        return face_idx != iface.face_idx;
    }

    QFace operator[](Side s) const {
        return QFace(m_cells, face_idx + s);
    }

    // Лайфках, начало структур совпадает
    QFace &operator*() { return *reinterpret_cast<QFace*>(this); }
    const QFace &operator*() const { return *reinterpret_cast<const QFace*>(this); }
};

/// @brief Интерфейс для итерирования по граням ячейки
struct FacesIts {
public:
    FaceIt m_beg, m_end;

    FacesIts(SoaCell* cells,
             index_t face_idx_beg,
             index_t face_idx_end,
             Direction dir = Direction::ANY)
            : m_beg(cells, face_idx_beg, face_idx_end, dir),
              m_end(cells, face_idx_end, face_idx_end, dir) { }

    FaceIt begin() const { return m_beg; }

    FaceIt end() const { return m_end; }

    QFace operator[](Side s) const {
        return m_beg[s];
    }
};

/// @brief Классная ячейка
struct QCell final {
private:
    SoaCell* m_cells;  //< Указатель на сетку
    index_t   cell_idx;  //< Индекс ячейки

public:
    QCell(SoaCell* cells, index_t cell_idx)
        : m_cells(cells), cell_idx(cell_idx) { }


    /// @brief Размерность ячейки
    inline int dim() const;

    /// @brief Адаптивная (AMR) ячейка?
    inline bool adaptive() const;

    /// @brief Ранг, которому принадлежит ячейка
    inline int rank() const;

    /// @brief Индекс ячейки на z-кривой
    inline index_t index() const;

    /// @brief Индекс новой ячейки (в алгоритмах)
    inline index_t next() const;

    template <typename T>
    inline T& operator()(Storable<T> type);


    inline const Vector3d& center() const;

    inline double volume() const;

    inline double volume(bool axial) const;

    inline double linear_size() const;

    double diameter() const;

    template <int dim>
    inline typename std::conditional<dim < 3, const SqQuad&, const SqCube&>::type
    get_vertices() const;

    FacesIts faces(Direction dir = Direction::ANY);
};

/// @brief Итератор по ячейкам SoaMesh
struct CellIt final {
private:
    SoaCell* m_cells;   //< Указатель на сетку
    index_t   cell_idx;  //< Индекс ячейки

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = QCell;
    using pointer    = QCell*;
    using reference  = QCell&;

    // Конструктор
    CellIt(SoaCell* cells, index_t cell_idx)
       : m_cells(cells), cell_idx(cell_idx) { }

    /// @brief Разыменование итератора (для цикла for)
    inline QCell& operator*() {
        // Грязноватый лайфхак
        return *reinterpret_cast<QCell*>(this);
    }

    /// @brief Инкремент
    inline CellIt &operator++() { ++cell_idx; return *this; }

    /// @brief Декремент
    inline CellIt &operator--() { --cell_idx; return *this; }

    /// @brief Итератор через step
    inline CellIt &operator+=(index_t step) {
        cell_idx += step; return *this;
    }

    /// @brief Итератор через step
    inline CellIt operator+(index_t step) const {
        return CellIt(m_cells, cell_idx + step);
    }

    /// @brief Оператор доступа как для указателя
    /// (интерфейс random access iterator)
    inline QCell operator[](index_t offset) const {
        return QCell(m_cells, cell_idx + offset);
    }

    /// @brief Расстояние между двумя ячейками
    inline ptrdiff_t operator-(const CellIt& cell) const {
        return cell_idx - cell.cell_idx;
    }

    /// @brief Оператор сравнения
    inline bool operator<(const CellIt& cell) const {
        return cell_idx < cell.cell_idx;
    }

    /// @brief Оператор сравнения
    inline bool operator!=(const CellIt& cell) const {
        return cell_idx != cell.cell_idx;
    }

    /// @brief Оператор сравнения
    inline bool operator==(const CellIt& cell) const {
        return cell_idx == cell.cell_idx;
    }
};


/// @brief Аналог AmrCell развернутый в структуру массивов,
/// теоретически, все функции тоже можно просто скопировать.
class SoaCell final {
public:
    // Копируется из EuMesh
    SoaCell(AmrStorage &locals);

    // Число ячеек
    index_t n_cells;

    index_t size() const { return n_cells; }

    // Общие данные ячеек

    int  dim;       ///< Размерность ячейки
    bool adaptive;  ///< Адаптивная ячейка?
    bool linear;    ///< Линейная ячейка?
    bool axial;     ///< Осевая симметрия, см описание класса

    // Тип Element

    std::vector<int> rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> index;  ///< Индекс элемента в локальном Storage (< 0 для неактивных, неопределенных элементов, элементов на удаление)
    std::vector<index_t> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)

    // Величины, связаные с адаптацией

    std::vector<int> flag;   ///< Желаемый флаг адаптации
    std::vector<int> level;  ///< Уровень адаптации (0 для базовой)
    std::vector<index_t> b_idx;  ///< Индекс среди базовых ячеек
    std::vector<index_t> z_idx;  ///< Индекс ячейки на z-кривой

    // Геометрия

    std::vector<Vector3d> center;      ///< Барицентр ячейки
    std::vector<double>   volume;      ///< Объем трехмерной / площадь двумерной ячейки
    std::vector<double>   volume_alt;  ///< Объем двумерной осесимметричной ячейки

    // Индексы граней и узлов

    std::vector<index_t>   face_begin;
    std::vector<index_t>   node_begin;

    /// @brief Аналог BFaces развернутый в структуру массивов
    struct SoaFace {
        std::vector<Boundary> boundary;  ///< Тип граничного условия
        std::vector<Adjacent> adjacent;  ///< Составной индекс смежной ячейки
        std::vector<Vector3d> normal;    ///< Внешняя нормаль к грани
        std::vector<Vector3d> center;    ///< Барицентр грани
        std::vector<double>   area;      ///< Площадь грани
        std::vector<double>   area_alt;  ///< "Альтернативная" площадь грани

        /// @brief Список индексов вершин в массиве вершин ячейки
        std::vector<std::array<int, BFace::max_size>> vertices;

        void resize(index_t n_faces) {
            boundary.resize(n_faces);
            adjacent.resize(n_faces);
            normal.resize(n_faces);
            center.resize(n_faces);
            area.resize(n_faces);
            area_alt.resize(n_faces);
            vertices.resize(n_faces);
        }

        /// @brief Является ли грань граничной?
        inline bool is_boundary(index_t face_idx) const {
            return boundary[face_idx] != Boundary::ORDINARY &&
                   boundary[face_idx] != Boundary::PERIODIC &&
                   boundary[face_idx] != Boundary::UNDEFINED;
        }

        /// @brief Является ли грань актуальной?
        inline bool is_actual(index_t face_idx) const {
            return boundary[face_idx] != Boundary::UNDEFINED;
        }

        /// @return 'true', если грань не актуальна
        inline bool is_undefined(index_t face_idx) const {
            return boundary[face_idx] == Boundary::UNDEFINED;
        }

        /// @brief Установить неопределенную грань
        inline void set_undefined(index_t face_idx) {
            boundary[face_idx] = Boundary::UNDEFINED;
            adjacent[face_idx].rank  = -1;
            adjacent[face_idx].index = -1;
            adjacent[face_idx].alien = -1;
        }

        /// @brief Внешняя нормаль грани на площадь
        inline Vector3d area_n(index_t face_idx) const { return area[face_idx] * normal[face_idx]; }

        /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
        inline double get_area(index_t face_idx, bool axial = false) const {
            return axial ? area_alt[face_idx] : area[face_idx];
        }

        inline Vector3d symm_point(index_t face_idx, const Vector3d& p) const {
            return p + 2.0 * (center[face_idx] - p).dot(normal[face_idx]) * normal[face_idx];
        }

        /// @brief Пропустить грань?
        /// @return 'true' если грань неопределена или
        /// не совпадает с выбраным направлением
        inline bool to_skip(index_t face_idx, Direction dir) const {
            if (boundary[face_idx] == Boundary::UNDEFINED) {
                return true;
            }
            switch (dir) {
                case Direction::ANY:
                    return false;
                case Direction::X:
                    return std::abs(normal[face_idx].x()) < 0.7;
                case Direction::Y:
                    return std::abs(normal[face_idx].y()) < 0.7;
                case Direction::Z:
                    return std::abs(normal[face_idx].z()) < 0.7;
                default:
                    return false;
            }
        }
    };

    // Грани ячеек
    SoaFace faces;

    // Вершины ячеек
    std::vector<Vector3d> verts;

    // Данные ячеек
    SoaStorage data;


    template <int dim>
    typename std::conditional<dim < 3, const SqQuad&, const SqCube&>::type
    get_vertices(index_t ic) const {
        if constexpr (dim < 3) {
            return *reinterpret_cast<const SqQuad*>(verts.data() + node_begin[ic]);
        }
        else {
            return *reinterpret_cast<const SqCube*>(verts.data() + node_begin[ic]);
        }
    }

    CellIt begin() { return {this, 0}; }

    CellIt end() { return {this, n_cells}; }

    QCell operator[](index_t cell_idx) {
        return {this, cell_idx};
    }

    /// @brief Актуальная ячейка?
    inline bool is_actual(index_t ic) const { return index[ic] >= 0; }

    /// @brief Ячейка к удалению
    inline bool is_undefined(index_t ic) const { return index[ic] < 0; }

    /// @brief Устанавливает index = -1 (ячейка вне сетки)
    inline void set_undefined(index_t ic) { index[ic] = -1; }



    void resize(index_t n_cells, index_t faces_per_cell, index_t nodes_per_cell);
};

class SoaMesh {
public:
    SoaCell cells;

    CellIt begin() { return {&cells, 0}; }

    CellIt end() { return {&cells, cells.size()}; }

    QCell operator[](index_t cell_idx) {
        return {&cells, cell_idx};
    }


    template <typename T>
    Storable<T> add_data(const std::string& name) {
        return cells.data.add<T>(name);
    }



    // Преобразовать из классической сетки
    SoaMesh(EuMesh &mesh);

    void to_eu_mesh(EuMesh &mesh) const;


    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        utils::range<index_t> range(cells.size());

        threads::for_each<n_tasks_per_thread>(range.begin(), range.end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each2(Func &&func, Args&&... args) {
        threads::for_each<n_tasks_per_thread>(begin(), end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, index_t&>>
    auto min(Func &&func, const Value &init)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        utils::range<index_t> range(cells.size());
        return threads::min<n_tasks_per_thread>(range.begin(), range.end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, QCell>>
    auto min2(Func &&func, const Value &init)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        return threads::min<n_tasks_per_thread>(begin(), end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, index_t&>>
    auto min(Func &&func)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        utils::range<index_t> range(cells.size());
        return threads::min<n_tasks_per_thread>(range.begin(), range.end(), std::forward<Func>(func));
    }
};

// Определения inline функций


inline int QCell::dim() const { return m_cells->dim; }

inline bool QCell::adaptive() const { return m_cells->adaptive; }

inline int QCell::rank() const { return m_cells->rank[cell_idx]; }

inline index_t QCell::index() const { return m_cells->index[cell_idx]; }

inline index_t QCell::next() const { return m_cells->next[cell_idx]; }

template <typename T>
inline T& QCell::operator()(Storable<T> type) {
    return m_cells->data(type)[cell_idx];
}

inline const Vector3d& QCell::center() const {
    return m_cells->center[cell_idx];
}

inline double QCell::volume() const {
    return m_cells->volume[cell_idx];
}

inline double QCell::volume(bool axial) const {
    return axial ? m_cells->volume_alt[cell_idx] : m_cells->volume[cell_idx];
}

inline double QCell::linear_size() const {
    return m_cells->dim < 3 ?
           std::sqrt(m_cells->volume[cell_idx]) :
           std::cbrt(m_cells->volume[cell_idx]);
}

template <int dim>
inline
typename std::conditional<dim < 3, const SqQuad&, const SqCube&>::type
QCell::get_vertices() const {
    return m_cells->get_vertices<dim>(cell_idx);
}

}