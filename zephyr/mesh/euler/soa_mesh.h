#pragma once

#include <vector>
#include <boost/next_prior.hpp>
#include <boost/mpl/next_prior.hpp>

#include <zephyr/utils/storage.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/utils/range.h>
#include <zephyr/mesh/primitives/bface.h>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>

namespace zephyr::mesh {
class SoaMesh;
class SoaCell;
class QCell;

using geom::generator::Rectangle;
using geom::generator::Cuboid;

using index_t = int;  // Для индексации примитивов

/// @brief Классная грань
class QFace final {
public:
    using Vector3d = zephyr::geom::Vector3d;
    using Boundary = zephyr::geom::Boundary;

    SoaCell* m_cells;
    index_t iface;

public:
    QFace(SoaCell* cell, index_t iface)
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
    SoaCell* m_cells;
    index_t iface;
    index_t face_end;
    Direction m_dir;

public:
    /// @brief Изолированная грань на стороне side,
    /// не позволяет обходить грани
    FaceIt(SoaCell* cells, index_t iface, index_t face_end, Direction dir = Direction::ANY)
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
    SoaCell* m_cells;
    index_t iface_beg, iface_end;
    Direction m_dir;

    FacesIts(SoaCell* cells,
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
    SoaCell* m_cells;  //< Указатель на сетку
    index_t  m_index;  //< Индекс ячейки

    // Нулевые значения валидны, при таких значениях нельзя перейти к соседней ячейке
    SoaCell* m_locals = nullptr;
    SoaCell* m_aliens = nullptr;

public:
    QCell(SoaCell* cells, index_t index,
          SoaCell* locals = nullptr,
          SoaCell* aliens = nullptr)
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
};

/// @brief Итератор по ячейкам SoaMesh
struct CellIt final {
private:
    SoaCell* m_cells;   //< Указатель на сетку
    index_t  m_index;  //< Индекс ячейки

    // Нулевые значения валидны, при таких значениях нельзя перейти к соседней ячейке
    SoaCell* m_locals = nullptr;
    SoaCell* m_aliens = nullptr;

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = index_t;
    using value_type = QCell;
    using pointer    = QCell*;
    using reference  = QCell&;

    // Конструктор
    CellIt(SoaCell* cells, index_t index,
        SoaCell* locals = nullptr,
        SoaCell* aliens = nullptr)
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


/// @brief Аналог AmrCell развернутый в структуру массивов,
/// теоретически, все функции тоже можно просто скопировать.
/// Три типа сеток.
/// 1. Двумерная AMR сетка, по 8 граней, по 9 вершин на ячейку.
/// 2. Трехмерная AMR сетка, по 24 грани, по 27 вершин на ячейку.
/// 3. Неструктурированная/произвольная сетка. Произвольное число
/// граней и вершин на ячейку, но вершины не склеиваются (не уникальны).
/// На сетки третьего типа пока забьем.
class SoaCell final {
public:
    SoaCell() = default;

    // Копируется из EuMesh
    SoaCell(AmrStorage &locals);

    // Число ячеек
    index_t m_size;

    inline bool empty() const { return m_size == 0; }

    /// @brief Число ячеек
    inline index_t size() const { return m_size; }
    inline index_t n_cells() const { return m_size; }

    void initialize(AmrStorage &locals);

    void initialize(const Rectangle& gen);

    void initialize(const Cuboid& gen);

    // Общие данные ячеек

    int dim;       ///< Размерность ячейки
    int faces_per_cell;  ///< 8 или 24
    int nodes_per_cell;  ///< 9 или 27

    bool adaptive;  ///< Адаптивная ячейка?
    bool linear;    ///< Линейная ячейка?
    bool axial;     ///< Осевая симметрия, см описание класса

    // Тип Element

    std::vector<index_t> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)

    std::vector<int> rank;       ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> index;  ///< Глобальный индекс элемента в локальном Storage (< 0 для неактивных, неопределенных элементов, элементов на удаление)

    // Величины, связанные с адаптацией

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


    int face_count(index_t ic) const {
        int count = 0;
        for (index_t iface: faces_range(ic)) {
            if (faces.is_actual(iface)) {
                ++count;
            }
        }
        return count;
    }

    inline int max_faces(index_t ic) const {
        return face_begin[ic + 1] - face_begin[ic];
    }

    inline int max_nodes(index_t ic) const {
        return node_begin[ic + 1] - node_begin[ic];
    }

    inline double get_volume(index_t ic, bool axial) const {
        return axial ? volume_alt[ic] : volume[ic];
    }

    inline double linear_size(index_t ic) const {
        return dim < 3 ? std::sqrt(volume[ic]) : std::cbrt(volume[ic]);
    }


    struct SoaAdjacent final {
        std::vector<int> rank;
        std::vector<index_t> index;
        std::vector<index_t> alien;

        void resize(index_t n_faces) {
            rank.resize(n_faces);
            index.resize(n_faces);
            alien.resize(n_faces);
        }
    };

    /// @brief Аналог BFaces развернутый в структуру массивов
    struct SoaFace final {
        std::vector<Boundary> boundary;  ///< Тип граничного условия
        std::vector<Vector3d> normal;    ///< Внешняя нормаль к грани
        std::vector<Vector3d> center;    ///< Барицентр грани
        std::vector<double>   area;      ///< Площадь грани
        std::vector<double>   area_alt;  ///< "Альтернативная" площадь грани

        /// @brief Список индексов вершин в массиве вершин ячейки
        std::vector<std::array<int, BFace::max_size>> vertices;

        SoaAdjacent adjacent;

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
        inline bool is_boundary(index_t iface) const {
            return boundary[iface] != Boundary::ORDINARY &&
                   boundary[iface] != Boundary::PERIODIC &&
                   boundary[iface] != Boundary::UNDEFINED;
        }

        /// @brief Является ли грань актуальной?
        inline bool is_actual(index_t iface) const {
            return boundary[iface] != Boundary::UNDEFINED;
        }

        /// @return 'true', если грань не актуальна
        inline bool is_undefined(index_t iface) const {
            return boundary[iface] == Boundary::UNDEFINED;
        }

        /// @brief Установить неопределенную грань
        inline void set_undefined(index_t iface) {
            boundary[iface] = Boundary::UNDEFINED;
            adjacent.rank[iface]  = -1;
            adjacent.alien[iface] = -1;
            adjacent.index[iface] = -1;
        }

        /// @brief Внешняя нормаль грани на площадь
        inline Vector3d area_n(index_t iface) const { return area[iface] * normal[iface]; }

        /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
        inline double get_area(index_t iface, bool axial = false) const {
            return axial ? area_alt[iface] : area[iface];
        }

        inline Vector3d symm_point(index_t iface, const Vector3d& p) const {
            return p + 2.0 * (center[iface] - p).dot(normal[iface]) * normal[iface];
        }


        /// @brief Пропустить грань?
        /// @return 'true' если грань неопределена или
        /// не совпадает с выбраным направлением
        inline bool to_skip(index_t iface, Direction dir) const {
            if (boundary[iface] == Boundary::UNDEFINED) {
                return true;
            }
            switch (dir) {
                case Direction::ANY:
                    return false;
                case Direction::X:
                    return std::abs(normal[iface].x()) < 0.7;
                case Direction::Y:
                    return std::abs(normal[iface].y()) < 0.7;
                case Direction::Z:
                    return std::abs(normal[iface].z()) < 0.7;
                default:
                    return false;
            }
        }

        /// Добавить грани, соответствующие ячейке
        void insert(index_t iface, CellType ctype, int count = -1);

        /// @brief Возвращает индекс из массива locals или aliens
        inline index_t neib_index(index_t iface) const {
            return adjacent.alien[iface] < 0 ? adjacent.index[iface] : adjacent.alien[iface];
        }

        /// @brief Возвращает ссылку на хранилище (locals или aliens), в котором
        /// хранится соседняя ячейка
        inline SoaCell& neib_storage(index_t iface, SoaCell& locals, SoaCell& aliens) const {
            return adjacent.alien[iface] < 0 ? locals : aliens;
        }

        inline std::tuple<const SoaCell&, index_t> get_neib(index_t iface, SoaCell& locals, SoaCell& aliens) const {
            if (adjacent.alien[iface] < 0) {
                return {locals, adjacent.index[iface]};
            }
            else {
                return {aliens, adjacent.alien[iface]};
            }
        }
    };

    void move_item(index_t ic);

    // Грани ячеек
    SoaFace faces;

    // Вершины ячеек
    std::vector<Vector3d> verts;

    // Данные ячеек
    SoaStorage data;


    template <int dim>
    std::conditional_t<dim < 3, const SqQuad&, const SqCube&>
    get_vertices(index_t ic) const {
        if constexpr (dim < 3) {
            return *reinterpret_cast<const SqQuad*>(verts.data() + node_begin[ic]);
        }
        else {
            return *reinterpret_cast<const SqCube*>(verts.data() + node_begin[ic]);
        }
    }

    template <int dim>
    std::conditional_t<dim < 3, SqQuad&, SqCube&>
    get_vertices(index_t ic) {
        if constexpr (dim < 3) {
            return *reinterpret_cast<SqQuad*>(verts.data() + node_begin[ic]);
        }
        else {
            return *reinterpret_cast<SqCube*>(verts.data() + node_begin[ic]);
        }
    }

    /// @brief Простая грань по стороне?
    template <int dim>
    inline bool simple_face(index_t ic, Side<dim> side) const {
        return faces.is_undefined(face_begin[ic] + side[1]);
    }

    /// @brief Простая грань по стороне?
    template <int dim>
    inline bool complex_face(index_t ic, Side<dim> side) const {
        return faces.is_actual(face_begin[ic] + side[1]);
    }

    CellIt begin() { return {this, 0}; }

    CellIt end() { return {this, size()}; }

    QCell operator[](index_t cell_idx) {
        return {this, cell_idx};
    }

    /// @brief Актуальная ячейка?
    inline bool is_actual(index_t ic) const { return index[ic] >= 0; }

    /// @brief Ячейка к удалению
    inline bool is_undefined(index_t ic) const { return index[ic] < 0; }

    /// @brief Устанавливает index = -1 (ячейка вне сетки)
    inline void set_undefined(index_t ic) { index[ic] = -1; }

    // ФУНКЦИИ AmrCell по сути

    // Конструкторы ячеек


    /// @brief Двумерная простая
    void add_cell(index_t ic, const Quad &quad);

    /// @brief Двумерная с осевой симметрией
    void add_cell(index_t ic, const Quad &quad, bool axial);

    /// @brief Двумерная криволинейная
    void add_cell(index_t ic, const SqQuad &quad);

    /// @brief Двумерная криволинейная с осевой симметрией
    void add_cell(index_t ic, const SqQuad &quad, bool axial);

    /// @brief Трехмерная простая
    void add_cell(index_t ic, const Cube &cube);

    /// @brief Трехмерная криволинейная ячейка
    void add_cell(index_t ic, const SqCube &cube);

    /// @brief Двумерная полигональная ячейка. Не адаптивная ячейка, может
    /// представлять четырехугольник, но вершины упорядочены иначе.
    void add_cell(index_t ic, const Polygon &poly);

    /// @brief Трехмрный многогранник. Не адаптивная ячейка, может
    /// представлять куб, но вершины и грани упорядочены иначе.
    void add_cell(index_t ic, Polyhedron poly);


    void print_info(index_t ic) const;

    void visualize(index_t ic, std::string filename) const;

    int check_geometry(index_t ic) const;

    int check_base_face_orientation(index_t ic) const;

    int check_base_vertices_order(index_t ic) const;

    int check_complex_faces(index_t ic) const;

    int check_connectivity(index_t ic) const;

    int check_connectivity(index_t ic, SoaCell& aliens) const;

    utils::range<index_t> faces_range(index_t ic) const {
        return utils::range(face_begin[ic], face_begin[ic + 1]);
    }

    utils::range<index_t> nodes_range(index_t ic) const {
        return utils::range(node_begin[ic], node_begin[ic + 1]);
    }


    void resize(index_t n_cells);
};

class SoaMesh {
public:
    SoaCell m_locals;
    SoaCell m_aliens;

    int m_max_level = 0;
    Distributor distributor;

    /// @brief Установить максимальный допустимый уровень адаптации (>= 0)
    void set_max_level(int max_level);

    int max_level() const;

    bool is_adaptive() const;

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
    Storable<T> add_data(const std::string& name) {
        auto res1 = m_locals.data.add<T>(name);
        auto res2 = m_aliens.data.add<T>(name);
        if (res1.idx != res2.idx) {
            throw std::runtime_error("add_data failed");
        }
        return res1;
    }



    // Преобразовать из классической сетки
    SoaMesh(EuMesh &mesh);

    SoaMesh(Generator& gen);

    SoaMesh(const Rectangle& gen);

    SoaMesh(const Cuboid& gen);

    void to_eu_mesh(EuMesh &mesh) const;


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

}