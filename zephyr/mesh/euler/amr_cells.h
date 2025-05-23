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

#include <zephyr/mesh/euler/amr_faces.h>

namespace zephyr::mesh {

class SoaMesh;

class AmrCells;

class QCell;
class CellIt;

using geom::generator::Strip;
using geom::generator::Rectangle;
using geom::generator::Cuboid;


/// @brief Аналог AmrCell развернутый в структуру массивов,
/// теоретически, все функции тоже можно просто скопировать.
/// Три типа сеток.
/// 1. Двумерная AMR сетка, по 8 граней, по 9 вершин на ячейку.
/// 2. Трехмерная AMR сетка, по 24 грани, по 27 вершин на ячейку.
/// 3. Неструктурированная/произвольная сетка. Произвольное число
/// граней и вершин на ячейку, но вершины не склеиваются или не уникальны.
/// На сетки третьего типа пока забьем.
class AmrCells final {
    using AmrVerts = std::vector<Vector3d>;

public:

    index_t m_size = 0;     ///< Число ячеек

    // Общие свойства ячеек

    int  dim      = -1;     ///< Размерность ячейки
    bool adaptive = false;  ///< Адаптивная ячейка?
    bool linear   = true;   ///< Линейная ячейка?
    bool axial    = false;  ///< Осевая симметрия?

    // Главные индексы ячейки

    std::vector<int>     rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
    std::vector<index_t> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)
    std::vector<index_t> index;  ///< Глобальный индекс элемента в локальном Storage
                                 /// ( < 0 для неопределенных элементов, элементов на удаление)

    // Величины, связанные с адаптацией

    std::vector<int> flag;       ///< Желаемый флаг адаптации
    std::vector<int> level;      ///< Уровень адаптации (0 для базовой)
    std::vector<index_t> b_idx;  ///< Индекс среди базовых ячеек
    std::vector<index_t> z_idx;  ///< Индекс ячейки на z-кривой

    // Геометрия ячеек

    std::vector<Vector3d> center;      ///< Барицентр ячейки
    std::vector<double>   volume;      ///< Объем трехмерной / площадь двумерной ячейки
    std::vector<double>   volume_alt;  ///< Объем двумерной осесимметричной ячейки

    // Грани и вершины ячеек

    std::vector<index_t> face_begin;   ///< Индекс первой грани ячейки
    std::vector<index_t> node_begin;   ///< Индекс первой вершины ячейки

    AmrFaces   faces;  ///< Массив граней ячеек
    AmrVerts   verts;  ///< Массив вершин ячеек

    SoaStorage data;   ///< SoA с данными ячеек





    AmrCells() = default;

    /// @brief Пустое хранилище?
    inline bool empty() const { return m_size == 0; }

    /// @brief Число ячеек
    inline index_t size() const { return m_size; }

    /// @brief Число ячеек
    inline index_t n_cells() const { return m_size; }



    void initialize(const Strip& gen);

    void initialize(const Rectangle& gen);

    void initialize(const Cuboid& gen);


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


    void move_item(index_t ic);

    /// @brief Скопировать все данные с индекса from на индекс to.
    void copy_data(index_t from, index_t to);

    /// @brief Скопировать все данные целиком с индекса from,
    /// в хранилище dst на индекс to
    void copy_data(index_t from, AmrCells* dst, index_t to);



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

    CellIt begin();

    CellIt end();

    QCell operator[](index_t cell_idx);

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

    int check_connectivity(index_t ic, AmrCells& aliens) const;

    range<index_t> faces_range(index_t ic) const {
        return utils::range(face_begin[ic], face_begin[ic + 1]);
    }

    range<index_t> nodes_range(index_t ic) const {
        return utils::range(node_begin[ic], node_begin[ic + 1]);
    }


    void resize(index_t n_cells);
};

} // namespace zephyr::mesh