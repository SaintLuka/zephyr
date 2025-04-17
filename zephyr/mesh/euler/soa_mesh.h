#pragma once

#include <vector>
#include <boost/next_prior.hpp>
#include <boost/mpl/next_prior.hpp>

#include <zephyr/utils/storage.h>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/utils/range.h>

namespace zephyr::mesh {

class SoaMesh {
public:
    int  dim;       ///< Размерность ячейки
    bool adaptive;  ///< Адаптивная ячейка?
    bool linear;    ///< Линейная ячейка?
    bool axial;     ///< Осевая симметрия, см описание класса

    struct Cell {
        std::vector<int> rank;   ///< Ранг процесса владельца (< 0 -- ошибка, не используется)
        std::vector<int> index;  ///< Индекс элемента в локальном Storage (< 0 для неактивных, неопределенных элементов, элементов на удаление)
        std::vector<int> next;   ///< Новый индекс в хранилище (в алгоритмах с перестановками)

        std::vector<Vector3d> center;      ///< Барицентр ячейки
        std::vector<double>   volume;      ///< Объем трехмерной / площадь двумерной ячейки
        std::vector<double>   volume_alt;  ///< Объем двумерной осесимметричной ячейки
        std::vector<size_t>   face_beg;
        std::vector<size_t>   face_end;
        std::vector<size_t>   node_beg;
        std::vector<size_t>   node_end;

        void resize(size_t n_cells) {
            rank.resize(n_cells);
            index.resize(n_cells);
            next.resize(n_cells);

            center.resize(n_cells);
            volume.resize(n_cells);
            volume_alt.resize(n_cells);
            face_beg.resize(n_cells);
            face_end.resize(n_cells);
            node_beg.resize(n_cells);
            node_end.resize(n_cells);
        }
    };

    struct Face {
        std::vector<Boundary> boundary;  ///< Тип граничного условия
        std::vector<Adjacent> adjacent;  ///< Составной индекс смежной ячейки
        std::vector<Vector3d> normal;    ///< Внешняя нормаль к грани
        std::vector<Vector3d> center;    ///< Барицентр грани
        std::vector<double>   area;      ///< Площадь грани
        std::vector<double>   area_alt;  ///< "Альтернативная" площадь грани

        /// @brief Список индексов вершин в массиве вершин ячейки
        std::vector<std::array<int, BFace::max_size>> vertices;

        void resize(size_t n_faces) {
            boundary.resize(n_faces);
            adjacent.resize(n_faces);
            normal.resize(n_faces);
            center.resize(n_faces);
            area.resize(n_faces);
            area_alt.resize(n_faces);
            vertices.resize(n_faces);
        }
    };

    Cell cell;
    Face face;
    std::vector<size_t> node;
    std::vector<Vector3d> vertex;

    utils::Storage data;

    const Vector3d& cell_center(size_t ic) const {
        return cell.center[ic];
    }

    double cell_volume(size_t ic) const {
        return cell.volume[ic];
    }

    double cell_size(size_t ic) const {
        return dim < 3 ? std::sqrt(cell.volume[ic]) : std::cbrt(cell.volume[ic]);
    }

    double cell_diameter(size_t ic) const {
        if (adaptive) {
            if (dim == 2) {
                SqQuad vertices = (SqQuad& ) vertex[cell.node_beg[ic]];
                return std::sqrt(std::min(
                        (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                        (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
            }
            else {
                SqCube vertices = (SqCube& ) vertex[cell.node_beg[ic]];
                return std::sqrt(std::min({
                        (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                        (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                        (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()}));
            }
        }
        else {
            throw std::runtime_error("Not implemented");
        }
    }

    double cell_volume(size_t ic, bool axial) const {
        return axial ? cell.volume_alt[ic] : cell.volume[ic];
    }

    inline utils::range<size_t> cell_faces(size_t ic) const {
        return utils::range<size_t>(cell.face_beg[ic], cell.face_end[ic]);
    }

    /// @brief Является ли грань граничной?
    inline bool is_boundary_face(size_t k) const {
        return face.boundary[k] != Boundary::ORDINARY &&
               face.boundary[k] != Boundary::PERIODIC &&
               face.boundary[k] != Boundary::UNDEFINED;
    }

    /// @brief Является ли грань актуальной?
    inline bool is_actual_face(size_t k) const {
        return face.boundary[k] != Boundary::UNDEFINED;
    }

    /// @return 'true', если грань не актуальна
    inline bool is_undefined_face(size_t k) const {
        return face.boundary[k] == Boundary::UNDEFINED;
    }

    /// @brief Установить неопределенную грань
    inline void set_undefined_face(size_t k) {
        face.boundary[k] = Boundary::UNDEFINED;
        face.adjacent[k].rank  = -1;
        face.adjacent[k].index = -1;
        face.adjacent[k].alien = -1;
    }

    inline Vector3d face_symm_point(size_t k, const Vector3d& p) const {
        return p + 2.0 * (face.center[k] - p).dot(face.normal[k]) * face.normal[k];
    }

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    inline double face_area(size_t k) const { return face.area[k]; }

    /// @brief Площадь/длина обычной грани или грани осесимметричной ячейки
    inline double face_area(size_t k, bool axial) const {
        return axial ? face.area_alt[k] : face.area[k];
    }

    // Преобразовать из классической сетки
    SoaMesh(EuMesh &mesh) {
        AmrStorage& locals = mesh.locals();

        if (locals.empty()) {
            throw std::runtime_error("SoaMesh: locals is empty");
        }

        dim = locals[0].dim;
        adaptive = locals[0].adaptive;
        linear = locals[0].linear;
        axial = locals[0].axial;

        int faces_per_cell = dim < 3 ? 8 : 24;
        int nodes_per_cell = dim < 3 ? 9 : 27;

        cell.resize(locals.size());
        face.resize(locals.size() * faces_per_cell);
        node.resize(locals.size() * nodes_per_cell);
        vertex.resize(locals.size() * nodes_per_cell);

        for (size_t ic = 0; ic < locals.size(); ++ic) {
            cell.rank[ic] = locals[ic].rank;
            cell.index[ic] = locals[ic].index;
            cell.next[ic] = locals[ic].next;

            cell.center[ic] = locals[ic].center;
            cell.volume[ic] = locals[ic].volume;
            cell.volume_alt[ic] = locals[ic].volume_alt;
            cell.face_beg[ic] = ic * faces_per_cell;
            cell.face_end[ic] = (ic + 1) * faces_per_cell;
            cell.node_beg[ic] = ic * nodes_per_cell;
            cell.node_end[ic] = (ic + 1) * nodes_per_cell;

            for (size_t jf = 0; jf < faces_per_cell; ++jf) {
                size_t f_idx = ic * faces_per_cell + jf;

                face.boundary[f_idx] = locals[ic].faces[jf].boundary;
                face.adjacent[f_idx] = locals[ic].faces[jf].adjacent;
                face.normal  [f_idx] = locals[ic].faces[jf].normal;
                face.center  [f_idx] = locals[ic].faces[jf].center;
                face.area    [f_idx] = locals[ic].faces[jf].area;
                face.area_alt[f_idx] = locals[ic].faces[jf].area_alt;
                face.vertices[f_idx] = locals[ic].faces[jf].vertices;
            }

            for (int jn = 0; jn < nodes_per_cell; ++jn) {
                size_t v_idx = ic * nodes_per_cell + jn;
                node[v_idx] = v_idx;
                vertex[v_idx] = locals[ic].vertices[jn];
            }
        }
    }

    template<int n_tasks_per_thread = 10, class Func, class... Args>
    void for_each(Func &&func, Args&&... args) {
        utils::range<size_t> range(cell.rank.size());

        threads::for_each<n_tasks_per_thread>(range.begin(), range.end(),
                std::forward<Func>(func), std::forward<Args>(args)...);
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, EuCell&>>
    auto min(Func &&func, const Value &init)
    -> typename std::enable_if<!std::is_void<Value>::value, Value>::type {
        utils::range<size_t> range(cell.rank.size());
        return threads::min<n_tasks_per_thread>(range.begin(), range.end(), init, std::forward<Func>(func));
    }

    /// @brief Параллельно по тредам посчитать минимум
    template<int n_tasks_per_thread = 10, class Func,
            typename Value = std::invoke_result_t<Func, EuCell&>>
    auto min(Func &&func)
    -> typename std::enable_if<std::is_arithmetic<Value>::value, Value>::type {
        utils::range<size_t> range(cell.rank.size());
        return threads::min<n_tasks_per_thread>(range.begin(), range.end(), std::forward<Func>(func));
    }
};

}