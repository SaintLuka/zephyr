#include <iostream>
#include <map>
#include <unordered_map>

#include <zephyr/geom/intersection.h>
#include <zephyr/geom/polyhedron.h>


namespace zephyr::geom {

namespace {
// Число вершин для стандартных многогранников
int default_face_size(CellType ctype) {
    switch (ctype) {
        case CellType::TETRA:      return 4;
        case CellType::PYRAMID:    return 5;
        case CellType::WEDGE:      return 6;
        case CellType::HEXAHEDRON: return 8;
        default:
            throw std::runtime_error("Has no default face size");
    }
}

// Нумерация индексов вершин для стандартных многогранников
std::vector<std::vector<int>> default_face_indices(CellType ctype) {
    switch (ctype) {
        case CellType::TETRA:
            return {{0, 2, 1},
                    {0, 1, 3},
                    {1, 2, 3},
                    {0, 3, 2}};

        case CellType::PYRAMID:
            return {{3, 2, 1, 0},
                    {0, 1, 4},
                    {1, 2, 4},
                    {2, 3, 4},
                    {0, 4, 3}};

        case CellType::WEDGE:
            return {{0, 3, 4, 1},
                    {1, 4, 5, 2},
                    {0, 2, 5, 3},
                    {0, 1, 2},
                    {3, 5, 4}};

        case CellType::HEXAHEDRON:
            return {{0, 3, 2, 1},
                    {4, 5, 6, 7},
                    {0, 3, 4, 7},
                    {1, 2, 6, 5},
                    {0, 1, 5, 4},
                    {2, 3, 7, 6}};

        default:
            throw std::runtime_error("Has no default face indices");
    }
}
}

// ============================================================================
//                         Вспомогательные структуры
// ============================================================================

// Пара (вершина, индекс)
struct point_t {
    Vector3d v;
    int index = -1;
};

// Пара индексов смежных вершин, можно поставить два одинаковых индекса,
// тогда это считается индексом вершины
struct edge_t {
    explicit edge_t(int v1) : m_v1(v1), m_v2(v1) {}

    edge_t(int v1, int v2) {
        if (v1 < v2) {
            m_v1 = v1;
            m_v2 = v2;
        } else {
            m_v1 = v2;
            m_v2 = v1;
        }
    }

    inline int v1() const { return m_v1; }

    inline int v2() const { return m_v2; }

    bool operator==(const edge_t &other) const {
        return m_v1 == other.m_v1 && m_v2 == other.m_v2;
    }

    bool operator<(const edge_t &rhs) const {
        return m_v1 < rhs.m_v1 || (m_v1 == rhs.m_v1 && m_v2 < rhs.m_v2);
    }

protected:
    // Нужно поддерживать упорядоченными, поэтому прямой доступ ограничен
    int m_v1, m_v2;
};

// Компаратор.
// Сортировщик для точек на плоскости, позволяет упорядочить
// точки в плоскости против часовой стрелки вокруг нормали n.
struct SortRule {
    SortRule(const Vector3d& c, const Vector3d& n, const Vector3d& v0)
            : c(c), n(n), v0(v0) {
        e_x = (v0 - c).normalized();
        e_y = n.cross(e_x);
    }

    // Оператор сравнение двух точек, сортирует
    bool operator()(const Vector3d& v1, const Vector3d& v2) {
        double x1 = (v1 - c).dot(e_x);
        double y1 = (v1 - c).dot(e_y);
        double x2 = (v2 - c).dot(e_x);
        double y2 = (v2 - c).dot(e_y);

        double phi1 = std::atan2(y1, x1);
        double phi2 = std::atan2(y2, x2);
        return phi1 < phi2;
    }

private:
    Vector3d c;   // Точка плоскости
    Vector3d n;   // Нормаль к плоскости
    Vector3d v0;  // Одна из точек (обычно первая)

    // Базис в плоскости
    Vector3d e_x, e_y;
};

// ============================================================================
//                         Небольшие функции функции
// ============================================================================

// Центр грани общего вида
template <int nv = -1>
Vector3d get_center(const std::vector<Vector3d>& vs,
                    const std::vector<int>& inds) {
    Vector3d c = Vector3d::Zero();
    for (int i: inds) {
        c += vs[i];
    }
    return c / inds.size();
}

// Центр треугольной грани
template <>
inline Vector3d get_center<3>(const std::vector<Vector3d>& vs,
                              const std::vector<int>& inds) {
    const double coeff = 1.0 / 3.0;
    return coeff * (vs[inds[0]] + vs[inds[1]] + vs[inds[2]]);
}

// Центр четырехугольной грани
template <>
inline Vector3d get_center<4>(const std::vector<Vector3d>& vs,
                              const std::vector<int>& inds) {
    return 0.25 * (vs[inds[0]] + vs[inds[1]] + vs[inds[2]] + vs[inds[3]]);
}

// Сортировка вершин
// face_inds -- изменяется
template <int nv = -1>
void sort_indices(const std::vector<Vector3d>& vs,
                  std::vector<int>& face_inds,
                  const Vector3d& face_c,
                  const Vector3d& cell_c) {
    // Примерная внешняя нормаль
    Vector3d nb = (face_c - cell_c).normalized();

    SortRule comp(face_c, nb, vs[face_inds[0]]);

    std::sort(face_inds.begin(), face_inds.end(),
              [&comp, &vs](int i, int j) -> bool {
                  return comp(vs[i], vs[j]);
              });
}

template <>
inline void sort_indices<3>(const std::vector<Vector3d>& vs,
                            std::vector<int>& face_inds,
                            const Vector3d& face_c,
                            const Vector3d& cell_c) {
    const Vector3d &v0 = vs[face_inds[0]];
    const Vector3d &v1 = vs[face_inds[0]];
    const Vector3d &v2 = vs[face_inds[0]];

    if ((v1 - v0).cross(v2 - v0).dot(face_c - cell_c) < 0.0) {
        std::swap(face_inds[1], face_inds[2]);
    }
}

template <>
inline void sort_indices<4>(const std::vector<Vector3d>& vs,
                            std::vector<int>& face_inds,
                            const Vector3d& face_c,
                            const Vector3d& cell_c) {
    const Vector3d &v0 = vs[face_inds[0]];
    const Vector3d &v1 = vs[face_inds[1]];
    const Vector3d &v2 = vs[face_inds[2]];
    const Vector3d &v3 = vs[face_inds[3]];

    // Базовое направление внешней нормали
    Vector3d nb = face_c - cell_c;

    // Проверим на перекрут
    bool ok0 = (v1 - v0).cross(v3 - v0).dot(nb) > 0.0;
    bool ok1 = (v2 - v1).cross(v0 - v1).dot(nb) > 0.0;
    bool ok2 = (v3 - v2).cross(v1 - v2).dot(nb) > 0.0;
    // bool ok3 = (v0 - v3).cross(v2 - v3).dot(nb) > 0.0;

    // Тут всего 6 вариантов, сделаю полный перебор
    // 1. Все OK.
    // 2. Все не ОК.
    // 3-6. Пара соседних ОК и пара не OK.

    // Правильный обход, ничего не делаем
    if (ok0 && ok1 && ok2) {
        return;
    }

    // Полностью неправильный, меняем противоположные
    if (!(ok0 || ok1 || ok2)) {
        // Неправильный обход
        std::swap(face_inds[1], face_inds[3]);
        return;
    }

    // Есть хотя бы один !OK
    if (ok0) {
        if (ok1) {
            // ok0 & ok1 & !ok2
            std::swap(face_inds[2], face_inds[3]);
        } else {
            // ok0 & !ok1 & !ok2
            std::swap(face_inds[1], face_inds[2]);
        }
    } else {
        if (ok1) {
            // !ok0 & ok1 & ok2
            std::swap(face_inds[0], face_inds[3]);
        } else {
            // !ok0 & !ok1 & ok2
            std::swap(face_inds[0], face_inds[1]);
        }
    }
}

// Ориентированная площадь, предполагаем, что вершины
// отсортированы
// c -- Центр, среднее вершин
template <int n_verts = -1>
Vector3d get_surface(const std::vector<Vector3d>& vs,
                     const std::vector<int>& inds, 
                     const Vector3d& face_c) {
    Vector3d S = Vector3d::Zero();
    
    int nv = inds.size();    
    for (int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        Vector3d v1 = vs[inds[i]] - face_c;
        Vector3d v2 = vs[inds[j]] - face_c;

        S += v1.cross(v2);
    }
    return 0.5 * S;
}

// Ориентированная площадка треугольной грани
template <>
Vector3d get_surface<3>(const std::vector<Vector3d>& vs,
                        const std::vector<int>& inds, 
                        const Vector3d& face_c) {
    const Vector3d &v0 = vs[inds[0]];
    return 0.5 * (vs[inds[1]] - v0).cross(vs[inds[2]] - v0);    
}

// Ориентированная площадка четырехугольной грани
template <>
Vector3d get_surface<4>(const std::vector<Vector3d>& vs,
                        const std::vector<int>& inds, 
                        const Vector3d& face_c) {
    return 0.5 * (vs[inds[2]] - vs[inds[0]]).cross(vs[inds[3]] - vs[inds[1]]);
}

// ============================================================================
//                         Непосредственно функции
// ============================================================================

Polyhedron::Polyhedron(const std::vector<Vector3d>& vertices,
           const std::vector<std::vector<int>>& face_indices) {
    build(vertices, face_indices);
}

Polyhedron::Polyhedron(CellType ctype, const std::vector<Vector3d>& vertices) {
    auto face_indices = default_face_indices(ctype);
    if (vertices.size() != default_face_size(ctype)) {
        throw std::runtime_error("Bad vertices size for cell type");
    }
    build(vertices, face_indices);
}

void Polyhedron::build(const std::vector<Vector3d>& vertices,
                       const std::vector<std::vector<int>>& face_indices) {
    verts = vertices;
    faces = face_indices;

    int n_faces = faces.size();

    faces_c.resize(n_faces);
    faces_s.resize(n_faces);

    // Определим центр, нужна опорная точка внутри многогранника
    m_center = Vector3d::Zero();
    for (auto& v: verts) {
        m_center += v;
    }
    m_center /= verts.size();

    // Определим характеристики граней, исправим порядок, если необходимо
    for (int i = 0; i < n_faces; ++i) {
        int nv = faces[i].size();

        assert(nv > 2 && "Two vertex polyhedron");

        for (int j = 0; j < nv; ++j) {
            assert(faces[i][j] < verts.size() && "Wrong polyhedron");
        }

        if (nv == 3) {
            // Центр треугольной грани
            faces_c[i] = get_center<3>(verts, faces[i]);

            // Сортировка индексов на грани
            sort_indices<3>(verts, faces[i], faces_c[i], m_center);

            // Внешняя нормаль
            faces_s[i] = get_surface<3>(verts, faces[i], faces_c[i]);
        } else if (nv == 4) {
            // Центр треугольной грани
            faces_c[i] = get_center<4>(verts, faces[i]);

            // Сортировка индексов на грани
            sort_indices<4>(verts, faces[i], faces_c[i], m_center);

            // Внешняя нормаль
            faces_s[i] = get_surface<4>(verts, faces[i], faces_c[i]);
        } else {
            // Сложная полигональная грань с числом вершин > 4

            // Центр квадратной грани
            faces_c[i] = get_center(verts, faces[i]);

            // Сортировка индексов на грани
            sort_indices(verts, faces[i], faces_c[i], m_center);

            // Внешняя нормаль и площадь
            faces_s[i] = get_surface(verts, faces[i], faces_c[i]);

        }
    }
}

Box Polyhedron::bbox() const {
    if (empty()) {
        return Box::Zero();
    }

    Box box = Box::Empty(3);
    for (auto& v: verts) {
        box.capture(v);
    }
    return box;
}

Vector3d Polyhedron::center() const {
    return m_center;
}

double Polyhedron::excircle_radius() const {
    double R = 0.0;
    for (auto& v: verts) {
        double dist = (v - m_center).norm();
        if (dist > R) {
            R = dist;
        }
    }
    return R;
}

double Polyhedron::face_area(int idx) const {
    switch (faces[idx].size()) {
        case 3:
            return get_surface<3>(verts, faces[idx], faces_c[idx]).norm();
        case 4:
            return get_surface<4>(verts, faces[idx], faces_c[idx]).norm();
        default:
            return get_surface(verts, faces[idx], faces_c[idx]).norm();
    }
}

Vector3d Polyhedron::face_center(int idx) const {
    return faces_c[idx];
}

Vector3d Polyhedron::face_normal(int idx) const {
    return faces_s[idx].normalized();
}

bool Polyhedron::inside(const Vector3d& p) const {
    for (int i = 0; i < n_faces(); ++i) {
        if ((p - faces_c[i]).dot(faces_s[i]) > 0.0) {
            return false;
        }
    }
    return true;
}

double Polyhedron::volume() const {
    const double coeff = 1.0 / 3.0;

    double sum = 0.0;
    for (int i = 0; i < n_faces(); ++i) {
        sum += faces_s[i].dot(faces_c[i] - m_center);
    }
    return coeff * sum;
}

Vector3d Polyhedron::centroid(double vol) const {
    return m_center;
}

void Polyhedron::replace_face(int idx, int v1, int v2, int v3, int v4) {
    faces[idx] = {v1, v2, v3, v4};
    faces_c[idx] = get_center<4>(verts, faces[idx]);
    faces_s[idx] = get_surface<4>(verts, faces[idx], faces_c[idx]);
}

void Polyhedron::add_face(int v1, int v2, int v3) {
    faces.push_back({v1, v2, v3});
    faces_c.push_back(get_center<3>(verts, faces.back()));
    faces_s.push_back(get_surface<3>(verts, faces.back(), faces_c.back()));
}

void Polyhedron::add_face(int v1, int v2, int v3, int v4) {
    faces.push_back({v1, v2, v3, v4});
    faces_c.push_back(get_center<4>(verts, faces.back()));
    faces_s.push_back(get_surface<4>(verts, faces.back(), faces_c.back()));
}

void Polyhedron::canonic() {
    for (int f_idx = 0; f_idx < n_faces(); ++f_idx) {
        int nv = faces[f_idx].size();
        if (nv <= 4) {
            continue;
        }

        std::vector<int> old = faces[f_idx];
        switch (nv) {
            case 5:
                add_face(old[3], old[4], old[0]);
                break;
            case 6:
                add_face(old[3], old[4], old[5], old[0]);
                break;
            case 7:
                add_face(old[3], old[4], old[5], old[6]);
                add_face(old[0], old[3], old[6]);
                break;
            case 8:
                add_face(old[3], old[4], old[5], old[7]);
                add_face(old[0], old[3], old[4], old[7]);
                break;
            default:
                throw std::runtime_error("Sorry, I'm too lazy");
        }
        replace_face(f_idx, old[0], old[1], old[2], old[3]);
    }
}

double Polyhedron::clip_volume(
        const std::function<bool(const Vector3d& p)>& inside,
        int n_points) const {
    return NAN;
}

// TODO: Описание алгоритма
Polyhedron Polyhedron::clip(const Vector3d& p, const Vector3d& n) const {
    // Плоскость
    obj::plane plane{.p=p, .n=n};

    // Пересечения рёбер с индексами вершин (i, j).
    // Или сама точка, тогда индекс (i, i).
    std::map<edge_t, point_t> edges;

    // Считаем вершины строго снаружи или внутри
    int count_inside  = 0;
    int count_outside = 0;

    // Положения / индикаторы вершин.
    // -1: под/внутри. 0: на плоскости. +1: снаружи
    std::vector<int> vs_pos(n_verts());
    for (int i = 0; i < n_verts(); ++i) {
        vs_pos[i] = plane.position(verts[i]);

        if (vs_pos[i] > 0) {
            ++count_outside;
        }
        else {
            if (vs_pos[i] < 0) {
                ++count_inside;
            }

            // Точки внутри и на пересечении сразу добавляем в массив,
            // эти вершины гарантированно войдут в итоговое отсечение
            edges[edge_t(i)] = {.v=verts[i], .index=-1};
        }
    }

    // Многогранник снаружи, нет пересечений
    if (count_inside == 0) {
        return Polyhedron::Empty();
    }

    // Многогранник целиком внутри, нет пересечений
    if (count_outside == 0) {
        return *this;
    }

    std::vector<edge_t> slice;
    std::vector<std::vector<edge_t>> parts;

    for (int i  = 0; i < n_faces(); ++i) {
        // Считаем вершины строго снаружи или внутри
        count_inside  = 0;
        count_outside = 0;
        for (auto j: faces[i]) {
            if (vs_pos[j] < 0) {
                ++count_inside;
            } else if (vs_pos[j] > 0) {
                ++count_outside;
            }
        }

        // Грань целиком снаружи, ничего не делаем
        if (count_inside == 0) {
            continue;
        }

        int nv = faces[i].size();

        if (count_outside == 0) {
            // Грань целиком внутри
            std::vector<edge_t> part;
            part.reserve(nv);
            for (int j = 0; j < nv; ++j) {
                part.emplace_back(edge_t(faces[i][j]));
            }
            parts.emplace_back(part);
        }
        else {
            // Грань пересекается, добавить нужную часть,
            // новые точки в массив slice

            std::vector<edge_t> part;
            part.reserve(nv + 1);

            for (int j = 0; j < nv; ++j) {
                int v_idx1 = faces[i][j];
                int v_idx2 = faces[i][(j + 1) % nv];

                if (vs_pos[v_idx1] <= 0) {
                    part.emplace_back(edge_t(v_idx1));
                }
                if (vs_pos[v_idx1] * vs_pos[v_idx2] < 0) {
                    // Всегда один порядок
                    edge_t edge(v_idx1, v_idx2);

                    // Пересечение уже найдено
                    if (edges.count(edge) > 0) {
                        part.push_back(edge);
                        continue;
                    }

                    // Найдем пересечение
                    obj::segment seg{verts[edge.v1()], verts[edge.v2()]};

                    Vector3d new_v = intersection2D::find_fast(plane, seg);
                    edges[edge] = {.v=new_v, .index=-1};

                    part.push_back(edge);
                    slice.push_back(edge);
                }
            }
            parts.emplace_back(std::move(part));
        }
    }

    // Вершины в Slice не сортируем, потому что build
    // сортирует вершины.

    // Помещаем всё в выходные массивы
    std::vector<Vector3d> out_verts(edges.size());

    int counter = 0;
    for (auto& edge: edges) {
        edge.second.index = counter;
        out_verts[counter] = edge.second.v;
        ++counter;
    }

    std::vector<std::vector<int>> out_faces;
    out_faces.reserve(parts.size() + 1); // parts + slice

    for (auto& face: parts) {
        std::vector<int> ids;
        ids.reserve(face.size());

        for (auto edge: face) {
            ids.push_back(edges[edge].index);
        }
        out_faces.push_back(ids);
    }

    // То же самое для slice
    if (slice.size() < 9) {
        // Записать целиком
        std::vector<int> ids;
        ids.reserve(slice.size());

        for (auto edge: slice) {
            ids.push_back(edges[edge].index);
        }
        out_faces.push_back(ids);
    }
    else if (slice.size() < 18) {
        // Черт, придется таки отсортировать вершины в Slice
        Vector3d slice_c = Vector3d::Zero();
        for (auto e: slice) {
            slice_c += edges[e].v;
        }
        slice_c /= slice.size();

        SortRule comp(slice_c, n, edges[slice[0]].v);
        std::sort(slice.begin(), slice.end(),
                [&comp, &edges](const edge_t& e1, const edge_t& e2) -> bool {
                    return comp(edges[e1].v, edges[e2].v);
                });

        // разбить на две части
        std::vector<int> ids1, ids2;
        ids1.reserve(slice.size() / 2 + 1);
        ids2.reserve(slice.size() / 2 + 1);

        for (int i = 0; i <= slice.size() / 2; ++i) {
            ids1.push_back(edges[slice[i]].index);
        }
        for (int i = slice.size() / 2; i < slice.size(); ++i) {
            ids2.push_back(edges[slice[i]].index);
        }
        ids2.push_back(ids1.front());

        out_faces.push_back(ids1);
        out_faces.push_back(ids2);

    } else {
        throw std::runtime_error("I'm too lazy");
    }

    return Polyhedron(out_verts, out_faces);
}

double Polyhedron::clip_volume(const Vector3d& p, const Vector3d& n) const {
    return NAN;
}


Vector3d Polyhedron::find_section(const Vector3d& n, double alpha) const {
    return {NAN, NAN, NAN};
}

double Polyhedron::volume_fraction(
        const std::function<bool(const Vector3d&)>& inside,
        int n_points) const {
    return NAN;
}


std::ostream& operator<<(std::ostream& os, const Polyhedron& poly) {
    os << "Polyhedron operator<<(): Not yet\n";
    return os;
}

// ============================================================================
//                          Заготовки многогранников
// ============================================================================

Polyhedron Polyhedron::Cuboid(double a, double b, double c) {
    std::vector<Vector3d> vs = {
            {-0.5 * a, -0.5 * b, -0.5 * c},
            {+0.5 * a, -0.5 * b, -0.5 * c},
            {+0.5 * a, +0.5 * b, -0.5 * c},
            {-0.5 * a, +0.5 * b, -0.5 * c},
            {-0.5 * a, -0.5 * b, +0.5 * c},
            {+0.5 * a, -0.5 * b, +0.5 * c},
            {+0.5 * a, +0.5 * b, +0.5 * c},
            {-0.5 * a, +0.5 * b, +0.5 * c},
    };

    return Polyhedron(CellType::HEXAHEDRON, vs);
}

Polyhedron Polyhedron::Pyramid() {
    double z = 1.0 / 3.0;
    double a = std::sqrt(1.0 - z * z) / std::sqrt(2.0);
    std::vector<Vector3d> vs = {
            {-a,  -a,  0.0},
            {+a,  -a,  0.0},
            {+a,  +a,  0.0},
            {-a,  +a,  0.0},
            {0.0, 0.0, 1.0},
    };
    return Polyhedron(CellType::PYRAMID, vs);
}

Polyhedron Polyhedron::Wedge() {
    double x = 0.6;
    double r = std::sqrt(1.0 - x * x);

    double a = r / 2.0;
    double b = std::sqrt(3.0) * r / 2.0;

    std::vector<Vector3d> vs = {
            {-x, -b, -a},
            {-x, +b, -a},
            {-x, .0, +r},
            {+x, -b, -a},
            {+x, +b, -a},
            {+x, .0, +r},
    };
    return Polyhedron(CellType::WEDGE, vs);
}

Polyhedron Polyhedron::Tetrahedron() {
    std::vector<Vector3d> vs = {
            {+std::sqrt(8.0) / 3.0, 0.0,                   -1.0 / 3.0},
            {-std::sqrt(2.0) / 3.0, +std::sqrt(2.0 / 3.0), -1.0 / 3.0},
            {-std::sqrt(2.0) / 3.0, -std::sqrt(2.0 / 3.0), -1.0 / 3.0},
            {0.0,                   0.0,                   1.0},
    };
    return Polyhedron(CellType::TETRA, vs);
}

Polyhedron Polyhedron::Octahedron() {
    double a = 1.0 / std::sqrt(2.0);
    std::vector<Vector3d> vs = {
            {-a,  -a,  0.0},
            {+a,  -a,  0.0},
            {+a,  +a,  0.0},
            {-a,  +a,  0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, -1.0},
    };

    std::vector<std::vector<int>> faces = {
            {0, 1, 4},
            {1, 2, 4},
            {2, 3, 4},
            {3, 0, 4},
            {1, 0, 5},
            {2, 1, 5},
            {3, 2, 5},
            {0, 3, 5}};

    return Polyhedron(vs, faces);
}

Polyhedron Polyhedron::Dodecahedron() {
    double len = 1.0 / std::sqrt(3.0);
    double phi = len * 0.5 * (std::sqrt(5.0) + 1.0);
    double inv = len * 0.5 * (std::sqrt(5.0) - 1.0);

    std::vector<Vector3d> vs = {
            {-len, -len, -len},
            {+len, -len, -len},
            {-len, +len, -len},
            {+len, +len, -len},
            {-len, -len, +len},
            {+len, -len, +len},
            {-len, +len, +len},
            {+len, +len, +len},

            {0.0,  -phi, -inv},
            {0.0,  +phi, -inv},
            {0.0,  -phi, +inv},
            {0.0,  +phi, +inv},

            {-inv, 0.0,  -phi},
            {+inv, 0.0,  -phi},
            {-inv, 0.0,  +phi},
            {+inv, 0.0,  +phi},

            {-phi, -inv, 0.0},
            {+phi, -inv, 0.0},
            {-phi, +inv, 0.0},
            {+phi, +inv, 0.0}
    };

    std::vector<std::vector<int>> faces = {
            {2, 3, 9,  12, 13},
            {0, 1, 8,  12, 13},
            {1, 3, 13, 17, 19},
            {0, 2, 12, 16, 18},
            {0, 4, 8,  10, 16},
            {4, 6, 14, 16, 18},
            {6, 7, 11, 14, 15},
            {3, 7, 9,  11, 19},
            {5, 7, 15, 17, 19},
            {4, 5, 10, 14, 15},
            {1, 5, 8,  10, 17},
            {2, 6, 9,  11, 18}
    };

    return Polyhedron(vs, faces);
}

Polyhedron Polyhedron:: Icosahedron() {
    double chi = 1.0 / std::sqrt(0.5 * (5.0 + std::sqrt(5.0)));
    double phi = chi * (0.5 * (1.0 + std::sqrt(5.0)));

    std::vector<Vector3d> vs{
            {-phi, -chi, 0.0},
            {+phi, -chi, 0.0},
            {-phi, +chi, 0.0},
            {+phi, +chi, 0.0},

            {-chi, 0.0,  -phi},
            {+chi, 0.0,  -phi},
            {-chi, 0.0,  +phi},
            {+chi, 0.0,  +phi},

            {0.0,  -phi, -chi},
            {0.0,  +phi, -chi},
            {0.0,  -phi, +chi},
            {0.0,  +phi, +chi},
    };

    std::vector<std::vector<int>> faces = {
            {0, 2, 4},
            {0, 2, 6},
            {1, 3, 5},
            {1, 3, 7},

            {4, 5, 8},
            {4, 5, 9},
            {6, 7, 10},
            {6, 7, 11},

            {0, 8, 10},
            {1, 8, 10},
            {2, 9, 11},
            {3, 9, 11},

            {3, 5, 9},
            {2, 6, 11},
            {2, 4, 9},
            {3, 7, 11},
            {0, 4, 8},
            {1, 5, 8},
            {0, 6, 10},
            {1, 7, 10},
    };

    return Polyhedron(vs, faces);
}

Polyhedron Polyhedron::TruncatedCube() {
    double A = 0.5;
    double a = 0.2;
    
    std::vector<Vector3d> vs = {
            {-a, -A, -A},
            {-A, -a, -A},
            {-A, -A, -a},
            
            {+a, -A, -A},
            {+A, -a, -A},
            {+A, -A, -a},
            
            {+a, +A, -A},
            {+A, +a, -A},
            {+A, +A, -a},
            
            {-a, +A, -A},
            {-A, +a, -A},
            {-A, +A, -a},
            
            {-a, -A, +A},
            {-A, -a, +A},
            {-A, -A, +a},
            
            {+a, -A, +A},
            {+A, -a, +A},
            {+A, -A, +a},
            
            {+a, +A, +A},
            {+A, +a, +A},
            {+A, +A, +a},
            
            {-a, +A, +A},
            {-A, +a, +A},
            {-A, +A, +a},
    };

    std::vector<std::vector<int>> faces = {
            {0,  2,  1},
            {3,  4,  5},
            {6,  8,  7},
            {9,  10, 11},
            {15, 17, 16},
            {14, 12, 13},
            {20, 18, 19},
            {22, 21, 23},

            {0,  2,  3,  5,  17, 15, 12, 14},
            {5,  4,  7,  8,  20, 19, 16, 17},
            {1,  2,  10, 11, 22, 23, 14, 13},
            {0,  1,  3,  4,  6,  7,  9,  10},
            {6,  8,  9,  11, 18, 20, 21, 23},
            {15, 16, 18, 19, 21, 22, 12, 13}
    };

    return Polyhedron(vs, faces);
}

} // namespace zephyr::geom
