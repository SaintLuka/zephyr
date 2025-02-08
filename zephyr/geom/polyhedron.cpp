#include <iostream>
#include <map>

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
    vs = vertices;
    fs = face_indices;

    fcs.resize(face_indices.size());
    fns.resize(face_indices.size());

    // Определим центр, нужна опорная точка внутри
    m_center = Vector3d::Zero();
    for (auto& v: vs) {
        m_center += v;
    }
    m_center /= vs.size();

    // Определим характеристики граней,
    // исправим порядок, если необходимо
    for (int i = 0; i < int(fs.size()); ++i) {
        int nv = fs[i].size();

        for (int j = 0; j < nv; ++j) {
            assert(fs[i][j] < vs.size() && "Wrong polyhedron");
        }

        Vector3d v0 = vs[fs[i][0]];
        Vector3d v1 = vs[fs[i][1]];
        Vector3d v2 = vs[fs[i][2]];

        if (nv == 3) {
            // Центр треугольной грани
            fcs[i] = (v0 + v1 + v2) / 3.0;

            // Базовое направление внешней нормали
            Vector3d nb = fcs[i] - m_center;

            // Внешняя нормаль
            fns[i] = (v1 - v0).cross(v2 - v0).normalized();

            // Меняем обход
            if (nb.dot(fns[i]) < 0.0) {
                fns[i] *= 1.0;
                std::swap(fs[i][1], fs[i][2]);
            }
        }
        else if (nv == 4) {
            Vector3d v3 = vs[fs[i][3]];

            // Центр квадратной грани
            fcs[i] = 0.25 * (v0 + v1 + v2 + v3);

            // Базовое направление внешней нормали
            Vector3d nb = fcs[i] - m_center;

            // Проверим на перекрут,
            bool ok0 = (v1 - v0).cross(v3 - v0).dot(nb) > 0.0;
            bool ok1 = (v2 - v1).cross(v0 - v1).dot(nb) > 0.0;
            bool ok2 = (v3 - v2).cross(v1 - v2).dot(nb) > 0.0;
            bool ok3 = (v0 - v3).cross(v2 - v3).dot(nb) > 0.0;

            // Тут всего 5 вариантов, сделаю полный перебор
            if (ok0 && ok1 && ok2) {
                // Правильный обход
                // Ничего не делаем
            }
            else if (!(ok0 || ok1 || ok2)) {
                // Неправильный обход
                std::swap(fs[i][1], fs[i][3]);
                std::swap(v1, v3);
            }
            else {
                // Есть хотя бы один !OK
                if (ok0) {
                    if (ok1) {
                        // ok0 & ok1 & !ok2
                        std::swap(fs[i][2], fs[i][3]);
                        std::swap(v2, v3);
                    } else {
                        // ok0 & !ok1 & !ok2
                        std::swap(fs[i][1], fs[i][2]);
                        std::swap(v1, v2);
                    }
                } else {
                    if (ok1) {
                        // !ok0 & ok1 & ok2
                        std::swap(fs[i][0], fs[i][3]);
                        std::swap(v0, v3);
                    } else {
                        // !ok0 & !ok1 & ok2
                        std::swap(fs[i][0], fs[i][1]);
                        std::swap(v0, v1);
                    }
                }
            }

            fns[i] = (v2 - v0).cross(v3 - v1).normalized();
        }
        else {
            // Полигональная грань, сложная сортировка
            throw std::runtime_error("Polyhedron complex face");
        }
    }
}

Box Polyhedron::bbox() const {
    if (empty()) {
        return {Vector3d::Zero(), Vector3d::Zero()};
    }

    const double max = +std::numeric_limits<double>::infinity();
    const double min = -std::numeric_limits<double>::infinity();

    Vector3d vmin = {max, max, min};
    Vector3d vmax = {min, min, max};
    for (int i = 0; i < n_verts(); ++i) {
        if (vs[i].x() < vmin.x()) {
            vmin.x() = vs[i].x();
        } else if (vs[i].x() > vmax.x()) {
            vmax.x() = vs[i].x();
        }

        if (vs[i].y() < vmin.y()) {
            vmin.y() = vs[i].y();
        } else if (vs[i].y() > vmax.y()) {
            vmax.y() = vs[i].y();
        }

        if (vs[i].z() < vmin.z()) {
            vmin.z() = vs[i].z();
        } else if (vs[i].z() > vmax.z()) {
            vmax.z() = vs[i].z();
        }
    }

    return {vmin, vmax};
}

Vector3d Polyhedron::center() const {
    return m_center;
}

double Polyhedron::face_area(int idx) const {
    const Vector3d &v0 = vs[fs[idx][0]];
    const Vector3d &v1 = vs[fs[idx][1]];
    const Vector3d &v2 = vs[fs[idx][2]];

    if (fs[idx].size() < 4) {
        return 0.5 * (v2 - v0).cross(v1 - v0).norm();
    } else {
        const Vector3d &v3 = vs[fs[idx][3]];
        return 0.5 * (v2 - v0).cross(v3 - v1).norm();
    }
}

Vector3d Polyhedron::face_center(int idx) const {
    return fcs[idx];
}

Vector3d Polyhedron::face_normal(int idx) const {
    return fns[idx];
}

bool Polyhedron::inside(const Vector3d& p) const {
    for (int i = 0; i < n_faces(); ++i) {
        if ((p - fcs[i]).dot(fns[i]) > 0.0) {
            return false;
        }
    }
    return true;
}

double Polyhedron::volume() const {
    const double coeff = 1.0 / 6.0;

    double sum = 0.0;
    for (int i = 0; i < n_faces(); ++i) {
        Vector3d ch = fcs[i] - m_center;

        int nv = fs[i].size();
        if (nv == 3) {
            const Vector3d &v0 = vs[fs[i][0]];
            const Vector3d &v1 = vs[fs[i][1]];
            const Vector3d &v2 = vs[fs[i][2]];

            sum += (v1 - v0).cross(v2 - v0).dot(ch);
        } else if (nv == 4) {
            const Vector3d &v0 = vs[fs[i][0]];
            const Vector3d &v1 = vs[fs[i][1]];
            const Vector3d &v2 = vs[fs[i][2]];
            const Vector3d &v3 = vs[fs[i][3]];

            sum += (v2 - v0).cross(v3 - v1).dot(ch);
        } else {
            throw std::runtime_error("Polyhedron volume: complex face");
        }
    }
    return coeff * sum;
}

Vector3d Polyhedron::centroid(double vol) const {
    return m_center;
}


double Polyhedron::clip_volume(
        const std::function<bool(const Vector3d& p)>& inside,
        int n_points) const {
    return NAN;
}


Polyhedron Polyhedron::clip(const Vector3d& p, const Vector3d& n) const {
    // Напишу какую-нибудь простую, но не очень быструю версию

    // Обходим грани.
    // Кажду грань или убираем, или добавляем целиком, или режем на две
    // части и добавляем часть.

    // Собираем полигон для среза.
    // При разрезании получается пара вершин, их в отдельный список.
    // Собираем все получившиеся вершины в один массив, удаляем дубликаты
    // и сортируем?

    // Массив из полигонов.

    obj::plane plane{.p=p, .n=n};

    // Индикаторы для вершин
    // -1: под/внутри. 0: на плоскости. +1: снаружи
    std::vector<int> vs_pos(n_verts());
    for (int i = 0; i < n_verts(); ++i) {
        vs_pos[i] = plane.position(vs[i]);
    }

    // Пересечения рёбер
    std::map<std::array<int, 2>, Vector3d> edges;

    // Срез плоскостью, поддерживаем в отсортированном виде
    // (отсортированном против часовой вокруг нормали плоскости)
    std::vector<Vector3d> slice;

    std::vector<std::vector<Vector3d>> faces;

    for (int i  = 0; i < n_faces(); ++i) {
        int nv = fs[i].size();

        // все вершины снаружи (+1 или 0)
        // все вершины внутри (-1 или 0)
        // есть -1 и +1 одновременно
        int c_minus = 0;
        int c_plus  = 0;
        for (auto j: fs[i]) {
            if (vs_pos[j] < 0) {
                ++c_minus;
            }
            else if (vs_pos[j] > 0) {
                ++c_plus;
            }
        }

        if (c_plus == 0) {
            // Грань целиком внутри
            std::vector<Vector3d> poly(nv);
            for (int j = 0; j < nv; ++j) {
                poly[j] = vs[fs[i][j]];
            }
            faces.emplace_back(std::move(poly));
            continue;
        }
        else if (c_minus == 0) {
            // Грань целиком снаружи, ничего не делаем
            continue;
        }
        else {
            // Грань пересекается, добавить нужную часть,
            // новые точки в массив slice

            std::vector<Vector3d> poly;
            for (int j = 0; j < nv; ++j) {
                int j2 = (j + 1) % nv;

                if (vs_pos[j] <= 0) {
                    poly.push_back(vs[j]);
                }
                else if (vs_pos[j] * vs_pos[j2] < 0) {
                    obj::segment seg {vs[j], vs[j2]};

                    Vector3d sec = intersection2D::find_fast(plane, seg);
                    poly.push_back(sec);
                    slice.push_back(sec);
                }
            }
            faces.emplace_back(std::move(poly));
            continue;
        }
    }


    std::vector<Vector3d> verts;
    std::vector<std::vector<int>> f_inds(faces.size());

    for (int i = 0; i < faces.size(); ++i) {
        auto& face = faces[i];
        std::cout << "  face.size: " << face.size() << "\n";

        for (auto& v: face) {
            f_inds[i].push_back(verts.size());
            verts.push_back(v);
        }
    }

    return Polyhedron(verts, f_inds);
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
    os << "Not yet\n";
    return os;
}

} // namespace zephyr::geom
