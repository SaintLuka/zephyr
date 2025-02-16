#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/geom.h>

namespace zephyr::geom {

AmrCell::AmrCell(const Quad& quad)
    : Element(0, 0), dim(2),
    adaptive(true), linear(true),
    vertices(quad), faces(CellType::AMR2D),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double area = quad.area();

    size = std::sqrt(area);
    center = quad.centroid(area);

    for (int i = 0; i < 4; ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqQuad& quad)
    : AmrCell(quad.reduce()) {
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

AmrCell::AmrCell(const Cube& cube)
    : Element(0, 0), dim(3),
    adaptive(true), linear(true),
    vertices(cube), faces(CellType::AMR3D),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double volume = cube.volume();
    size = std::cbrt(volume);
    center = cube.centroid(volume);

    for (int i = 0; i < 6; ++i) {
        Quad vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]],
                vertices[faces[i].vertices[2]],
                vertices[faces[i].vertices[3]]
        };

        faces[i].area     = vs.area();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqCube& cube)
    : AmrCell(cube.reduce()) {
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

AmrCell::AmrCell(const Polygon& poly)
        : Element(0, 0), dim(2),
          adaptive(false), linear(true),
          vertices(poly), faces(CellType::POLYGON, poly.size()),
          b_idx(-1), z_idx(-1), level(0), flag(0) {

    double area = poly.area();

    size   = std::sqrt(area);
    center = poly.centroid(area);

    for (int i = 0; i < poly.size(); ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(Polyhedron poly)
        : Element(0, 0), dim(3),
          adaptive(false), linear(true),
          vertices(poly), faces(CellType::POLYHEDRON),
          b_idx(-1), z_idx(-1), level(0), flag(0)  {

    double vol = poly.volume();

    size   = std::cbrt(vol);
    center = poly.centroid(vol);

    if (poly.n_faces() > BFaces::max_count) {
        throw std::runtime_error("Polygon has > 24 faces");
    }

    for (int i = 0; i < poly.n_faces(); ++i) {
        int n_verts = poly.face_indices(i).size();
        if (n_verts < 4) {
            faces[i].vertices = {poly.face_indices(i)[0],
                                 poly.face_indices(i)[1],
                                 poly.face_indices(i)[2],
                                 -1
            };
        }
        else if (n_verts < 5) {
            faces[i].vertices = {poly.face_indices(i)[0],
                                 poly.face_indices(i)[1],
                                 poly.face_indices(i)[2],
                                 poly.face_indices(i)[3]
            };
        }
        else {
            // Хардкор, полигональная грань
            faces[i].set_polygonal();
            for (int j = 0; j < n_verts; ++j) {
                faces[i].set_poly_vertex(j, poly.face_indices(i)[j]);
            }
        }

        faces[i].area     = poly.face_area(i);
        faces[i].center   = poly.face_center(i);
        faces[i].normal   = poly.face_normal(i);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

double AmrCell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

Polygon AmrCell::polygon() const {
    if (dim > 2) {
        throw std::runtime_error("AmrCell::polygon() error #1");
    }

    if (adaptive) {
        Polygon poly;
        poly.reserve(8);

        poly += vertices.vs<-1, -1>();
        if (!linear && faces[Side::BOTTOM1].is_actual()) {
            poly += vertices.vs<0, -1>();
        }

        poly += vertices.vs<+1, -1>();
        if (!linear && faces[Side::RIGHT1].is_actual()) {
            poly += vertices.vs<+1, 0>();
        }

        poly += vertices.vs<+1, +1>();
        if (!linear && faces[Side::TOP1].is_actual()) {
            poly += vertices.vs<0, +1>();
        }

        poly += vertices.vs<-1, +1>();
        if (!linear && faces[Side::LEFT1].is_actual()) {
            poly += vertices.vs<-1, 0>();
        }
        return poly;
    } else {
        return Polygon(&vertices[0], vertices.count());
    }
}

Polyhedron AmrCell::polyhedron() const {
    throw std::runtime_error("AmrCell::polyhedron");
}

inline double min(double x, double y, double z) {
    return std::min(x, std::min(y, z));
}

double AmrCell::incircle_diameter() const {
    if (adaptive) {
        if (dim == 2) {
            return std::sqrt(std::min(
                    (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                    (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
        }
        else {
            return std::sqrt(min(
                    (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                    (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                    (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()));
        }
    }
    else {
        assert(dim == 2);
        int n = vertices.count();
        // Диаметр описаной окружности вокруг правильного многоугольника
        // с площадью size^2.
        return 2.0 * size / std::sqrt(n * std::tan(M_PI / n));
    }
}

double AmrCell::approx_vol_fraction(const std::function<double(const Vector3d &)> &inside) const {
    if (dim < 3) {
        if (adaptive) {
            int sum = 0;
            // Угловые точки, вес = 1
            if (inside(vertices.vs<-1, -1>()))
                sum += 1;
            if (inside(vertices.vs<-1, +1>()))
                sum += 1;
            if (inside(vertices.vs<+1, +1>()))
                sum += 1;
            if (inside(vertices.vs<+1, -1>()))
                sum += 1;

            // Ребра, вес = 2
            if (inside(vertices.vs<0, -1>()))
                sum += 2;
            if (inside(vertices.vs<0, +1>()))
                sum += 2;
            if (inside(vertices.vs<-1, 0>()))
                sum += 2;
            if (inside(vertices.vs<+1, 0>()))
                sum += 2;

            // Центр, вес = 3
            if (inside(vertices.vs<0, 0>()))
                sum += 4;

            if (sum == 0) {
                return 0.0;
            }
            else if (sum == 16) {
                return 1.0;
            }
            else {
                return 0.0625 * sum; // sum / 16.0
            }
        }
        else {
            // Не адаптивная ячейка
            int sum = 0;
            int count = 0;
            for (; count < BVertices::max_count; ++count) {
                const Vector3d &v = vertices[count];

                // Дошли до конца
                if (v.hasNaN()) {
                    break;
                } else {
                    // Вершины многоугольника, вес 2
                    if (inside(v)) {
                        sum += 2;
                    }
                }
            }
            // Центр многоугольника, вес равен числу вершин
            if (inside(center)) {
                sum += count;
            }
            return count < 1 ? 0.0 : sum / (3.0 * count);
        }
    }
    else {
        // Трехмерная ячейка
        throw std::runtime_error("AmrCell::approx_vol_fraction #1");
    }
}

double AmrCell::volume_fraction(const std::function<double(const Vector3d &)> &inside, int n_points) const {
    if (dim < 3) {
        if (adaptive) {
            if (linear) {
                return vertices.as2D().reduce().volume_fraction(inside, n_points);
            }
            else {
                return vertices.as2D().volume_fraction(inside, n_points);
            }
        }
        else {
            // Полигон
            int count = vertices.count();
            int N = n_points / count + 1;

            double res = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;
                Triangle tri(center, vertices[i], vertices[j]);
                res += tri.volume_fraction(inside, N) * tri.area();
            }

            return res / volume();
        }
    }
    else {
        // Трехмерная ячейка
        throw std::runtime_error("AmrCell::volume_fraction #1");
    }
}

bool AmrCell::const_function(const std::function<double(const Vector3d&)>& func) const {
    double value = func(center);
    if (dim < 3) {
        if (adaptive) {            
            // Угловые точки
            if (func(vertices.vs<-1, -1>()) != value) { return false; }
            if (func(vertices.vs<-1, +1>()) != value) { return false; }               
            if (func(vertices.vs<+1, +1>()) != value) { return false; }               
            if (func(vertices.vs<+1, -1>()) != value) { return false; }

            // Ребра
            if (func(vertices.vs<0, -1>()) != value) { return false; }                
            if (func(vertices.vs<0, +1>()) != value) { return false; }               
            if (func(vertices.vs<-1, 0>()) != value) { return false; }               
            if (func(vertices.vs<+1, 0>()) != value) { return false; }               

            // Центр
            return func(vertices.vs<0, 0>()) == value;
        }
        else {
            // Не адаптивная ячейка
            int i = 0;
            for (; i < BVertices::max_count; ++i) {
                if (vertices[i].hasNaN()) {
                    // Дошли до конца
                    return true;
                } else {
                    if (func(vertices[i]) != value) {
                        return false;
                    }
                }
            }
            return true;
        }
    }
    else {
        // Трехмерная ячейка
        throw std::runtime_error("AmrCell::const_function #1");
    }    
}

double AmrCell::integrate_low(const std::function<double(const Vector3d&)>& func, int n_points) const {
    if (dim < 3) {
        if (adaptive) {
            if (linear) {
                return vertices.as2D().reduce().integrate_low(func, n_points);
            }
            else {
                return vertices.as2D().integrate_low(func, n_points);
            }
        }
        else {
            // Полигон
            int count = vertices.count();
            int N = n_points / count + 1;

            double sum = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;
                Triangle tri(center, vertices[i], vertices[j]);
                sum += tri.integrate_low(func, N);
            }

            return sum;
        }
    }
    else {
        // Трехмерная ячейка
        throw std::runtime_error("AmrCell::volume_fraction #1");
    }
}

void AmrCell::mark_actual_nodes(int mark) {
    static_assert(BNodes::max_count == BVertices::max_count);

    nodes.clear();
    for (const BFace &face: faces) {
        for (int k = 0; k < BFace::max_size; ++k) {
            int iv = face.vertices[k];
            if (iv >= 0) {
                nodes[iv] = -13;
            }
            else {
                break;
            }
        }
    }
}

} // namespace zephyr::geom