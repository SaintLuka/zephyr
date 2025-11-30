#include <iostream>

#include <zephyr/geom/geom.h>
#include <zephyr/math/funcs.h>
#include <zephyr/mesh/euler/amr_cells.h>

using namespace zephyr::geom;

namespace zephyr::mesh {

AmrCells::AmrCells(int dim, bool adaptive, bool axial)
        : m_size(0) {

    // На единицу больше числа ячеек
    face_begin = {0};
    node_begin = {0};

    set_dimension(dim);
    set_adaptive(adaptive);
    set_axial(axial);
    set_linear(true);
}

AmrCells AmrCells::same() const {
    AmrCells cells(m_dim, m_adaptive, m_axial);
    cells.data = data.same();
    return cells;
}

void AmrCells::set_dimension(int dim) {
    if (!empty() && dim != m_dim) {
        throw std::runtime_error("Can't change dimension. Mesh is not empty.");
    }

    if (dim == 2) {
        m_dim = 2;
    }
    else if (dim == 3) {
        m_dim = 3;
        m_axial = false;
        m_linear = true;
    }
    else {
        throw std::runtime_error("Mesh dimension should be equal 2 or 3");
    }
}

void AmrCells::set_adaptive(bool adaptive) {
    if (!empty() && adaptive != m_adaptive) {
        throw std::runtime_error("Can't change 'adaptive'. Mesh is not empty.");
    }

    m_adaptive = adaptive;
}

void AmrCells::set_axial(bool axial) {
    if (!empty() && axial != m_axial) {
        throw std::runtime_error("Can't change symmetry. Mesh is not empty.");
    }

    m_axial = axial;
    if (axial) {
        m_dim = 2;
    }
}

void AmrCells::set_linear(bool linear) {
    m_linear = linear;
}

int AmrCells::face_count(index_t ic) const {
    if (m_adaptive) {
        int count = 0;
        for (index_t iface: faces_range(ic)) {
            if (faces.is_actual(iface)) {
                ++count;
            }
        }
        return count;
    }
    else {
        return max_face_count(ic);
    }
}

double AmrCells::hx(index_t ic) const {
    z_assert(m_adaptive, "Not adaptive mesh, can't get 'hx' for cell");
    return (mapping<2>(ic).vs<+1, 0>() - mapping<2>(ic).vs<-1, 0>()).norm();
}

double AmrCells::hy(index_t ic) const {
    z_assert(m_adaptive, "Not adaptive mesh, can't get 'hy' for cell");
    return (mapping<2>(ic).vs<0, +1>() - mapping<2>(ic).vs<0, -1>()).norm();
}

double AmrCells::hz(index_t ic) const {
    z_assert(m_adaptive, "Not adaptive mesh, can't get 'hz' for cell");
    z_assert(m_dim == 3, "Two dimensional mesh, can't get 'hz' for cell");
    return (mapping<3>(ic).vs<0, 0, +1>() - mapping<3>(ic).vs<0, 0, -1>()).norm();
}

double AmrCells::incircle_diameter(index_t ic) const {
    if (m_adaptive) {
        if (m_dim == 2) {
            const SqQuad &vertices = mapping<2>(ic);
            return std::sqrt(std::min(
                    (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                    (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
        } else {
            const SqCube &vertices = mapping<3>(ic);
            return std::sqrt(math::min(
                    (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                    (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                    (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()));
        }
    } else {
        assert(m_dim == 2);
        int n = node_count(ic);
        // Диаметр вписанной окружности внутрь правильного многоугольника
        // с площадью volume.
        return 2.0 * std::sqrt(volume[ic] / (n * std::tan(M_PI / n)));
    }
}

Box AmrCells::bbox(index_t ic) const {
    // TODO: Сделать оптимальный код для декартовых сеток, и не только здесь
    Box box = Box::Empty(m_dim);
    for (index_t iv: nodes_range(ic)) {
        box.capture(verts[iv]);
    }
    return box;
}

Polygon AmrCells::polygon(index_t ic) const {
    if (m_dim > 2) {
        throw std::runtime_error("AmrCell::polygon() error #1");
    }

    if (m_adaptive) {
        Polygon poly;
        poly.reserve(8);

        const SqQuad& vertices = mapping<2>(ic);

        poly += vertices.vs<-1, -1>();
        if (!m_linear && complex_face(ic, Side2D::BOTTOM)) {
            poly += vertices.vs<0, -1>();
        }

        poly += vertices.vs<+1, -1>();
        if (!m_linear && complex_face(ic, Side2D::RIGHT)) {
            poly += vertices.vs<+1, 0>();
        }

        poly += vertices.vs<+1, +1>();
        if (!m_linear && complex_face(ic, Side2D::TOP)) {
            poly += vertices.vs<0, +1>();
        }

        poly += vertices.vs<-1, +1>();
        if (!m_linear && complex_face(ic, Side2D::LEFT)) {
            poly += vertices.vs<-1, 0>();
        }

        return poly;
    } else {
        return {vertices_data(ic), node_count(ic)};
    }
}

Polyhedron AmrCells::polyhedron(index_t ic) const {
    if (m_dim < 3) {
        throw std::runtime_error("AmrCell::polyhedron() error #1");
    }

    if (m_adaptive && m_linear) {
        // Пока я умею делать только кубы
        const auto& map = mapping<3>(ic);
        std::vector<Vector3d> vs = {
            map.vs<-1, -1, -1>(),
            map.vs<+1, -1, -1>(),
            map.vs<+1, +1, -1>(),
            map.vs<-1, +1, -1>(),
            map.vs<-1, -1, +1>(),
            map.vs<+1, -1, +1>(),
            map.vs<+1, +1, +1>(),
            map.vs<-1, +1, +1>(),
        };
        return Polyhedron(CellType::HEXAHEDRON, vs);
    }

    throw std::runtime_error("AmrCell::polyhedron() error #2");
}

double AmrCells::approx_vol_fraction(index_t ic, const InFunction &inside) const {
    if (m_dim < 3) {
        if (m_adaptive) {
            const SqQuad& vertices = mapping<2>(ic);

            int sum = 0;
            // Угловые точки, вес = 1
            if (inside(vertices.vs<-1, -1>())) sum += 1;
            if (inside(vertices.vs<-1, +1>())) sum += 1;
            if (inside(vertices.vs<+1, +1>())) sum += 1;
            if (inside(vertices.vs<+1, -1>())) sum += 1;

            // Ребра, вес = 2
            if (inside(vertices.vs<0, -1>())) sum += 2;
            if (inside(vertices.vs<0, +1>())) sum += 2;
            if (inside(vertices.vs<-1, 0>())) sum += 2;
            if (inside(vertices.vs<+1, 0>())) sum += 2;

            // Центр, вес = 3
            if (inside(vertices.vs<0, 0>())) sum += 4;

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
            int count = node_count(ic);

            int sum = 0;
            for (auto i: nodes_range(ic)) {
                // Вершины многоугольника, вес 2
                if (inside(verts[i])) {
                    sum += 2;
                }
            }
            // Центр многоугольника, вес равен числу вершин
            if (inside(center[ic])) {
                sum += count;
            }
            return sum < 1 ? 0.0 : sum / (3.0 * count);
        }
    }
    else {
        // Трехмерная ячейка
        if (m_adaptive) {
            const SqCube& vertices = mapping<3>(ic);

            int sum = 0;
            // Угловые точки, вес = 1
            if (inside(vertices.vs<-1, -1, -1>())) sum += 1;
            if (inside(vertices.vs<+1, -1, -1>())) sum += 1;
            if (inside(vertices.vs<-1, +1, -1>())) sum += 1;
            if (inside(vertices.vs<+1, +1, -1>())) sum += 1;
            if (inside(vertices.vs<-1, -1, +1>())) sum += 1;
            if (inside(vertices.vs<+1, -1, +1>())) sum += 1;
            if (inside(vertices.vs<-1, +1, +1>())) sum += 1;
            if (inside(vertices.vs<+1, +1, +1>())) sum += 1;

            // Ребра, вес = 2
            if (inside(vertices.vs<-1, -1, 0>())) sum += 2;
            if (inside(vertices.vs<+1, -1, 0>())) sum += 2;
            if (inside(vertices.vs<-1, +1, 0>())) sum += 2;
            if (inside(vertices.vs<+1, +1, 0>())) sum += 2;
            if (inside(vertices.vs<-1, 0, -1>())) sum += 2;
            if (inside(vertices.vs<+1, 0, -1>())) sum += 2;
            if (inside(vertices.vs<-1, 0, +1>())) sum += 2;
            if (inside(vertices.vs<+1, 0, +1>())) sum += 2;
            if (inside(vertices.vs<0, -1, -1>())) sum += 2;
            if (inside(vertices.vs<0, +1, -1>())) sum += 2;
            if (inside(vertices.vs<0, -1, +1>())) sum += 2;
            if (inside(vertices.vs<0, +1, +1>())) sum += 2;

            // Грани, вес = 4
            if (inside(vertices.vs<-1, 0, 0>())) sum += 4;
            if (inside(vertices.vs<+1, 0, 0>())) sum += 4;
            if (inside(vertices.vs<0, -1, 0>())) sum += 4;
            if (inside(vertices.vs<0, +1, 0>())) sum += 4;
            if (inside(vertices.vs<0, 0, -1>())) sum += 4;
            if (inside(vertices.vs<0, 0, +1>())) sum += 4;

            // Центр, вес = 8
            if (inside(vertices.vs<0, 0, 0>())) sum += 8;

            if (sum == 0) {
                return 0.0;
            }
            else if (sum == 64) {
                return 1.0;
            }
            else {
                return 0.015625 * sum; // sum / 64.0
            }
        }
        else {
            throw std::runtime_error("Approx volume fraction error #1");
        };
    }
}

double AmrCells::volume_fraction(index_t ic, const InFunction &inside, int n_points) const {
    if (m_dim < 3) {
        if (m_adaptive) {
            if (m_linear) {
                return mapping<2>(ic).reduce().volume_fraction(inside, n_points);
            }
            else {
                return mapping<2>(ic).volume_fraction(inside, n_points);
            }
        }
        else {
            // Полигон
            int count = node_count(ic);
            int N = n_points / count + 1;

            double res = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;

                index_t I = node_begin[ic] + i;
                index_t J = node_begin[ic] + j;

                Triangle tri(center[ic], verts[I], verts[J]);
                res += tri.volume_fraction(inside, N) * tri.area();
            }

            return res / volume[ic];
        }
    }
    else {
        if (m_adaptive) {
            return mapping<3>(ic).reduce().volume_fraction(inside, n_points);
        }
        // Трехмерный многогранник
        throw std::runtime_error("AmrCell::volume_fraction #1");
    }
}

bool AmrCells::const_function(index_t ic, const SpFunction& func) const {
    double value = func(center[ic]);
    if (m_dim < 3) {
        if (m_adaptive) {
            const SqQuad& vertices = mapping<2>(ic);

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
            for (auto i: nodes_range(ic)) {
                if (func(verts[i]) != value) {
                    return false;
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

double AmrCells::integrate_low(index_t ic, const SpFunction& func, int n_points) const {
    if (m_dim < 3) {
        if (m_adaptive) {
            if (m_linear) {
                return mapping<2>(ic).reduce().integrate_low(func, n_points);
            }
            else {
                return mapping<2>(ic).integrate_low(func, n_points);
            }
        }
        else {
            // Полигон
            int count = node_count(ic);
            int N = n_points / count + 1;

            double sum = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;

                index_t I = node_begin[ic] + i;
                index_t J = node_begin[ic] + j;

                Triangle tri(center[ic], verts[I], verts[J]);
                sum += tri.integrate_low(func, N) * tri.area();
            }

            return sum;
        }
    }
    else {
        // Трехмерная ячейка
        if (m_adaptive) {
            return mapping<3>(ic).reduce().integrate_low(func, n_points);
        }
        throw std::runtime_error("AmrCell::volume_fraction #1");
    }
}

void AmrCells::move_item(index_t from, index_t to) {
    rank[to] = rank[from];
    next[to] = to;
    index[to] = to;

    flag[to] = flag[from];
    level[to] = level[from];
    b_idx[to] = b_idx[from];
    z_idx[to] = z_idx[from];

    center[to] = center[from];
    volume[to] = volume[from];

    volume_alt[to] = volume_alt[from];

    for (index_t i = 0; i < face_begin[from + 1] - face_begin[from]; ++i) {
        index_t iface = face_begin[from] + i;
        index_t jface = face_begin[to] + i;

        faces.boundary[jface] = faces.boundary[iface];
        faces.normal  [jface] = faces.normal  [iface];
        faces.center  [jface] = faces.center  [iface];
        faces.area    [jface] = faces.area    [iface];
        faces.area_alt[jface] = faces.area_alt[iface];
        faces.vertices[jface] = faces.vertices[iface];

        faces.adjacent.rank [jface] = faces.adjacent.rank[iface];
        faces.adjacent.index[jface] = faces.adjacent.index[iface];
        faces.adjacent.alien[jface] = faces.adjacent.alien[iface];
        faces.adjacent.basic[jface] = to;
    }

    for (index_t i = 0; i < node_begin[from + 1] - node_begin[from]; ++i) {
        index_t jv = node_begin[to] + i;
        index_t iv = node_begin[from] + i;
        verts[jv] = verts[iv];
    }

    set_undefined(from);
}

void AmrCells::copy_data(index_t from, index_t to) {
    copy_data(from, this, to);
}

void AmrCells::copy_data(index_t from, AmrCells* dst, index_t to) const {
    data.copy_data(from, &dst->data, to);
}

void AmrCells::copy_geom(index_t ic, AmrCells& cells,
        index_t jc, index_t face_beg, index_t node_beg) const {

    cells.rank [jc] = rank [ic];
    cells.next [jc] = next [ic];
    cells.index[jc] = index[ic];

    cells.flag [jc] = flag [ic];
    cells.level[jc] = level[ic];
    cells.b_idx[jc] = b_idx[ic];
    cells.z_idx[jc] = z_idx[ic];

    cells.center[jc] = center[ic];
    cells.volume[jc] = volume[ic];

    cells.volume_alt[jc] = volume_alt[ic];

    cells.face_begin[jc] = face_beg;
    cells.face_begin[jc + 1] = face_beg + max_face_count(ic);

    for (index_t i = 0; i < max_face_count(ic); ++i) {
        index_t iface = face_begin[ic] + i;
        index_t jface = cells.face_begin[jc] + i;

        cells.faces.boundary[jface] = faces.boundary[iface];
        cells.faces.normal  [jface] = faces.normal  [iface];
        cells.faces.center  [jface] = faces.center  [iface];
        cells.faces.area    [jface] = faces.area    [iface];
        cells.faces.area_alt[jface] = faces.area_alt[iface];
        cells.faces.vertices[jface] = faces.vertices[iface];

        cells.faces.adjacent.rank [jface] = faces.adjacent.rank [iface];
        cells.faces.adjacent.index[jface] = faces.adjacent.index[iface];
        cells.faces.adjacent.alien[jface] = faces.adjacent.alien[iface];
        cells.faces.adjacent.basic[jface] = index[ic];
    }

    cells.node_begin[jc] = node_beg;
    cells.node_begin[jc + 1] = node_beg + max_node_count(ic);

    for (index_t i = 0; i < max_node_count(ic); ++i) {
        index_t iv = node_begin[ic] + i;
        index_t jv = cells.node_begin[jc] + i;
        cells.verts[jv] = verts[iv];
    }
}

void AmrCells::copy_geom_basic(index_t ic, AmrCells& cells,
        index_t jc, index_t face_beg, index_t node_beg) const {

    cells.rank [jc] = rank [ic];
    cells.index[jc] = index[ic];

    cells.face_begin[jc] = face_beg;
    cells.face_begin[jc + 1] = face_beg + max_face_count(ic);

    for (index_t i = 0; i < max_face_count(ic); ++i) {
        index_t iface = face_begin[ic] + i;
        index_t jface = cells.face_begin[jc] + i;

        cells.faces.boundary[jface] = faces.boundary[iface];

        cells.faces.adjacent.rank [jface] = faces.adjacent.rank [iface];
        cells.faces.adjacent.index[jface] = faces.adjacent.index[iface];
        cells.faces.adjacent.alien[jface] = faces.adjacent.alien[iface];
    }

    cells.node_begin[jc] = node_beg;
    cells.node_begin[jc + 1] = node_beg + max_node_count(ic);
}

void AmrCells::clear() {
    resize(0, 0, 0);
}

void AmrCells::resize_amr(index_t n_cells) {
    if (!m_adaptive) {
        throw std::runtime_error("Resize of unstructured mesh");
    }

    index_t n_faces = n_cells * (m_dim == 2 ? 8 : 24);
    index_t n_nodes = n_cells * (m_dim == 2 ? 9 : 27);
    
    resize(n_cells, n_faces, n_nodes);
}

void AmrCells::reserve_amr(index_t n_cells) {
    if (!m_adaptive) {
        throw std::runtime_error("Resize of unstructured mesh");
    }

    index_t n_faces = n_cells * (m_dim == 2 ? 8 : 24);
    index_t n_nodes = n_cells * (m_dim == 2 ? 9 : 27);

    reserve(n_cells, n_faces, n_nodes);
}

void AmrCells::resize_cells(index_t n_cells) {
    m_size = n_cells;

    // Поля ячеек по числу ячеек, логично
    next.resize(n_cells, -1);
    rank.resize(n_cells, -1);
    index.resize(n_cells, -1);

    center.resize(n_cells);
    volume.resize(n_cells);
    volume_alt.resize(n_cells);

    // +1 для заключительной
    face_begin.resize(n_cells + 1);
    node_begin.resize(n_cells + 1);

    flag.resize(n_cells);
    b_idx.resize(n_cells);
    z_idx.resize(n_cells);
    level.resize(n_cells);

    // Поля данных только для ячеек
    data.resize(n_cells);
}

void AmrCells::reserve_cells(index_t n_cells) {
    // Поля ячеек по числу ячеек, логично
    next.reserve(n_cells);
    rank.reserve(n_cells);
    index.reserve(n_cells);

    center.reserve(n_cells);
    volume.reserve(n_cells);
    volume_alt.reserve(n_cells);

    // +1 для заключительной
    face_begin.reserve(n_cells + 1);
    node_begin.reserve(n_cells + 1);

    flag.reserve(n_cells);
    b_idx.reserve(n_cells);
    z_idx.reserve(n_cells);
    level.reserve(n_cells);

    // Поля данных только для ячеек
    data.reserve(n_cells);
}

void AmrCells::shrink_to_fit_cells() {
    // Поля ячеек по числу ячеек, логично
    next.shrink_to_fit();
    rank.shrink_to_fit();
    index.shrink_to_fit();

    center.shrink_to_fit();
    volume.shrink_to_fit();
    volume_alt.shrink_to_fit();

    // +1 для заключительной
    face_begin.shrink_to_fit();
    node_begin.shrink_to_fit();

    flag.shrink_to_fit();
    b_idx.shrink_to_fit();
    z_idx.shrink_to_fit();
    level.shrink_to_fit();

    // Поля данных только для ячеек
    data.shrink_to_fit();
}

void AmrCells::resize(index_t n_cells, index_t n_faces, index_t n_nodes) {
    resize_cells(n_cells);
    faces.resize(n_faces);
    verts.resize(n_nodes);
}

void AmrCells::reserve(index_t n_cells, index_t n_faces, index_t n_nodes) {
    reserve_cells(n_cells);
    faces.reserve(n_faces);
    verts.reserve(n_nodes);
}

void AmrCells::shrink_to_fit() {
    shrink_to_fit_cells();
    faces.shrink_to_fit();
    verts.shrink_to_fit();
}

void AmrCells::set_cell(index_t ic, const Quad& quad) {
    assert(m_dim == 2);
    assert(m_adaptive);
    assert(m_linear);
    assert(!m_axial);

    rank[ic] = -1;
    index[ic] = -1;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = quad.area();
    center[ic] = quad.centroid(volume[ic]);

    face_begin[ic] = 8 * ic;
    face_begin[ic + 1] = 8 * (ic + 1);
    faces.insert(face_begin[ic], CellType::AMR2D);

    node_begin[ic] = 9 * ic;
    node_begin[ic + 1] = 9 * (ic + 1);
    mapping<2>(ic) = SqQuad(quad);

    for (index_t iface = face_begin[ic]; iface < face_begin[ic] + Side2D::count(); ++iface) {
        Line vs = {
                verts[node_begin[ic] + faces.vertices[iface][0]],
                verts[node_begin[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.center();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::ORDINARY;
    }
}

void AmrCells::set_cell(index_t ic, const Quad& quad, bool axial) {
    assert(m_dim == 2);
    assert(m_adaptive);
    assert(m_linear);
    assert(m_axial == axial);

    rank[ic] = -1;
    index[ic] = -1;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = quad.area();
    center[ic] = quad.centroid(volume[ic]);

    volume[ic]     = quad.area();
    volume_alt[ic] = quad.volume_as();
    center[ic]     = quad.centroid_as(volume_alt[ic]);

    volume[ic] = quad.area();
    if (!axial) {
        center[ic]     = quad.centroid(volume[ic]);
    }
    else {
        volume_alt[ic] = quad.volume_as();
        center[ic]     = quad.centroid_as(volume_alt[ic]);
    }

    face_begin[ic] = 8 * ic;
    face_begin[ic + 1] = 8 * (ic + 1);
    faces.insert(face_begin[ic], CellType::AMR2D);

    node_begin[ic] = 9 * ic;
    node_begin[ic + 1] = 9 * (ic + 1);
    mapping<2>(ic) = SqQuad(quad);

    for (auto i: nodes_range(ic)) {
        if (verts[i].z() != 0.0) {
            throw std::runtime_error("AmrCells add axial cell, vertex.z != 0.0");
        }
    }

    for (index_t iface = face_begin[ic]; iface < face_begin[ic] + Side2D::count(); ++iface) {
        Line vs = {
                verts[node_begin[ic] + faces.vertices[iface][0]],
                verts[node_begin[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.centroid(axial);
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::ORDINARY;

        if (axial) {
            faces.area_alt[iface] = vs.area_as();
        }
    }
}

void AmrCells::set_cell(index_t ic, const SqQuad& quad) {
    set_cell(ic, quad.reduce());
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

void AmrCells::set_cell(index_t ic, const SqQuad& quad, bool axial) {
    set_cell(ic, quad.reduce(), axial);
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

void AmrCells::set_cell(index_t ic, const Cube& cube) {
    assert(m_dim == 3);
    assert(m_adaptive);
    assert(m_linear);
    assert(!m_axial);

    rank[ic] = -1;
    index[ic] = -1;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = cube.volume();
    center[ic] = cube.centroid(volume[ic]);

    face_begin[ic] = 24 * ic ;
    face_begin[ic + 1] = 24 * (ic + 1);
    faces.insert(face_begin[ic], CellType::AMR3D);

    node_begin[ic] = 27 * ic;
    node_begin[ic + 1] = 27 * (ic + 1);
    mapping<3>(ic) = SqCube(cube);

    for (index_t iface = face_begin[ic]; iface < face_begin[ic] + Side3D::count(); ++iface) {
        Quad vs = {
                verts[node_begin[ic] + faces.vertices[iface][0]],
                verts[node_begin[ic] + faces.vertices[iface][1]],
                verts[node_begin[ic] + faces.vertices[iface][2]],
                verts[node_begin[ic] + faces.vertices[iface][3]]
        };

        faces.area[iface]     = vs.area();
        faces.center[iface]   = vs.center();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::ORDINARY;
    }
}

void AmrCells::set_cell(index_t ic, const SqCube& cube) {
    set_cell(ic, cube.reduce());
}

void AmrCells::push_back(const geom::Line &line) {
    throw std::runtime_error("NO WAY ACPBWER");
}

void AmrCells::push_back(const Polygon& poly) {
    assert(m_dim == 2);
    assert(!m_adaptive);
    assert(m_linear);

    index_t ic = size();

    resize_cells(ic + 1);

    rank[ic] = 0;
    index[ic] = ic;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = poly.area();
    center[ic] = poly.centroid(volume[ic]);

    if (m_axial) {
        volume_alt[ic] = poly.volume_as();
    }

    int n_nodes = poly.size();
    int n_faces = poly.size();

    faces.resize(faces.size() + n_faces);
    verts.resize(verts.size() + n_nodes);

    face_begin[ic + 1] = face_begin[ic] + n_faces;
    faces.insert(face_begin[ic], CellType::POLYGON, n_faces);

    node_begin[ic + 1] = node_begin[ic] + n_nodes;
    for (int i = 0; i < n_nodes; ++i) {
        verts[node_begin[ic] + i] = poly[i];
    }

    for (index_t iface = face_begin[ic]; iface < face_begin[ic] + n_faces; ++iface) {
        Line vs = {
                verts[node_begin[ic] + faces.vertices[iface][0]],
                verts[node_begin[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.centroid();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::ORDINARY;

        if (m_axial) {
            faces.area_alt[iface] = vs.area_as();
        }
    }
}

void AmrCells::push_back(const Polyhedron& poly) {
    if (poly.need_simplify(AmrFaces::max_vertices)) {
        Polyhedron simple_poly = poly;
        simple_poly.simplify_faces(AmrFaces::max_vertices);
        push_back_impl(simple_poly);
    }
    else {
        push_back_impl(poly);
    }
}

void AmrCells::push_back_impl(const Polyhedron& poly) {
    assert(m_dim == 3);
    assert(!m_adaptive);
    assert(m_linear);
    assert(!m_axial);

    index_t ic = size();

    resize_cells(ic + 1);

    rank[ic] = 0;
    index[ic] = ic;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = poly.volume();
    center[ic] = poly.centroid(volume[ic]);

    volume_alt[ic] = NAN;

    // Зададим вершины многогранника
    int n_nodes = poly.n_verts();
    verts.resize(verts.size() + n_nodes);

    node_begin[ic + 1] = node_begin[ic] + n_nodes;
    for (int i = 0; i < n_nodes; ++i) {
        verts[node_begin[ic] + i] = poly.vertex(i);
    }

    // Определим грани многогранника
    int n_faces = poly.n_faces();
    faces.resize(faces.size() + n_faces);

    face_begin[ic + 1] = face_begin[ic] + n_faces;
    for (int i = 0; i < poly.n_faces(); ++i) {
        index_t iface = face_begin[ic] + i;

        faces.area[iface] = poly.face_area(i);
        faces.center[iface] = poly.face_center(i);
        faces.normal[iface] = poly.face_normal(i);
        faces.boundary[iface] = Boundary::ORDINARY;

        faces.area_alt[iface] = NAN;

        // Выставить вершины грани
        int n_verts = poly.face_indices(i).size();
        if (n_verts > AmrFaces::max_vertices) {
            std::string message = "Can't add polyhedron with " +
                    std::to_string(n_verts) + " vertices per face.";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }

        faces.vertices[iface].fill(-1);
        for (int j = 0; j < n_verts; ++j) {
            faces.vertices[iface][j] = poly.face_indices(i)[j];
        }

        faces.adjacent.basic[iface] = ic;
    }
}

} // namespace zephyr::mesh