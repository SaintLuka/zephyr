#include <filesystem>
#include <iostream>
#include <fstream>
#include <cassert>

#include <zephyr/geom/geom.h>
#include <zephyr/math/funcs.h>
#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/raw/raw_cells.h>

using namespace zephyr::geom;
using zephyr::utils::mpi;

namespace zephyr::mesh {

RawCells::RawCells(MeshOpts options) :
    dim_(options.dim),
    adaptive_(options.adaptive),
    linear_(options.linear),
    axial_(options.axial),
    verts(options.nodes) {

}

RawCells RawCells::same() const {
    RawCells cells(options());
    cells.data = data.same();
    return cells;
}

MeshOpts RawCells::options() const {
    return MeshOpts {
        .dim      = dim_,
        .adaptive = adaptive_,
        .linear   = linear_,
        .axial    = axial_,
        .nodes    = verts.has_nodes()
    };
}

int RawCells::face_count(index_t ic) const {
    if (adaptive_) {
        int count = 0;
        for (index_t iface: faces.range(ic)) {
            if (faces.is_actual(iface)) {
                ++count;
            }
        }
        return count;
    }
    else {
        return faces.max_count(ic);
    }
}

double RawCells::hx(index_t ic) const {
    z_assert(adaptive_, "Not adaptive mesh, can't get 'hx' for cell");
    return (verts.mapping<2>(ic).vs<+1, 0>() - verts.mapping<2>(ic).vs<-1, 0>()).norm();
}

double RawCells::hy(index_t ic) const {
    z_assert(adaptive_, "Not adaptive mesh, can't get 'hy' for cell");
    return (verts.mapping<2>(ic).vs<0, +1>() - verts.mapping<2>(ic).vs<0, -1>()).norm();
}

double RawCells::hz(index_t ic) const {
    z_assert(adaptive_, "Not adaptive mesh, can't get 'hz' for cell");
    z_assert(dim_ == 3, "Two dimensional mesh, can't get 'hz' for cell");
    return (verts.mapping<3>(ic).vs<0, 0, +1>() - verts.mapping<3>(ic).vs<0, 0, -1>()).norm();
}

double RawCells::incircle_diameter(index_t ic) const {
    if (adaptive_) {
        if (dim_ == 2) {
            const SqQuad &vertices = verts.mapping<2>(ic);
            return std::sqrt(std::min(
                    (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                    (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
        } else {
            const SqCube &vertices = verts.mapping<3>(ic);
            return std::sqrt(math::min(
                    (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                    (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                    (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()));
        }
    } else {
        if (dim_ == 2) {
            int n = verts.count(ic);
            // Диаметр вписанной окружности внутрь правильного многоугольника
            // с площадью volume.
            return 2.0 * std::sqrt(volume[ic] / (n * std::tan(M_PI / n)));
        }
        else {
            // Найдем минимальное расстояние до грани, умножим на два
            double r = std::numeric_limits<double>::infinity();
            for (auto j: faces.range(ic)) {
                double dist = std::abs((faces.center[j] - center[ic]).dot(faces.normal[j]));
                r = std::min(r, dist);
            }
            return 2.0 * r;
        }
    }
}

Box RawCells::bbox(index_t ic) const {
    // TODO: Сделать оптимальный код для декартовых сеток, и не только здесь
    Box box = Box::Empty(dim_);
    for (index_t iv: verts.range(ic)) {
        box.capture(verts[iv]);
    }
    return box;
}

Polygon RawCells::polygon(index_t ic) const {
    if (dim_ > 2) {
        throw std::runtime_error("RawCell::polygon() error #1");
    }

    if (adaptive_) {
        std::vector<Vector3d> poly;
        poly.reserve(8);

        const SqQuad& vertices = verts.mapping<2>(ic);

        poly.push_back(vertices.vs<-1, -1>());
        if (!linear_ && faces.is_complex(ic, Side2D::BOTTOM)) {
            poly.push_back(vertices.vs<0, -1>());
        }

        poly.push_back(vertices.vs<+1, -1>());
        if (!linear_ && faces.is_complex(ic, Side2D::RIGHT)) {
            poly.push_back(vertices.vs<+1, 0>());
        }

        poly.push_back(vertices.vs<+1, +1>());
        if (!linear_ && faces.is_complex(ic, Side2D::TOP)) {
            poly.push_back( vertices.vs<0, +1>());
        }

        poly .push_back(vertices.vs<-1, +1>());
        if (!linear_ && faces.is_complex(ic, Side2D::LEFT)) {
            poly.push_back(vertices.vs<-1, 0>());
        }

        return Polygon(std::move(poly));
    }
    return Polygon(std::span(verts.coords_data(ic), verts.count(ic)));
}

Polyhedron RawCells::polyhedron(index_t ic) const {
    if (dim_ < 3) {
        throw std::runtime_error("RawCell::polyhedron() error #1");
    }

    if (adaptive_ && linear_) {
        // Пока я умею делать только кубы
        const auto& map = verts.mapping<3>(ic);
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

    throw std::runtime_error("RawCell::polyhedron() error #2");
}

double RawCells::approx_vol_fraction(index_t ic, const InFunction &inside) const {
    if (dim_ < 3) {
        if (adaptive_) {
            const SqQuad& vertices = verts.mapping<2>(ic);

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
            int count = verts.count(ic);

            int sum = 0;
            for (auto i: verts.range(ic)) {
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
        if (adaptive_) {
            const SqCube& vertices = verts.mapping<3>(ic);

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
        }
    }
}

double RawCells::volume_fraction(index_t ic, const InFunction &inside, int n_points) const {
    if (dim_ < 3) {
        if (adaptive_) {
            if (linear_) {
                return verts.mapping<2>(ic).reduce().volume_fraction(inside, n_points);
            }
            else {
                return verts.mapping<2>(ic).volume_fraction(inside, n_points);
            }
        }
        else {
            // Полигон
            int count = verts.count(ic);
            int N = n_points / count + 1;

            double res = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;

                index_t I = verts.offsets[ic] + i;
                index_t J = verts.offsets[ic] + j;

                Triangle tri(center[ic], verts[I], verts[J]);
                res += tri.volume_fraction(inside, N) * tri.area();
            }

            return res / volume[ic];
        }
    }
    else {
        if (adaptive_) {
            return verts.mapping<3>(ic).reduce().volume_fraction(inside, n_points);
        }
        // Трехмерный многогранник
        throw std::runtime_error("RawCell::volume_fraction #1");
    }
}

bool RawCells::const_function(index_t ic, const SpFunction& func) const {
    double value = func(center[ic]);
    if (dim_ < 3) {
        if (adaptive_) {
            const SqQuad& vertices = verts.mapping<2>(ic);

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
            for (auto i: verts.range(ic)) {
                if (func(verts[i]) != value) {
                    return false;
                }
            }
            return true;
        }
    }
    else {
        // Трехмерная ячейка
        throw std::runtime_error("RawCell::const_function #1");
    }
}

double RawCells::integrate_low(index_t ic, const SpFunction& func, int n_points) const {
    if (dim_ < 3) {
        if (adaptive_) {
            if (linear_) {
                return verts.mapping<2>(ic).reduce().integrate_low(func, n_points);
            }
            else {
                return verts.mapping<2>(ic).integrate_low(func, n_points);
            }
        }
        else {
            // Полигон
            int count = verts.count(ic);
            int N = n_points / count + 1;

            double sum = 0.0;
            for (int i = 0; i < count; ++i) {
                int j = (i + 1) % count;

                index_t I = verts.offsets[ic] + i;
                index_t J = verts.offsets[ic] + j;

                Triangle tri(center[ic], verts[I], verts[J]);
                sum += tri.integrate_low(func, N) * tri.area();
            }

            return sum;
        }
    }
    else {
        // Трехмерная ячейка
        if (adaptive_) {
            return verts.mapping<3>(ic).reduce().integrate_low(func, n_points);
        }
        throw std::runtime_error("RawCell::volume_fraction #1");
    }
}

void RawCells::move_item(index_t from, index_t to) {
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

    for (index_t i = 0; i < faces.offsets[from + 1] - faces.offsets[from]; ++i) {
        index_t iface = faces.offsets[from] + i;
        index_t jface = faces.offsets[to] + i;

        faces.boundary[jface] = faces.boundary[iface];
        faces.normal  [jface] = faces.normal  [iface];
        faces.center  [jface] = faces.center  [iface];
        faces.area    [jface] = faces.area    [iface];
        faces.area_alt[jface] = faces.area_alt[iface];
        faces.vertices[jface] = faces.vertices[iface];

        faces.adjacent.rank [jface] = faces.adjacent.rank[iface];
        faces.adjacent.index[jface] = faces.adjacent.index[iface];
        faces.adjacent.ghost[jface] = faces.adjacent.ghost[iface];
        faces.adjacent.basic[jface] = to;
        faces.adjacent.rotation[jface] = faces.adjacent.rotation[iface];
    }

    for (index_t i = 0; i < verts.offsets[from + 1] - verts.offsets[from]; ++i) {
        index_t jv = verts.offsets[to] + i;
        index_t iv = verts.offsets[from] + i;
        verts[jv] = verts[iv];
    }

    set_undefined(from);
}

void RawCells::copy_data(index_t from, index_t to) {
    copy_data(from, this, to);
}

void RawCells::copy_data(index_t from, RawCells* dst, index_t to) const {
    data.copy_data(from, &dst->data, to);
}

void RawCells::copy_geom(index_t ic, RawCells& cells,
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

    cells.faces.offsets[jc] = face_beg;
    cells.faces.offsets[jc + 1] = face_beg + faces.max_count(ic);

    for (index_t i = 0; i < faces.max_count(ic); ++i) {
        index_t iface = faces.offsets[ic] + i;
        index_t jface = cells.faces.offsets[jc] + i;

        cells.faces.boundary[jface] = faces.boundary[iface];
        cells.faces.normal  [jface] = faces.normal  [iface];
        cells.faces.center  [jface] = faces.center  [iface];
        cells.faces.area    [jface] = faces.area    [iface];
        cells.faces.area_alt[jface] = faces.area_alt[iface];
        cells.faces.vertices[jface] = faces.vertices[iface];

        cells.faces.adjacent.rank [jface] = faces.adjacent.rank [iface];
        cells.faces.adjacent.index[jface] = faces.adjacent.index[iface];
        cells.faces.adjacent.ghost[jface] = faces.adjacent.ghost[iface];
        cells.faces.adjacent.basic[jface] = index[ic];
        cells.faces.adjacent.rotation[jface] = faces.adjacent.rotation[iface];
    }

    cells.verts.offsets[jc] = node_beg;
    cells.verts.offsets[jc + 1] = node_beg + verts.max_count(ic);

    z_assert(cells.has_nodes() == has_nodes(), "Different style cells");

    bool unique_nodes = cells.verts.has_nodes();
    for (index_t i = 0; i < verts.max_count(ic); ++i) {
        index_t iv = verts.offsets[ic] + i;
        index_t jv = cells.verts.offsets[jc] + i;
        cells.verts[jv] = verts[iv];
        if (unique_nodes) {
            cells.verts.rank [jv] = verts.rank [iv];
            cells.verts.index[jv] = verts.index[iv];
            cells.verts.ghost[jv] = verts.ghost[iv];
        }
    }
}

void RawCells::copy_geom_basic(index_t ic, RawCells& cells,
        index_t jc, index_t face_beg, index_t node_beg) const {

    cells.rank [jc] = rank [ic];
    cells.index[jc] = index[ic];

    cells.faces.offsets[jc] = face_beg;
    cells.faces.offsets[jc + 1] = face_beg + faces.max_count(ic);

    for (index_t i = 0; i < faces.max_count(ic); ++i) {
        index_t iface = faces.offsets[ic] + i;
        index_t jface = cells.faces.offsets[jc] + i;

        cells.faces.boundary[jface] = faces.boundary[iface];

        cells.faces.adjacent.rank [jface] = faces.adjacent.rank [iface];
        cells.faces.adjacent.index[jface] = faces.adjacent.index[iface];
        cells.faces.adjacent.ghost[jface] = faces.adjacent.ghost[iface];
        cells.faces.adjacent.rotation[jface] = faces.adjacent.rotation[iface];
    }

    cells.verts.offsets[jc] = node_beg;
    cells.verts.offsets[jc + 1] = node_beg + verts.max_count(ic);
}

void RawCells::clear() {
    resize(0, 0, 0);
}

void RawCells::resize_amr(index_t n_cells) {
    if (!adaptive_) {
        throw std::runtime_error("Resize of unstructured mesh");
    }

    index_t n_faces = n_cells * (dim_ == 2 ? 8 : 24);
    index_t n_nodes = n_cells * (dim_ == 2 ? 9 : 27);
    
    resize(n_cells, n_faces, n_nodes);
}

void RawCells::reserve_amr(index_t n_cells) {
    if (!adaptive_) {
        throw std::runtime_error("Resize of unstructured mesh");
    }

    index_t n_faces = n_cells * (dim_ == 2 ? 8 : 24);
    index_t n_nodes = n_cells * (dim_ == 2 ? 9 : 27);

    reserve(n_cells, n_faces, n_nodes);
}

void RawCells::resize_cells(index_t n_cells) {
    m_size = n_cells;

    // Поля ячеек по числу ячеек, логично
    next.resize(n_cells, -1);
    rank.resize(n_cells, -1);
    index.resize(n_cells, -1);

    center.resize(n_cells);
    volume.resize(n_cells);
    volume_alt.resize(n_cells);

    flag.resize(n_cells);
    b_idx.resize(n_cells);
    z_idx.resize(n_cells);
    level.resize(n_cells);

    // Поля данных только для ячеек
    data.resize(n_cells);

    // +1 для заключительной
    faces.offsets.resize(n_cells + 1);
    verts.offsets.resize(n_cells + 1);
}

void RawCells::reserve_cells(index_t n_cells) {
    // Поля ячеек по числу ячеек, логично
    next.reserve(n_cells);
    rank.reserve(n_cells);
    index.reserve(n_cells);

    center.reserve(n_cells);
    volume.reserve(n_cells);
    volume_alt.reserve(n_cells);

    flag.reserve(n_cells);
    b_idx.reserve(n_cells);
    z_idx.reserve(n_cells);
    level.reserve(n_cells);

    // Поля данных только для ячеек
    data.reserve(n_cells);

    // +1 для заключительной
    faces.offsets.reserve(n_cells + 1);
    verts.offsets.reserve(n_cells + 1);
}

void RawCells::shrink_to_fit_cells() {
    // Поля ячеек по числу ячеек, логично
    next.shrink_to_fit();
    rank.shrink_to_fit();
    index.shrink_to_fit();

    center.shrink_to_fit();
    volume.shrink_to_fit();
    volume_alt.shrink_to_fit();

    flag.shrink_to_fit();
    b_idx.shrink_to_fit();
    z_idx.shrink_to_fit();
    level.shrink_to_fit();

    // Поля данных только для ячеек
    data.shrink_to_fit();

    // +1 для заключительной
    faces.offsets.shrink_to_fit();
    verts.offsets.shrink_to_fit();
}

void RawCells::resize(index_t n_cells, index_t n_faces, index_t n_nodes) {
    if_debug(adaptive_) {
        z_assert(n_faces == (dim_ < 3 ? 8 : 24) * n_cells, "bad sizes");
        z_assert(n_nodes == (dim_ < 3 ? 9 : 27) * n_cells, "bad sizes");
    }
    resize_cells(n_cells);
    faces.resize(n_cells, n_faces);
    verts.resize(n_cells, n_nodes);
}

void RawCells::reserve(index_t n_cells, index_t n_faces, index_t n_nodes) {
    reserve_cells(n_cells);
    faces.reserve(n_cells, n_faces);
    verts.reserve(n_cells, n_nodes);
}

void RawCells::shrink_to_fit() {
    shrink_to_fit_cells();
    faces.shrink_to_fit();
    verts.shrink_to_fit();
}

memory_t RawCells::memory_usage() const {
    memory_t mem;
    mem.add(next);
    mem.add(rank);
    mem.add(index);
    mem.add(center);
    mem.add(volume);
    mem.add(volume_alt);
    mem.add(flag);
    mem.add(b_idx);
    mem.add(z_idx);
    mem.add(level);
    return mem;
}

void RawCells::set_cell(index_t ic, const Quad& quad) {
    assert(dim_ == 2);
    assert(adaptive_);
    assert(linear_);
    assert(!axial_);

    rank[ic] = -1;
    index[ic] = -1;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = quad.area();
    center[ic] = quad.centroid(volume[ic]);

    faces.offsets[ic] = 8 * ic;
    faces.offsets[ic + 1] = 8 * (ic + 1);
    faces.insert(faces.offsets[ic], CellType::AMR2D);

    verts.offsets[ic] = 9 * ic;
    verts.offsets[ic + 1] = 9 * (ic + 1);
    verts.mapping<2>(ic) = SqQuad(quad);

    for (index_t iface = faces.offsets[ic]; iface < faces.offsets[ic] + Side2D::count(); ++iface) {
        Line vs = {
                verts[verts.offsets[ic] + faces.vertices[iface][0]],
                verts[verts.offsets[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.center();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::INNER;
    }
}

void RawCells::set_cell(index_t ic, const Quad& quad, bool axial) {
    assert(dim_ == 2);
    assert(adaptive_);
    assert(linear_);
    assert(axial_ == axial);

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

    faces.offsets[ic] = 8 * ic;
    faces.offsets[ic + 1] = 8 * (ic + 1);
    faces.insert(faces.offsets[ic], CellType::AMR2D);

    verts.offsets[ic] = 9 * ic;
    verts.offsets[ic + 1] = 9 * (ic + 1);
    verts.mapping<2>(ic) = SqQuad(quad);

    for (auto i: verts.range(ic)) {
        if (verts[i].z() != 0.0) {
            throw std::runtime_error("RawCells add axial cell, vertex.z != 0.0");
        }
    }

    for (index_t iface = faces.offsets[ic]; iface < faces.offsets[ic] + Side2D::count(); ++iface) {
        Line vs = {
                verts[verts.offsets[ic] + faces.vertices[iface][0]],
                verts[verts.offsets[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.centroid(axial);
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::INNER;

        if (axial) {
            faces.area_alt[iface] = vs.area_as();
        }
    }
}

void RawCells::set_cell(index_t ic, const SqQuad& quad) {
    set_cell(ic, quad.reduce());
    //std::cerr << "Nonlinear RawCells is not supported\n";
}

void RawCells::set_cell(index_t ic, const SqQuad& quad, bool axial) {
    set_cell(ic, quad.reduce(), axial);
    //std::cerr << "Nonlinear RawCells is not supported\n";
}

void RawCells::set_cell(index_t ic, const Cube& cube) {
    assert(dim_ == 3);
    assert(adaptive_);
    assert(linear_);
    assert(!axial_);

    rank[ic] = -1;
    index[ic] = -1;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = cube.volume();
    center[ic] = cube.centroid(volume[ic]);

    faces.offsets[ic] = 24 * ic ;
    faces.offsets[ic + 1] = 24 * (ic + 1);
    faces.insert(faces.offsets[ic], CellType::AMR3D);

    verts.offsets[ic] = 27 * ic;
    verts.offsets[ic + 1] = 27 * (ic + 1);
    verts.mapping<3>(ic) = SqCube(cube);

    for (index_t iface = faces.offsets[ic]; iface < faces.offsets[ic] + Side3D::count(); ++iface) {
        Quad vs = {
                verts[verts.offsets[ic] + faces.vertices[iface][0]],
                verts[verts.offsets[ic] + faces.vertices[iface][1]],
                verts[verts.offsets[ic] + faces.vertices[iface][2]],
                verts[verts.offsets[ic] + faces.vertices[iface][3]]
        };

        faces.area[iface]     = vs.area();
        faces.center[iface]   = vs.center();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::INNER;
    }
}

void RawCells::set_cell(index_t ic, const SqCube& cube) {
    set_cell(ic, cube.reduce());
}

void RawCells::push_back(const geom::Line &line) {
    throw std::runtime_error("NO WAY ACPBWER");
}

void RawCells::push_back(const Polygon& poly) {
    assert(dim_ == 2);
    assert(!adaptive_);
    assert(linear_);

    index_t ic = n_cells();

    resize_cells(ic + 1);

    rank[ic] = 0;
    index[ic] = ic;

    flag[ic] = 0;
    b_idx[ic] = -1;
    z_idx[ic] = -1;
    level[ic] = -1;

    volume[ic] = poly.area();
    center[ic] = poly.centroid(volume[ic]);

    if (axial_) {
        volume_alt[ic] = poly.volume_as();
    }

    int n_nodes = poly.size();
    int n_faces = poly.size();

    faces.resize(ic + 1, faces.n_faces() + n_faces);
    verts.resize(ic + 1, verts.n_verts() + n_nodes);

    faces.offsets[ic + 1] = faces.offsets[ic] + n_faces;
    faces.insert(faces.offsets[ic], CellType::POLYGON, n_faces);

    verts.offsets[ic + 1] = verts.offsets[ic] + n_nodes;
    for (int i = 0; i < n_nodes; ++i) {
        verts[verts.offsets[ic] + i] = poly[i];
    }

    for (index_t iface = faces.offsets[ic]; iface < faces.offsets[ic] + n_faces; ++iface) {
        Line vs = {
                verts[verts.offsets[ic] + faces.vertices[iface][0]],
                verts[verts.offsets[ic] + faces.vertices[iface][1]]
        };

        faces.area[iface]     = vs.length();
        faces.center[iface]   = vs.centroid();
        faces.normal[iface]   = vs.normal(center[ic]);
        faces.boundary[iface] = Boundary::INNER;

        if (axial_) {
            faces.area_alt[iface] = vs.area_as();
        }
    }
}

void RawCells::push_back(const Polyhedron& poly) {
    if (poly.need_simplify(RawFaces::max_vertices)) {
        Polyhedron simple_poly = poly;
        simple_poly.simplify_faces(RawFaces::max_vertices);
        push_back_impl(simple_poly);
    }
    else {
        push_back_impl(poly);
    }
}

void RawCells::push_back_impl(const Polyhedron& poly) {
    assert(dim_ == 3);
    assert(!adaptive_);
    assert(linear_);
    assert(!axial_);

    index_t ic = n_cells();

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
    verts.resize(ic + 1, verts.n_verts() + n_nodes);

    verts.offsets[ic + 1] = verts.offsets[ic] + n_nodes;
    for (int i = 0; i < n_nodes; ++i) {
        verts[verts.offsets[ic] + i] = poly.vertex(i);
    }

    // Определим грани многогранника
    int n_faces = poly.n_faces();
    faces.resize(ic + 1, faces.n_faces() + n_faces);

    faces.offsets[ic + 1] = faces.offsets[ic] + n_faces;
    for (int i = 0; i < poly.n_faces(); ++i) {
        index_t iface = faces.offsets[ic] + i;

        faces.area[iface] = poly.face_area(i);
        faces.center[iface] = poly.face_center(i);
        faces.normal[iface] = poly.face_normal(i);
        faces.boundary[iface] = Boundary::INNER;

        faces.area_alt[iface] = NAN;

        // Выставить вершины грани
        int n_verts = poly.face_indices(i).size();
        if (n_verts > RawFaces::max_vertices) {
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

namespace fs = std::filesystem;

template <typename T>
void save_vector(
    const fs::path& root,
    std::ofstream& file,
    const fs::path& path,
    const std::string& tab,
    const std::string& name,
    const std::vector<T>& arr,
    const std::string& end=""
) {
    // Сохранить бинарные данные
    fs::path fullname = root / path / (name + (mpi::single() ? "" : (".pt" + mpi::srank())) + ".bin");
    std::ofstream binary(fullname, std::ios::binary | std::ios::trunc | std::ios::out);
    binary.write(reinterpret_cast<const char*>(arr.data()), arr.size() * sizeof(T));
    binary.close();

    // Собрать размеры со всех процессов
    std::vector<size_t> sizes;
    if (!mpi::single()) {
        mpi::all_gather(arr.size(), sizes);
    }

    // Все побочные процессы завершают функцию
    if (!mpi::master()) return;

    file << tab << "\"" << name << "\": {\n";
    file << tab << "  \"sizeof\": " << sizeof(T) << ",\n";

    // Относительное имя файла
    std::string filename = (path / name).string();
    if (mpi::single()) {
        file << tab << "  \"data\": ";
        file << "{ \"size\": " << arr.size() << ", \"file\": \"" << filename << ".bin\" }\n";
    }
    else {
        file << tab << "  \"data\": [\n";
        for (int i = 0; i < sizes.size() - 1; ++i) {
            file << tab << "    { \"size\": " << sizes[i] << ", \"file\": \"" << filename << ".pt" << i << ".bin\" },\n";
        }
        file << tab << "    { \"size\": " << sizes.back() << ", \"file\": \"" << filename << ".pt" << sizes.size() << ".bin\" }\n";
        file << tab << "  ]\n";
    }
    file << tab << "}" << end;
}

inline void save_buffer(
    const fs::path& root,
    std::ofstream& file,
    const fs::path& path,
    const std::string& tab,
    const std::string& var,
    const Storage& storage,
    const std::string& end=""
) {
    // Получить ссылку на буфер
    if (!storage.contain(var)) {
        throw std::runtime_error("No such variable \"" + var + "\" in Storage");
    }
    const utils::Buffer& buffer = storage[var];

    // Сохранить бинарные данные
    fs::path fullname = root / path / (var + (mpi::single() ? "" : (".pt" + mpi::srank())) + ".bin");
    std::ofstream binary(fullname, std::ios::binary | std::ios::trunc | std::ios::out);
    binary.write(reinterpret_cast<const char*>(buffer.data()), buffer.byte_size());
    binary.close();

    // Собрать размеры со всех процессов
    std::vector<size_t> sizes;
    if (!mpi::single()) {
        mpi::all_gather(buffer.size(), sizes);
    }

    // Все побочные процессы завершают функцию
    if (!mpi::master()) return;

    file << tab << "\"" << var << "\": {\n";
    file << tab << "  \"count\": " << buffer.count() << ",\n";
    file << tab << "  \"sizeof\": " << buffer.element_size() << ",\n";

    // Относительное имя файла
    std::string filename = (path / var).string();
    if (mpi::single()) {
        file << tab << "  \"data\": ";
        file << "{ \"size\": " << buffer.size() << ", \"file\": \"" << filename << ".bin\" }\n";
    }
    else {
        file << tab << "  \"data\": [\n";
        for (int i = 0; i < sizes.size() - 1; ++i) {
            file << tab << "    { \"size\": " << sizes[i] << ", \"file\": \"" << filename << ".pt" << i << ".bin\" },\n";
        }
        file << tab << "    { \"size\": " << sizes.back() << ", \"file\": \"" << filename << ".pt" << sizes.size() << ".bin\" }\n";
        file << tab << "  ]\n";
    }
    file << tab << "}" << end;
}

void RawCells::backup(const std::filesystem::path& root, std::ofstream& file,
    const std::string& tab, const std::vector<std::string>& variables) const {

    if (!fs::exists(root) && !fs::is_directory(root)) {
        throw std::runtime_error("RawCells::backup(): directory " + root.string() + "\" doesn't exist");
    }

    if (mpi::master() && !file.is_open()) {
        throw std::runtime_error("RawCells::backup(): cannot open output file.");
    }

    const std::string tab2 = tab + "  ";
    const std::string tab3 = tab + "    ";

    if (mpi::master()) {
        file << std::boolalpha;
        file << tab << "\"dim\":      " << dim_ << ",\n";
        file << tab << "\"adaptive\": " << adaptive_ << ",\n";
        file << tab << "\"axial\":    " << axial_ << ",\n";
        file << tab << "\"linear\":   " << linear_ << ",\n";
    }

    if (mpi::single()) {
        file << tab << "\"size\":     " << n_cells() << ",\n";
    }
    else {
        auto sizes = mpi::all_gather(n_cells());
        if (mpi::master()) {
            file << tab << "\"sizes\": [";
            for (int i = 0; i < sizes.size() - 1; ++i) {
                file << sizes[i] << ", ";
            }
            file << sizes.back() << "],\n";
        }
    }

    save_vector(root, file, "cells", tab, "rank",  rank, ",\n");
    save_vector(root, file, "cells", tab, "index", index, ",\n");
    save_vector(root, file, "cells", tab, "level", level, ",\n");
    save_vector(root, file, "cells", tab, "b_idx", b_idx, ",\n");
    save_vector(root, file, "cells", tab, "z_idx", z_idx, ",\n");
    save_vector(root, file, "cells", tab, "center", center, ",\n");
    save_vector(root, file, "cells", tab, "volume", volume, ",\n");
    if (!volume_alt.empty()) {
        save_vector(root, file, "cells", tab, "volume_alt", volume_alt, ",\n");
    }

    if (mpi::master()) {
        fs::create_directory(root / "cells/faces");
        file << tab << "\"faces\": {\n";
    }

    save_vector(root, file, "cells/faces", tab2, "offsets", faces.offsets, ",\n");

    if (mpi::master()) {
        fs::create_directory(root / "cells/faces/adjacent");
        file << tab2 << "\"adjacent\": {\n";
    }

    save_vector(root, file, "cells/faces/adjacent", tab3, "rank", faces.adjacent.rank, ",\n");
    save_vector(root, file, "cells/faces/adjacent", tab3, "index", faces.adjacent.index, ",\n");
    save_vector(root, file, "cells/faces/adjacent", tab3, "ghost", faces.adjacent.ghost, ",\n");
    save_vector(root, file, "cells/faces/adjacent", tab3, "basic", faces.adjacent.basic, "\n");

    if (mpi::master()) {
        file << tab2 << "},\n"; // mesh.cells.faces.adjacent
    }

    save_vector(root, file, "cells/faces", tab2, "boundary", faces.boundary, ",\n");
    save_vector(root, file, "cells/faces", tab2, "normal", faces.normal, ",\n");
    save_vector(root, file, "cells/faces", tab2, "center", faces.center, ",\n");
    save_vector(root, file, "cells/faces", tab2, "area", faces.area, ",\n");
    if (!faces.area_alt.empty()) {
        save_vector(root, file, "cells/faces", tab2, "area_alt", faces.area_alt, ",\n");
    }
    save_vector(root, file, "cells/faces", tab2, "vertices", faces.vertices, "\n");

    if (mpi::master()) {
        file << tab << "},\n"; // mesh.cells.faces
    }

    throw std::runtime_error("Backup nodes error");
    // Старая версия сохранения узлов, нужна новая
    //save_vector(root, file, "cells", tab, "verts", verts, ",\n");

    if (mpi::master()) {
        fs::create_directory(root / "cells/data");
        file << tab << "\"data\": {\n";
    }

    for (int i = 0; i < variables.size() - 1; ++i) {
        save_buffer(root, file, "cells/data", tab2, variables[i], data, ",\n");
    }
    save_buffer(root, file, "cells/data", tab2, variables.back(), data, "\n");

    if (mpi::master()) {
        file << tab << "}\n"; // mesh.cells.data
    }
}

} // namespace zephyr::mesh