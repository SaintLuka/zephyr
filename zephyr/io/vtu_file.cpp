#include <fstream>
#include <filesystem>

#include <zephyr/geom/side.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_node.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr::io {

using namespace zephyr::geom;
using namespace zephyr::mesh;

namespace {

class DistNodes {
    // Могут быть указатели на реальные узлы
    AmrNodes* local_nodes_ = nullptr;
    AmrNodes* ghost_nodes_ = nullptr;

    // Или буферы данных
    std::optional<std::vector<mesh::index_t>> vert_index_buffer_ = std::nullopt;
    std::optional<std::vector<Vector3d>> local_coords_buffer_ = std::nullopt;

    bool has_no_ghosts_{true};
    mesh::index_t n_nodes_{0};
    mesh::index_t n_local_nodes_{0};
    mesh::index_t n_ghost_nodes_{0};

    std::span<const mesh::index_t> vert_index_{}; // по классике verts.index
    std::span<const mesh::index_t> vert_ghost_{}; // по классике verts.ghost
    std::span<const Vector3d> local_coords_{};
    std::span<const Vector3d> ghost_coords_{};

public:
    DistNodes() = default;

    // Конструктор из реальных узлов
    DistNodes(const AmrCells& cells, AmrNodes& local_nodes, AmrNodes& ghost_nodes)
            : local_nodes_(&local_nodes), ghost_nodes_(&ghost_nodes) {

        vert_index_ = cells.verts.index;
        vert_ghost_ = cells.verts.ghost;
        local_coords_ = local_nodes_->coord;
        ghost_coords_ = ghost_nodes_->coord;

        n_local_nodes_ = local_coords_.size();
        n_ghost_nodes_ = ghost_coords_.size();
        n_nodes_ = n_local_nodes_ + n_ghost_nodes_;
        has_no_ghosts_ = vert_ghost_.empty() || ghost_coords_.empty();

        z_assert(!vert_index_.empty(), "DistNodes empty vert_index #1");
    }

    // Конструктор из вспомогательных массивов
    DistNodes(std::vector<mesh::index_t>&& vert_index, std::vector<Vector3d>&& local_coords)
        : vert_index_buffer_(vert_index), local_coords_buffer_(local_coords) {

        vert_index_ = *vert_index_buffer_;
        vert_ghost_ = {};
        local_coords_ = *local_coords_buffer_;
        ghost_coords_ = {};

        n_local_nodes_ = local_coords_.size();
        n_ghost_nodes_ = ghost_coords_.size();
        n_nodes_ = n_local_nodes_ + n_ghost_nodes_;
        has_no_ghosts_ = vert_ghost_.empty() || ghost_coords_.empty();

        z_assert(!vert_index_.empty(), "DistNodes empty vert_index #2");
    }

    bool empty() const { return n_nodes_ == 0; }

    bool true_nodes() const { return local_nodes_ != nullptr && ghost_nodes_ != nullptr; }

    AmrNodes& local_nodes() const { return *local_nodes_; }

    AmrNodes& ghost_nodes() const { return *ghost_nodes_; }

    index_t n_nodes() const { return n_nodes_; }

    index_t n_local_nodes() const { return n_local_nodes_; }

    index_t n_ghost_nodes() const { return n_ghost_nodes_; }

    Vector3d coord(index_t vert_idx) const {
        if (has_no_ghosts_ || vert_ghost_[vert_idx] < 0) {
            return local_coords_[vert_index_[vert_idx]];
        }
        else {
            return ghost_coords_[vert_ghost_[vert_idx]];
        }
    }

    index_t index(index_t vert_idx) const {
        if (has_no_ghosts_ || vert_ghost_[vert_idx] < 0) {
            return vert_index_[vert_idx];
        }
        else {
            return n_local_nodes_ + vert_ghost_[vert_idx];
        }
    }

    std::vector<Vector3d> points() const {
        std::vector<Vector3d> res;
        res.reserve(n_nodes_);
        std::ranges::copy(local_coords_, std::back_inserter(res));
        std::ranges::copy(ghost_coords_, std::back_inserter(res));
        return res;
    }
};

// ====================================================================================================================
//                                                        VTU FORMAT
// ====================================================================================================================

// Полное описание геометрии неструктурированной сетки в формате VTU
class VtuStructure {
public:
    VtuStructure(const AmrCells &cells, bool polyhedral);

    VtuStructure(const AmrCells &cells, const DistNodes& nodes, bool polyhedral);

    void write_header(std::ofstream &file, bool node_data, const Variables &variables) const;

    void write_primitives(std::ofstream &file) const;

private:
    void fill(const AmrCells& cells, const DistNodes &nodes, bool polyhedral);

    // Двумерная адаптивная сетка в виде простых квадратов
    void fill_adaptive_hex_2D(const AmrCells& cells, const DistNodes &nodes = {});

    // Двумерная адаптивная сетка в виде полигонов
    void fill_adaptive_poly_2D(const AmrCells& cells, const DistNodes& nodes = {});

    // Трёхмерная адаптивная сетка в виде простых шестигранников
    void fill_adaptive_hex_3D(const AmrCells& cells, const DistNodes &nodes = {});

    // Двумерная полигональная сетка или трёхмерная сетка из базовых примитивов
    void fill_poly_classic(const AmrCells& cells, const DistNodes &nodes = {});

    // Трёхмерная сетка из нестандартных многогранников
    void fill_polyfaces_3D(const AmrCells& cells, const DistNodes &nodes = {});

    // Массивы данных для VTU
    std::vector<Vector3d> points;
    std::vector<index_t>  connectivity;
    std::vector<index_t>  offsets;
    std::vector<type_t>   types;

    // Для сеток из многогранников
    std::vector<index_t> faces;
    std::vector<index_t> face_offsets;
};

std::vector<index_t> arange(index_t size) {
    std::vector<index_t> arr(size);
    for (index_t i = 0; i < size; ++i) arr[i] = i;
    return arr;
}

std::vector<index_t> uniform_offsets(index_t size, index_t step) {
    std::vector<index_t> offsets(size);
    for (index_t i = 0; i < size; ++i) {
        offsets[i] = step * (i + 1);
    }
    return offsets;
}

type_t nv_to_type2D(int nv) {
    switch (nv) {
        case 3:  return VTK_TRIANGLE;
        case 4:  return VTK_QUAD;
        default: return VTK_POLYGON;
    }
}

type_t nv_to_type3D(int nv) {
    switch (nv) {
        case 4:  return VTK_TETRA;
        case 5:  return VTK_PYRAMID;
        case 6:  return VTK_WEDGE;
        default: return VTK_HEXAHEDRON;
    }
}

std::vector<Vector3d> collect_points(const AmrCells& cells, index_t n_points) {
    std::vector<Vector3d> points(n_points);
    index_t iv = 0;
    for (mesh::index_t ic = 0; ic < cells.n_cells(); ++ic) {
        const Vector3d* vertices = cells.verts.coords_data(ic);

        int nv = cells.verts.count(ic);
        for (int j = 0; j < nv; ++j) {
            points[iv++] = vertices[j];
        }
    }
    return points;
}

VtuStructure::VtuStructure(const AmrCells &cells, bool polyhedral) {
    static constexpr DistNodes nodes;
    fill(cells, nodes, polyhedral);
}

VtuStructure::VtuStructure(const AmrCells &cells, const DistNodes& nodes, bool polyhedral) {
    fill(cells, nodes, polyhedral);
}

void VtuStructure::fill(const AmrCells& cells, const DistNodes &nodes, bool polyhedral) {
    if (cells.adaptive()) {
        if (cells.dim() < 3) {
            if (!polyhedral) {
                fill_adaptive_hex_2D(cells, nodes);
                return;
            }
            else {
                fill_adaptive_poly_2D(cells, nodes);
                return;
            }
        }
        else {
            fill_adaptive_hex_3D(cells, nodes);
            return;
        }
    }
    else {
        if (cells.dim() == 2) {
            fill_poly_classic(cells, nodes);
            return;
        }
        else {
            if (!polyhedral) {
                fill_poly_classic(cells, nodes);
                return;
            }
            else {
                fill_polyfaces_3D(cells, nodes);
                return;
            }
        }
    }
}

void VtuStructure::fill_adaptive_hex_2D(const AmrCells& cells, const DistNodes &nodes) {
    index_t n_cells = cells.n_cells();
    index_t n_points = 4 * cells.n_cells();

    types.resize(n_cells, VTK_QUAD);
    if (nodes.empty()) {
        points.resize(n_points);
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            const SqQuad& vertices = cells.verts.mapping<2>(ic);
            points[4 * ic + 0] = vertices.vs<-1, -1>();
            points[4 * ic + 1] = vertices.vs<+1, -1>();
            points[4 * ic + 2] = vertices.vs<+1, +1>();
            points[4 * ic + 3] = vertices.vs<-1, +1>();
        }
        connectivity = arange(n_points);
    }
    else {
        points = nodes.points();
        connectivity.resize(n_points);
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            index_t beg = cells.verts.offsets[ic];
            connectivity[4 * ic + 0] = nodes.index(beg + SqQuad::iss<-1, -1>());
            connectivity[4 * ic + 1] = nodes.index(beg + SqQuad::iss<+1, -1>());
            connectivity[4 * ic + 2] = nodes.index(beg + SqQuad::iss<+1, +1>());
            connectivity[4 * ic + 3] = nodes.index(beg + SqQuad::iss<-1, +1>());
        }
    }
    offsets = uniform_offsets(n_cells, 4);
}

void VtuStructure::fill_adaptive_poly_2D(const AmrCells& cells, const DistNodes& nodes) {
    index_t n_cells = cells.n_cells();
    types.resize(n_cells);
    offsets.resize(n_cells);

    index_t n_points = 0;
    for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
        int nv = 4;
        if (cells.faces.is_complex(ic, Side2D::L)) { nv += 1; }
        if (cells.faces.is_complex(ic, Side2D::R)) { nv += 1; }
        if (cells.faces.is_complex(ic, Side2D::B)) { nv += 1; }
        if (cells.faces.is_complex(ic, Side2D::T)) { nv += 1; }

        n_points += nv;
        offsets[ic] = n_points;
        types[ic] = nv_to_type2D(nv);
    }

    if (nodes.empty()) {
        points.resize(n_points);
        index_t iv = 0;
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            const SqQuad& vertices = cells.verts.mapping<2>(ic);

            points[iv++] = vertices.vs<-1, -1>();
            if (cells.faces.is_complex(ic, Side2D::B)) { points[iv++] = vertices.vs<0, -1>(); }

            points[iv++] = vertices.vs<+1, -1>();
            if (cells.faces.is_complex(ic, Side2D::R)) { points[iv++] = vertices.vs<+1, 0>(); }

            points[iv++] = vertices.vs<+1, +1>();
            if (cells.faces.is_complex(ic, Side2D::T)) { points[iv++] = vertices.vs<0, +1>(); }

            points[iv++] = vertices.vs<-1, +1>();
            if (cells.faces.is_complex(ic, Side2D::L)) { points[iv++] = vertices.vs<-1, 0>(); }
        }

        connectivity = arange(n_points);
    }
    else {
        points = nodes.points();

        index_t iv = 0;
        connectivity.resize(n_points);
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            index_t beg = cells.verts.offsets[ic];

            connectivity[iv++] = nodes.index(beg + SqQuad::iss<-1, -1>());
            if (cells.faces.is_complex(ic, Side2D::B)) { connectivity[iv++] = nodes.index(beg + SqQuad::iss<0, -1>()); }

            connectivity[iv++] = nodes.index(beg + SqQuad::iss<+1, -1>());
            if (cells.faces.is_complex(ic, Side2D::R)) { connectivity[iv++] = nodes.index(beg + SqQuad::iss<+1, 0>()); }

            connectivity[iv++] = nodes.index(beg + SqQuad::iss<+1, +1>());
            if (cells.faces.is_complex(ic, Side2D::T)) { connectivity[iv++] = nodes.index(beg + SqQuad::iss<0, +1>()); }

            connectivity[iv++] = nodes.index(beg + SqQuad::iss<-1, +1>());
            if (cells.faces.is_complex(ic, Side2D::L)) { connectivity[iv++] = nodes.index(beg + SqQuad::iss<-1, 0>()); }
        }
    }
}

void VtuStructure::fill_adaptive_hex_3D(const AmrCells& cells, const DistNodes &nodes) {
    index_t n_cells = cells.n_cells();
    index_t n_points = 8 * cells.n_cells();

    types.resize(n_cells, VTK_HEXAHEDRON);

    if (nodes.empty()) {
        points.resize(n_points);
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            const SqCube& vertices = cells.verts.mapping<3>(ic);
            points[8 * ic + 0] = vertices.vs<-1, -1, -1>();
            points[8 * ic + 1] = vertices.vs<+1, -1, -1>();
            points[8 * ic + 2] = vertices.vs<+1, +1, -1>();
            points[8 * ic + 3] = vertices.vs<-1, +1, -1>();
            points[8 * ic + 4] = vertices.vs<-1, -1, +1>();
            points[8 * ic + 5] = vertices.vs<+1, -1, +1>();
            points[8 * ic + 6] = vertices.vs<+1, +1, +1>();
            points[8 * ic + 7] = vertices.vs<-1, +1, +1>();
        }
        connectivity = arange(n_points);
    }
    else {
        points = nodes.points();
        connectivity.resize(n_points);
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            index_t beg = cells.verts.offsets[ic];
            connectivity[8 * ic + 0] = nodes.index(beg + SqCube::iss<-1, -1, -1>());
            connectivity[8 * ic + 1] = nodes.index(beg + SqCube::iss<+1, -1, -1>());
            connectivity[8 * ic + 2] = nodes.index(beg + SqCube::iss<+1, +1, -1>());
            connectivity[8 * ic + 3] = nodes.index(beg + SqCube::iss<-1, +1, -1>());
            connectivity[8 * ic + 4] = nodes.index(beg + SqCube::iss<-1, -1, +1>());
            connectivity[8 * ic + 5] = nodes.index(beg + SqCube::iss<+1, -1, +1>());
            connectivity[8 * ic + 6] = nodes.index(beg + SqCube::iss<+1, +1, +1>());
            connectivity[8 * ic + 7] = nodes.index(beg + SqCube::iss<-1, +1, +1>());
        }
    }
    offsets = uniform_offsets(n_cells, 8);
}

void VtuStructure::fill_poly_classic(const AmrCells& cells, const DistNodes &nodes) {
    index_t n_cells = cells.n_cells();

    auto nv_to_type = cells.dim() == 2 ? nv_to_type2D : nv_to_type3D;

    index_t n_points = 0;
    types.resize(n_cells);
    offsets.resize(n_cells);
    for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
        int nv = cells.verts.count(ic);
        n_points += nv;
        offsets[ic] = n_points;
        types[ic] = nv_to_type(nv);
    }

    if (nodes.empty()) {
        // Вершины с дублированием
        points = collect_points(cells, n_points);
        connectivity = arange(n_points);
    }
    else {
        points = nodes.points();
        connectivity.resize(n_points);
        index_t offset = 0;
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            index_t beg = cells.verts.offsets[ic];
            for (int j = 0; j < offsets[ic] - offset; ++j) {
                connectivity[offset + j] = nodes.index(beg + j);
            }
            offset = offsets[ic];
        }
    }
}

void VtuStructure::fill_polyfaces_3D(const AmrCells& cells, const DistNodes &nodes) {
    // Данные многогранников
    index_t n_cells = cells.n_cells();

    index_t n_points = 0;
    types.resize(n_cells);
    offsets.resize(n_cells);

    for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
        int nv = cells.verts.count(ic);
        n_points += nv;
        offsets[ic] = n_points;
        types[ic] = VTK_POLYHEDRON;
    }

    if (nodes.empty()) {
        // Вершины с дублированием
        points = collect_points(cells, n_points);
        connectivity = arange(n_points);
    }
    else {
        points = nodes.points();
        connectivity.resize(n_points);
        index_t offset = 0;
        for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
            index_t beg = cells.verts.offsets[ic];
            for (int j = 0; j < offsets[ic] - offset; ++j) {
                connectivity[offset + j] = nodes.index(beg + j);
            }
            offset = offsets[ic];
        }
    }

    // Данные многогранников

    // Faces
    bool unique = !nodes.empty();
    for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
        faces.push_back(cells.face_count(ic));

        // Массив для описания грани
        for (auto iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) continue;

            int nv = cells.faces.n_vertices(iface);
            faces.push_back(nv);
            for (int j = 0; j < nv; ++j) {
                auto loc_i = cells.faces.vertices[iface][j];
                index_t node_idx = cells.verts.offsets[ic] + loc_i;
                if (unique) {
                    node_idx = nodes.index(node_idx);
                }
                faces.push_back(node_idx);
            }
        }
    }

    // FaceOffsets
    index_t offset = 0;
    for (mesh::index_t ic = 0; ic < n_cells; ++ic) {
        int n_fverts = 0;
        for (auto iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface)) continue;

            // Допускаются грани с числом вершин до 8
            n_fverts += cells.faces.n_vertices(iface) + 1;
        }
        offset += n_fverts + 1;

        face_offsets.push_back(offset);
    }
}

void VtuStructure::write_header(std::ofstream &file, bool node_data, const Variables &variables) const {
    index_t n_cells = types.size();
    index_t n_nodes = points.size();

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    offset_t byte_offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    file << "      </Points>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(points);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(connectivity);

    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(offsets);

    file << "        <DataArray type=\"" << VtkType::get<type_t>() << "\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(types);

    if (!faces.empty()) {
        file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"faces\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + buffer_size(faces);

        file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"faceoffsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + buffer_size(face_offsets);
    }

    file << "      </Cells>\n";

    // CellData
    file << "      <CellData>\n";
    for (auto &field: variables.list()) {
        if (!field.cell_data()) continue;

        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";

        byte_offset += sizeof(datasize_t) + n_cells * field.size();
    }
    file << "      </CellData>\n";

    // PointData
    if (node_data) {
        file << "      <PointData>\n";
        for (auto &field: variables.list()) {
            if (!field.node_data()) continue;

            file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
            if (!field.is_scalar()) {
                file << "\" NumberOfComponents=\"" << field.n_components();
            }
            file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";

            byte_offset += sizeof(datasize_t) + n_nodes * field.size();
        }
        file << "      </PointData>\n";
    }

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void VtuStructure::write_primitives(std::ofstream &file) const {
    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    write_buffer(file, points);       // PointsCoords
    write_buffer(file, connectivity); // Cells.Connectivity
    write_buffer(file, offsets);      // Cells.Offsets
    write_buffer(file, types);        // Cells.Types

    // Данные многогранников
    if (!faces.empty()) {
        write_buffer(file, faces);         // Cells.Faces
        write_buffer(file, face_offsets);  // Cells.FaceOffsets
    }
}

void write_cells_data(std::ofstream &file, AmrCells &cells, const Variables &variables) {
    std::vector<char> temp;

    const index_t n_cells = cells.n_cells();

    for (auto &field: variables.list()) {
        if (!field.cell_data()) continue;

        index_t field_size = field.size();
        datasize_t data_size = n_cells * field_size;

        file.write((char*) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        index_t counter = 0;
        for (auto& cell: cells) {
            field.write(cell, temp.data() + counter * field_size);
            ++counter;
        }

        file.write(temp.data(), data_size);
    }
}

void write_nodes_data(std::ofstream &file, const DistNodes &nodes, const Variables &variables) {
    if (!nodes.true_nodes()) return;

    std::vector<char> temp;

    const index_t n_nodes = nodes.n_nodes();

    for (auto &field: variables.list()) {
        if (!field.node_data()) continue;

        index_t field_size = field.size();
        datasize_t data_size = n_nodes * field_size;

        file.write((char*) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        index_t counter = 0;
        for (index_t i = 0; i < nodes.n_local_nodes(); ++i) {
            EuNode node(&nodes.local_nodes(), i);
            field.write(node, temp.data() + counter * field_size);
            ++counter;
        }
        for (index_t i = 0; i < nodes.n_ghost_nodes(); ++i) {
            EuNode node(&nodes.ghost_nodes(), i);
            field.write(node, temp.data() + counter * field_size);
            ++counter;
        }

        file.write(temp.data(), data_size);
    }
}

} // anonymous namespace

// ====================================================================================================================
//                                                       VTU FILE
// ====================================================================================================================

VtuFile::VtuFile(std::string_view filename, const Variables &variables, const VtuOptions& options) :
    filename(filename), variables(variables), options(options) {
}

void VtuFile::save(EuMesh &mesh) const {
    save(filename, mesh, variables, options);
}

void VtuFile::save(AmrCells &cells) const {
    save(filename, cells, variables, options);
}

void save_with_nodes(std::string_view filename, AmrCells &locals,
    const DistNodes& nodes, const Variables &variables, const VtuOptions& options) {
    std::string fullname = add_extension(filename, ".vtu");
    create_directories(fullname);

    std::ofstream file(fullname, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << fullname << "'\n";
        return;
    }

    VtuStructure formatter(locals, nodes, options.polyhedral);
    formatter.write_header(file, nodes.true_nodes(), variables);
    formatter.write_primitives(file);
    write_cells_data(file, locals, variables);
    write_nodes_data(file, nodes, variables);

    file.close();
}

void VtuFile::save(std::string_view filename, AmrCells &locals,
                   const Variables &variables, const VtuOptions& options) {

    DistNodes nodes;
    if (options.unique_nodes) {
        auto [inc_verts, inc_nodes] = AmrNodes::generate<false>(locals);
        nodes = DistNodes(
            std::move(inc_verts.index),
            std::move(inc_nodes.coord));
    }
    save_with_nodes(filename, locals, nodes, variables, options);
}

void VtuFile::save(std::string_view filename, EuMesh &mesh,
                   const Variables &variables, const VtuOptions& options) {

    DistNodes nodes;
    if (options.unique_nodes) {
        if (mesh.has_nodes()) {
            nodes = DistNodes(
                mesh.local_cells(),
                mesh.local_nodes(),
                mesh.ghost_nodes());
        }
        else {
            auto [inc_verts, inc_nodes] = AmrNodes::generate<false>(mesh.local_cells());
            nodes = DistNodes(
                std::move(inc_verts.index),
                std::move(inc_nodes.coord));
        }
    }
    save_with_nodes(filename, mesh.local_cells(), nodes, variables, options);
}

} // namespace zephyr::io