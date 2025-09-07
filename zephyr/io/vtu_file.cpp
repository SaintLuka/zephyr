#include <fstream>
#include <filesystem>

#include <zephyr/io/vtu_file.h>
#include <zephyr/mesh/side.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::io {

using namespace zephyr::geom;
using namespace zephyr::mesh;

namespace fs = std::filesystem;

/// ===========================================================================
///             Реализация записи в виде набора статических функций
/// ===========================================================================

inline bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

inline std::string byteorder() {
    return is_big_endian() ? "BigEndian" : "LittleEndian";
}

namespace {
// Тип для задания размера массивов в бинарном файле, я не уверен, что его можно менять
using datasize_t = std::uint32_t;
using offset_t   = std::uint32_t;  // Тип смещения массивов в байтах
using index_t    = std::uint32_t;  // Тип нумерации примитивов (ячеек, вершин и т.д.)
using type_t     = std::uint8_t;   // VTK тип примитива/ячейки

using byte_ptr = char*;
}

// Обработчик ячейки
struct Handler {
    bool hex_only, polyhedral;

    explicit Handler(bool hex_only = false, bool polyhedral = false)
        : hex_only(hex_only), polyhedral(polyhedral) {
    }

    // VTK тип ячейки
    type_t type(const EuCell& cell) const {
        if (cell.adaptive()) {
            // Адаптивная ячейка
            if (cell.dim() < 3) {
                if (hex_only) {
                    return 9;  // VTK_QUAD
                } else {
                    if (cell.complex_face(Side2D::L) || cell.complex_face(Side2D::R) ||
                        cell.complex_face(Side2D::B) || cell.complex_face(Side2D::T)) {
                        return 7;  // VTK_POLYGON
                    } else {
                        return 9;  // VTK_QUAD
                    }
                }
            } else {
                return 12;  // VTK_HEXAHEDRON

                // В будущем использовать TRIQUADRIC HEXAHEDRON
                // return 29; // VTK_TRIQUADRIC_HEXAHEDRON
            }
        }
        else {
            // Обычная эйлерова ячейка

            if (cell.dim() < 3) {
                // Двумерный полигон
                switch (cell.node_count()) {
                    case 3:  return 5;  // VTK_TRIANGLE
                    case 4:  return 9;  // VTK_QUAD
                    default: return 7;  // VTK_POLYGON
                }
            } else {
                // Многогранник общего вида
                if (polyhedral) {
                    return 42;  // VTK_POLYHEDRON
                }

                // Один из доступных примитивов
                switch (cell.node_count()) {
                    case 4:  return 10;  // VTK_TETRA
                    case 5:  return 14;  // VTK_PYRAMID
                    case 6:  return 13;  // VTK_WEDGE
                    default: return 12;  // VTK_HEXAHEDRON
                }
            }
        }
    }

    // Количество вершин элемента
    index_t n_points(const EuCell& cell) const {
        if (cell.adaptive()) {
            // Адаптивная ячейка
            if (cell.dim() < 3) {
                if (hex_only) {
                    return 4;
                } else {
                    index_t n = 4;
                    if (cell.complex_face(Side2D::LEFT))   { n += 1; }
                    if (cell.complex_face(Side2D::RIGHT))  { n += 1; }
                    if (cell.complex_face(Side2D::BOTTOM)) { n += 1; }
                    if (cell.complex_face(Side2D::TOP))    { n += 1; }
                    return n;
                }
            } else {
                return 8;
                // return hex_only ? 8 : 27; // Для VTK_TRIQUADRIC_HEXAHEDRON
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент
            return cell.node_count();
        }
    }

    // Число граней ячейки + число вершин на каждой грани
    index_t n_fverts(const EuCell& cell) const {
        if (!polyhedral) { return 0; }

        index_t res = 0;
        for (auto& face: cell.faces()) {
            if (face.is_undefined()) continue;

            // Допускаются грани с числом вершин до 8
            res += face.n_vertices() + 1;
        }
        return res;
    }

    // Записать в файл координаты вершин элемента
    void write_points(std::ofstream &file, const EuCell& cell) const {
        if (cell.adaptive()) {
            // Адаптивная ячейка
            if (cell.dim() < 3) {
                if (hex_only) {
                    const SqQuad& vertices = cell.mapping<2>();
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    file.write((byte_ptr) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    return;
                } else {
                    // Сохраняем как полигон
                    const SqQuad& vertices = cell.mapping<2>();

                    file.write((byte_ptr) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    if (cell.complex_face(Side2D::BOTTOM)) {
                        file.write((byte_ptr) vertices.vs<0, -1>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    if (cell.complex_face(Side2D::RIGHT)) {
                        file.write((byte_ptr) vertices.vs<+1, 0>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    if (cell.complex_face(Side2D::TOP)) {
                        file.write((byte_ptr) vertices.vs<0, +1>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    if (cell.complex_face(Side2D::LEFT)) {
                        file.write((byte_ptr) vertices.vs<-1, 0>().data(), 3 * sizeof(double));
                    }
                    return;
                }
            } else {
                const SqCube& vertices = cell.mapping<3>();

                // Сохраняем как простой шестигранник (VTK_HEXAHEDRON)
                file.write((byte_ptr) vertices.vs<-1, -1, -1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<+1, -1, -1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<+1, +1, -1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<-1, +1, -1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<-1, -1, +1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<+1, -1, +1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<+1, +1, +1>().data(), 3 * sizeof(double));
                file.write((byte_ptr) vertices.vs<-1, +1, +1>().data(), 3 * sizeof(double));
                return;
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент

            // В данном случае подразумеваем, что вершины пронумерованы
            // в соответствии с соглашениями, принятыми в формате VTK

            index_t n_vertices = cell.node_count();
            file.write((byte_ptr) cell.vertices_data(), 3 * n_vertices * sizeof(double));
        }
    }

    // Записать порядок вершин элемента
    void write_connectivity(std::ofstream &file, const EuCell& cell, index_t &counter) const {
        index_t n = n_points(cell);
        for (index_t i = 0; i < n; ++i) {
            index_t val = counter++;
            file.write((byte_ptr) &val, sizeof(index_t));
        }
    }
};

void write_mesh_header(
        std::ofstream &file, AmrCells &cells,
        const Variables &variables, bool hex_only, bool polyhedral
) {
    const Handler handler(hex_only, polyhedral);

    const index_t n_cells = cells.n_cells();

    // Количество вершин
    index_t n_points = 0;
    index_t n_fverts = 0;
    for (auto& cell: cells) {
        n_points += handler.n_points(cell);
        n_fverts += handler.n_fverts(cell);
    }

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    offset_t byte_offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    file << "      </Points>\n";
    byte_offset += sizeof(datasize_t) + 3 * n_points * sizeof(double);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_points * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<type_t>() << "\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(type_t);

    if (polyhedral) {
        file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"faces\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + (n_cells + n_fverts) * sizeof(index_t);

        file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"faceoffsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + n_cells * sizeof(index_t);
    }

    file << "      </Cells>\n";

    // CellData
    file << "      <CellData>\n";

    for (auto &field: variables.list()) {
        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";

        byte_offset += sizeof(datasize_t) + n_cells * field.size();
    }

    file << "      </CellData>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_mesh_primitives(
        std::ofstream &file, AmrCells &cells,
        bool hex_only, bool polyhedral
) {
    const Handler handler(hex_only, polyhedral);

    const index_t n_cells = cells.n_cells();

    // Количество вершин
    index_t n_points = 0;
    index_t n_fverts = 0;
    for (auto& cell: cells) {
        n_points += handler.n_points(cell);
        n_fverts += handler.n_fverts(cell);
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    datasize_t data_size = 3 * n_points * sizeof(double);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto &cell: cells) {
        handler.write_points(file, cell);
    }

    // Cells
    // Connectivity
    data_size = n_points * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t counter = 0;
    for (auto &cell: cells) {
        handler.write_connectivity(file, cell, counter);
    }

    // Offsets
    data_size = n_cells * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t p_index = 0;
    for (auto& cell: cells) {
        p_index += handler.n_points(cell);
        file.write((byte_ptr) &p_index, sizeof(index_t));
    }

    // Types
    data_size = n_cells * sizeof(type_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto& cell: cells) {
        type_t type = handler.type(cell);
        file.write((byte_ptr) &type, sizeof(type_t));
    }

    // Данные многогранников
    if (polyhedral) {
        // Faces
        data_size = (n_cells + n_fverts) * sizeof(index_t);
        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        counter = 0;
        for (auto &cell: cells) {
            index_t nf = cell.face_count();
            file.write((byte_ptr) &nf, sizeof(index_t));

            // Массив для описания грани
            std::array<index_t, AmrFaces::max_vertices + 1> some_face;
            for (auto &face: cell.faces()) {
                if (face.is_undefined()) { continue; }

                some_face[0] = face.n_vertices();
                for (int j = 0; j < some_face[0]; ++j) {
                    some_face[j + 1] = face.vertex_index(j);
                }
                file.write((byte_ptr) some_face.data(), (some_face[0] + 1) * sizeof(index_t));
            }
            counter += handler.n_points(cell);
        }

        // FaceOffsets
        data_size = n_cells * sizeof(index_t);
        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        p_index = 0;
        for (auto& cell: cells) {
            p_index += (1 + handler.n_fverts(cell));
            file.write((byte_ptr) &p_index, sizeof(index_t));
        }
    }
}

void write_cells_data(
        std::ofstream &file, AmrCells &cells,
        const Variables &variables
) {
    std::vector<char> temp;

    const index_t n_cells = cells.n_cells();

    for (auto &field: variables.list()) {
        index_t field_size = field.size();
        datasize_t data_size = n_cells * field_size;

        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        index_t counter = 0;
        for (auto& cell: cells) {
            field.write(cell, temp.data() + counter * field_size);
            ++counter;
        }

        file.write((byte_ptr) temp.data(), data_size);
    }
}

/// ===========================================================================
///                            Функции класса
/// ===========================================================================

VtuFile::VtuFile(
        const std::string &filename,
        const Variables &variables,
        bool hex_only, bool polyhedral, bool unique_nodes) :
    filename(filename), variables(variables),
    hex_only(hex_only), polyhedral(polyhedral) {
}

void VtuFile::save(EuMesh &mesh) const {
    if (unique_nodes) {
        throw std::runtime_error("VtuFile save unique nodes: Not implemented");
        //mesh.collect_nodes();
    }
    //if (mesh.has_nodes()) {
    //    save(mesh.locals(), mesh.nodes());
    //}
    //else {
        save(mesh.locals());
    //}
}

void VtuFile::save(AmrCells &cells) const {
    save(filename, cells, variables, hex_only, polyhedral, unique_nodes);
}

void VtuFile::save(
        const std::string &filename, AmrCells &locals,
        const Variables &variables, bool hex_only, bool polyhedral, bool unique_nodes
) {
    // Создать директории, если указано сложное имя filename
    fs::path file_path(filename);
    fs::path dir_path = file_path.parent_path();

    if (!dir_path.empty()) {
        fs::create_directories(dir_path);
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_mesh_header(file, locals, variables, hex_only, polyhedral);
    write_mesh_primitives(file, locals, hex_only, polyhedral);
    write_cells_data(file, locals, variables);

    file.close();
}

} // namespace zephyr::io