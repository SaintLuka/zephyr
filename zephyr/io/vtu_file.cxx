#include <fstream>

#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/mov_node.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/mov_cell.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr::io {

using namespace zephyr::geom;

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

namespace {
// Тип для задания размерва массивов в бинарном файле,
// я не уверен, что его можно менять
using datasize_t = std::uint32_t;
using offset_t   = std::uint32_t;  // Тип смещения массивов в байтах
using index_t    = std::uint32_t;  // Тип нумерации примитивов (ячеек, вершин и т.д.)
using type_t     = std::uint8_t;   // VTK тип примитива/ячейки)

using byte_ptr = char*;
}

inline index_t count_cells(AmrStorage &cells, const Filter &filter) {
    if (filter.is_trivial()) {
        return cells.size();
    }
    index_t count = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            ++count;
        }
    }
    return count;
}

struct Handler {
    bool hex_only, polyhedral;

    Handler(bool hex_only = false, bool polyhedral = false)
            : hex_only(hex_only), polyhedral(polyhedral) {
    }

    /// @brief VTK тип ячейки
    type_t type(const geom::AmrCell& cell) {
        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    return 9; // VTK_QUAD
                } else {
                    auto& faces = cell.faces;

                    if (faces[Side::LEFT1].is_actual() ||
                        faces[Side::RIGHT1].is_actual() ||
                        faces[Side::BOTTOM1].is_actual() ||
                        faces[Side::TOP1].is_actual()) {
                        return 7; // VTK_POLYGON
                    } else {
                        return 9; // VTK_QUAD
                    }
                }
            } else {
                return 12; // VTK_HEXAHEDRON

                // В будущем использовать TRIQUADRIC HEXAHEDRON
                // return 29; // VTK_TRIQUADRIC_HEXAHEDRON
            }
        }
        else {
            // Обычная эйлерова ячейка

            if (cell.dim < 3) {
                // Двумерный полигон
                index_t n = cell.vertices.count();

                switch (n) {
                    case 3:
                        return 5; // VTK_TRIANGLE
                    case 4:
                        return 9; // VTK_QUAD
                    default:
                        return 7; // VTK_POLYGON
                }
            } else {
                // Многогранник общего вида
                if (polyhedral) {
                    return 42;     // VTK_POLYHEDRON
                }

                // Один из доступных примитивов
                index_t n = cell.vertices.count();
                switch (n) {
                    case 4:
                        return 10; // VTK_TETRA
                    case 5:
                        return 14; // VTK_PYRAMID
                    case 6:
                        return 13; // VTK_WEDGE
                    default:
                        return 12; // VTK_HEXAHEDRON
                }
            }
        }
    }

    /// @brief VTK тип ячейки
    static type_t type(const geom::MovCell& cell) {
        if (cell.dim < 3) {
            // Двумерный полигон
            switch (cell.nodes.count()) {
                case 3:
                    return 5; // VTK_TRIANGLE
                case 4:
                    return 9; // VTK_QUAD
                default:
                    return 7; // VTK_POLYGON
            }
        } else {
            // Один из доступных примитивов
            switch (cell.nodes.count()) {
                case 4:
                    return 10; // VTK_TETRA
                case 5:
                    return 14; // VTK_PYRAMID
                case 6:
                    return 13; // VTK_WEDGE
                default:
                    return 12; // VTK_HEXAHEDRON
            }
        }
    }

    /// @brief Количество вершин элемента
    index_t n_points(geom::AmrCell& cell) {
        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    return 4;
                } else {
                    auto& faces = cell.faces;

                    index_t n = 4;
                    if (faces[Side::LEFT1].is_actual()) {
                        n += 1;
                    }
                    if (faces[Side::RIGHT1].is_actual()) {
                        n += 1;
                    }
                    if (faces[Side::BOTTOM1].is_actual()) {
                        n += 1;
                    }
                    if (faces[Side::TOP1].is_actual()) {
                        n += 1;
                    }
                    return n;
                }
            } else {
                return 8;
                // return hex_only ? 8 : 27; // Для VTK_TRIQUADRIC_HEXAHEDRON
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент
            return cell.vertices.count();
        }
    }

    /// @brief Количество вершин элемента
    static index_t n_points(geom::MovCell& cell) {
        return cell.nodes.count();
    }

    /// @brief Число граней ячейки + число вершин на каждой грани
    index_t n_fverts(geom::AmrCell& cell) {
        if (!polyhedral) { return 0; }

        index_t res = 0;
        for (auto& face: cell.faces) {
            if (face.is_undefined()) continue;

            // Допускаются грани с числом вершин до 9
            res += face.poly_size() + 1;
        }
        return res;
    }

    /// @brief Записать в файл координаты вершин элемента
    void write_points(std::ofstream &file, const geom::AmrCell& cell) {
        auto& vertices = cell.vertices;

        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    file.write((byte_ptr) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    file.write((byte_ptr) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    return;
                } else {
                    // Сохраняем как полигон
                    auto& faces = cell.faces;

                    file.write((byte_ptr) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    if (faces[Side::BOTTOM1].is_actual()) {
                        file.write((byte_ptr) vertices.vs<0, -1>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    if (faces[Side::RIGHT1].is_actual()) {
                        file.write((byte_ptr) vertices.vs<+1, 0>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    if (faces[Side::TOP1].is_actual()) {
                        file.write((byte_ptr) vertices.vs<0, +1>().data(), 3 * sizeof(double));
                    }

                    file.write((byte_ptr) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    if (faces[Side::LEFT1].is_actual()) {
                        file.write((byte_ptr) vertices.vs<-1, 0>().data(), 3 * sizeof(double));
                    }
                    return;
                }
            } else {
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

            index_t n_vertices = vertices.count();

            file.write((byte_ptr) vertices[0].data(), 3 * n_vertices * sizeof(double));
        }
    }

    /// @brief Записать порядок вершин элемента
    void write_connectivity(std::ofstream &file, geom::AmrCell& cell, index_t &counter) {
        index_t n = n_points(cell);
        for (index_t i = 0; i < n; ++i) {
            index_t val = counter++;
            file.write((byte_ptr) &val, sizeof(index_t));
        }
    }
    
    /// @brief Вариант вызывается при заданых nodes.
    void write_connectivity2(std::ofstream &file, geom::AmrCell& cell, index_t &counter) {
        auto& nodes = cell.nodes;

        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    const index_t n = 4;
                    std::array<index_t, n> vals = {
                            index_t(nodes[SqQuad::iss<-1, -1>()]),
                            index_t(nodes[SqQuad::iss<+1, -1>()]),
                            index_t(nodes[SqQuad::iss<+1, +1>()]),
                            index_t(nodes[SqQuad::iss<-1, +1>()])
                    };

                    file.write((byte_ptr) vals.data(), n * sizeof(index_t));
                    counter += n;
                    return;
                } else {
                    // Сохраняем как полигон
                    auto &faces = cell.faces;

                    index_t n = 0;
                    std::array<index_t, 8> vals;

                    vals[n++] = nodes[SqQuad::iss<-1, -1>()];
                    if (faces[Side::BOTTOM1].is_actual()) {
                        vals[n++] = nodes[SqQuad::iss<0, -1>()];
                    }

                    vals[n++] = nodes[SqQuad::iss<+1, -1>()];
                    if (faces[Side::RIGHT1].is_actual()) {
                        vals[n++] = nodes[SqQuad::iss<+1, 0>()];
                    }

                    vals[n++] = nodes[SqQuad::iss<+1, +1>()];
                    if (faces[Side::TOP1].is_actual()) {
                        vals[n++] = nodes[SqQuad::iss<0, +1>()];
                    }

                    vals[n++] = nodes[SqQuad::iss<-1, +1>()];
                    if (faces[Side::LEFT1].is_actual()) {
                        vals[n++] = nodes[SqQuad::iss<-1, 0>()];
                    }

                    file.write((byte_ptr) vals.data(), n * sizeof(index_t));
                    counter += n;
                    return;
                }
            } else {
                // Сохраняем как простой шестигранник (VTK_HEXAHEDRON)
                const index_t n = 8;
                std::array<index_t, n> vals = {
                        index_t(nodes[SqCube::iss<-1, -1, -1>()]),
                        index_t(nodes[SqCube::iss<+1, -1, -1>()]),
                        index_t(nodes[SqCube::iss<+1, +1, -1>()]),
                        index_t(nodes[SqCube::iss<-1, +1, -1>()]),
                        index_t(nodes[SqCube::iss<-1, -1, +1>()]),
                        index_t(nodes[SqCube::iss<+1, -1, +1>()]),
                        index_t(nodes[SqCube::iss<+1, +1, +1>()]),
                        index_t(nodes[SqCube::iss<-1, +1, +1>()]),

                };
                file.write((byte_ptr) vals.data(), n * sizeof(index_t));
                counter += n;
                return;
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент

            // В данном случае подразумеваем, что вершины пронумерованы
            // в соответствии с соглашениями, принятыми в формате VTK
            
            index_t n = n_points(cell);
            for (index_t i = 0; i < n; ++i) {
                index_t val = nodes[i];
                file.write((byte_ptr) &val, sizeof(index_t));
            }
            counter += n;
        }
    }

    /// @brief Записать порядок вершин элемента
    static void write_connectivity(std::ofstream &file, geom::MovCell& cell, index_t &counter) {
        index_t n = cell.nodes.count();
        for (index_t i = 0; i < n; ++i) {
            index_t val = cell.nodes[i];
            file.write((byte_ptr) &val, sizeof(index_t));
            ++counter;
        }
    }
};

void write_mesh_header(
        std::ofstream &file, AmrStorage &cells, index_t n_cells,
        const Variables &variables, bool hex_only, bool polyhedral,
        const Filter &filter
) {
    Handler handler(hex_only, polyhedral);

    // Количество вершин
    index_t n_points = 0;
    index_t n_fverts = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell);
            n_fverts += handler.n_fverts(cell);
        }
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
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
        if (!field.is_eu_cell()) {
            continue;
        }

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

void write_mesh_header(
        std::ofstream &file, AmrStorage &cells,
        const std::vector<Vector3d>& nodes,
        const Variables &variables, bool hex_only) {

    Handler handler(hex_only);

    // Количество вершин
    index_t n_nodes = nodes.size();
    index_t n_cells = cells.size();
    index_t n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += handler.n_points(cell);
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    offset_t byte_offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    file << "      </Points>\n";
    byte_offset += sizeof(datasize_t) + 3 * n_nodes * sizeof(double);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_connectivity * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<type_t>() << "\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(type_t);

    file << "      </Cells>\n";

    // CellData
    file << "      <CellData>\n";
    for (auto &field: variables.list()) {
        if (!field.is_eu_cell()) {
            continue;
        }

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

void write_mesh_header(
        std::ofstream &file, CellStorage &cells,
        NodeStorage& nodes, const Variables &variables
) {
    index_t n_nodes = nodes.size();
    index_t n_cells = cells.size();
    index_t n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += Handler::n_points(cell);
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    offset_t byte_offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    file << "      </Points>\n";
    byte_offset += sizeof(datasize_t) + 3 * n_nodes * sizeof(double);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_connectivity * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<index_t>() << "\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(index_t);

    file << "        <DataArray type=\"" << VtkType::get<type_t>() << "\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + n_cells * sizeof(type_t);

    file << "      </Cells>\n";

    // CellData
    file << "      <CellData>\n";
    for (auto &field: variables.list()) {
        if (!field.is_lag_cell()) {
            continue;
        }

        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";

        byte_offset += sizeof(datasize_t) + n_cells * field.size();
    }
    file << "      </CellData>\n";

    // PointData
    file << "      <PointData>\n";
    for (auto &field: variables.list()) {
        if (!field.is_node()) {
            continue;
        }

        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";

        byte_offset += sizeof(datasize_t) + n_nodes * field.size();
    }
    file << "      </PointData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_mesh_primitives(
        std::ofstream &file, AmrStorage &cells, index_t n_cells,
        const Variables &variables, bool hex_only, bool polyhedral,
        const Filter &filter
) {
    Handler handler(hex_only, polyhedral);

    // Количество вершин
    index_t n_points = 0;
    index_t n_fverts = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell);
            n_fverts += handler.n_fverts(cell);
        }
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    datasize_t data_size = 3 * n_points * sizeof(double);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto &cell: cells) {
        if (filter(cell)) {
            handler.write_points(file, cell);
        }
    }

    // Cells
    // Connectivity
    data_size = n_points * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t counter = 0;
    for (auto &cell: cells) {
        if (filter(cell)) {
            handler.write_connectivity(file, cell, counter);
        }
    }

    // Offsets
    data_size = n_cells * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t p_index = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            p_index += handler.n_points(cell);
            file.write((byte_ptr) &p_index, sizeof(index_t));
        }
    }

    // Types
    data_size = n_cells * sizeof(type_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto& cell: cells) {
        if (filter(cell)) {
            type_t type = handler.type(cell);
            file.write((byte_ptr) &type, sizeof(type_t));
        }
    }

    // Данные многогранников
    if (polyhedral) {
        // Faces
        data_size = (n_cells + n_fverts) * sizeof(index_t);
        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        counter = 0;
        for (auto &cell: cells) {
            if (!filter(cell)) { continue; }

            index_t nf = cell.faces.count();
            file.write((byte_ptr) &nf, sizeof(index_t));

            // Массив для описания грани
            std::array<index_t, 10> some_face;
            for (auto &face: cell.faces) {
                if (face.is_undefined()) { continue; }

                some_face[0] = face.poly_size();
                for (int j = 0; j < some_face[0]; ++j) {
                    some_face[j + 1] = face.get_poly_vertex(j);
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
            if (filter(cell)) {
                p_index += (1 + handler.n_fverts(cell));
                file.write((byte_ptr) &p_index, sizeof(index_t));
            }
        }
    }
}

void write_mesh_primitives(
        std::ofstream &file, AmrStorage &cells,
        const std::vector<Vector3d>& nodes,
        const Variables &variables, bool hex_only
) {
    Handler handler(hex_only);

    // Количество вершин
    index_t n_nodes = nodes.size();
    index_t n_cells = cells.size();
    index_t n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += handler.n_points(cell);
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    datasize_t data_size = 3 * n_nodes * sizeof(double);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));
    file.write((byte_ptr) nodes.data(), data_size);

    // Cells
    // Connectivity
    data_size = n_connectivity * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t counter = 0;
    for (auto &cell: cells) {
        handler.write_connectivity2(file, cell, counter);
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
}

void write_mesh_primitives(
        std::ofstream &file, CellStorage &cells,
        NodeStorage &nodes, const Variables &variables
) {
    // Количество вершин
    index_t n_nodes = nodes.size();
    index_t n_cells = cells.size();
    index_t n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += Handler::n_points(cell);
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    datasize_t data_size = 3 * n_nodes * sizeof(double);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto &node: nodes) {
        double *data = node.coords.data();
        file.write((byte_ptr) data, 3 * sizeof(double));
    }

    // Cells
    // Connectivity
    data_size = n_connectivity * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t counter = 0;
    for (auto &cell: cells) {
        Handler::write_connectivity(file, cell, counter);
    }

    // Offsets
    data_size = n_cells * sizeof(index_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    index_t p_index = 0;
    for (auto& cell: cells) {
        p_index += Handler::n_points(cell);
        file.write((byte_ptr) &p_index, sizeof(index_t));
    }

    // Types
    data_size = n_cells * sizeof(type_t);
    file.write((byte_ptr) &data_size, sizeof(datasize_t));

    for (auto& cell: cells) {
        type_t type = Handler::type(cell);
        file.write((byte_ptr) &type, sizeof(type_t));
    }
}

void write_cells_data(
        std::ofstream &file, AmrStorage &cells, index_t n_cells,
        const Variables &variables, const Filter &filter
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        if (!field.is_eu_cell()) {
            continue;
        }

        index_t field_size = field.size();
        datasize_t data_size = n_cells * field_size;

        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        index_t counter = 0;
        for (auto& cell: cells) {
            if (filter(cell)) {
                field.write(cell, temp.data() + counter * field_size);
                ++counter;
            }
        }

        file.write((byte_ptr) temp.data(), data_size);
    }
}

void write_cells_data(
        std::ofstream &file, AmrStorage &cells,
        const Variables &variables
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        if (!field.is_eu_cell()) {
            continue;
        }

        index_t field_size = field.size();
        datasize_t data_size = cells.size() * field_size;

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

void write_cells_data(
        std::ofstream &file, CellStorage &cells,
        const Variables &variables
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        if (!field.is_lag_cell()) {
            continue;
        }

        index_t field_size = field.size();
        datasize_t data_size = cells.size() * field_size;

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

void write_nodes_data(
        std::ofstream &file, NodeStorage &nodes,
        const Variables &variables
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        if (!field.is_node()) {
            continue;
        }

        index_t field_size = field.size();
        datasize_t data_size = nodes.size() * field_size;

        file.write((byte_ptr) &data_size, sizeof(datasize_t));

        temp.resize(data_size);

        index_t counter = 0;
        for (auto& node: nodes) {
            field.write(node, temp.data() + counter * field_size);
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
        bool hex_only, bool polyhedral) :
    filename(filename), variables(variables),
    filter(), hex_only(hex_only), polyhedral(polyhedral) {
}

void VtuFile::save(mesh::EuMesh &mesh) {
    if (unique_nodes) {
        mesh.collect_nodes();
    }
    if (mesh.has_nodes()) {
        save(mesh.locals(), mesh.nodes());
    }
    else {
        save(mesh.locals());
    }
}

void VtuFile::save(mesh::LaMesh &mesh) {
    save(mesh.locals(), mesh.nodes());
}

void VtuFile::save(AmrStorage &cells) {
    save(filename, cells, variables, hex_only, polyhedral, filter);
}

void VtuFile::save(AmrStorage &cells, const std::vector<Vector3d>& nodes) {
    save(filename, cells, nodes, variables, hex_only);
}

void VtuFile::save(CellStorage &cells, NodeStorage &nodes) {
    save(filename, cells, nodes, variables, filter);
}

void VtuFile::save(
    const std::string &filename,
    AmrStorage &cells, const Variables &variables,
    bool hex_only, bool polyhedral, const Filter &filter
) {
    index_t n_cells = count_cells(cells, filter);
    if (n_cells < 1) {
        return;
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_mesh_header(file, cells, n_cells, variables, hex_only, polyhedral, filter);
    write_mesh_primitives(file, cells, n_cells, variables, hex_only, polyhedral, filter);
    write_cells_data(file, cells, n_cells, variables, filter);

    file.close();
}


void VtuFile::save(
        const std::string &filename, AmrStorage &cells,
        const std::vector<geom::Vector3d>& nodes,
        const Variables &variables,
        bool hex_only, bool polyhedral) {

    if (cells.empty()) {
        return;
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_mesh_header(file, cells, nodes, variables, hex_only);
    write_mesh_primitives(file, cells, nodes, variables, hex_only);
    write_cells_data(file, cells, variables);

    file.close();
}

void VtuFile::save(
        const std::string &filename,
        CellStorage &cells, NodeStorage &nodes,
        const Variables &variables, const Filter &filter
) {
    if (cells.empty()) {
        return;
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_mesh_header(file, cells, nodes, variables);
    write_mesh_primitives(file, cells, nodes, variables);
    write_cells_data(file, cells, variables);
    write_nodes_data(file, nodes, variables);

    file.close();
}

} // namespace zephyr::io