#include <fstream>

#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/mov_node.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/mov_cell.h>

#include <zephyr/io/vtu_file.h>

namespace zephyr { namespace io {

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

inline size_t count_cells(AmrStorage &cells, const Filter &filter) {
    if (filter.is_trivial()) {
        return cells.size();
    }
    size_t count = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            ++count;
        }
    }
    return count;
}

struct Handler {
    bool hex_only;

    Handler(bool hex_only = false)
            : hex_only(hex_only) {
    }

    /// @brief VTK тип ячейки
    int type(const geom::AmrCell& cell) {
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
            int n = cell.vertices.count();

            if (cell.dim < 3) {
                // Двумерный полигон
                switch (n) {
                    case 3:
                        return 5; // VTK_TRIANGLE
                    case 4:
                        return 9; // VTK_QUAD
                    default:
                        return 7; // VTK_POLYGON
                }
            } else {
                // Один из доступных примитивов
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
    static int type(const geom::MovCell& cell) {
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
    int n_points(geom::AmrCell& cell) {
        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    return 4;
                } else {
                    auto& faces = cell.faces;

                    int n = 4;
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
    static int n_points(geom::MovCell& cell) {
        return cell.nodes.count();
    }

    /// @brief Записать в файл координаты вершин элемента
    void write_points(std::ofstream &file, const geom::AmrCell& cell) {
        auto& vertices = cell.vertices;

        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    file.write((char *) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    file.write((char *) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    file.write((char *) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    file.write((char *) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    return;
                } else {
                    // Сохраняем как полигон
                    auto& faces = cell.faces;

                    file.write((char *) vertices.vs<-1, -1>().data(), 3 * sizeof(double));
                    if (faces[Side::BOTTOM1].is_actual()) {
                        file.write((char *) vertices.vs<0, -1>().data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices.vs<+1, -1>().data(), 3 * sizeof(double));
                    if (faces[Side::RIGHT1].is_actual()) {
                        file.write((char *) vertices.vs<+1, 0>().data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices.vs<+1, +1>().data(), 3 * sizeof(double));
                    if (faces[Side::TOP1].is_actual()) {
                        file.write((char *) vertices.vs<0, +1>().data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices.vs<-1, +1>().data(), 3 * sizeof(double));
                    if (faces[Side::LEFT1].is_actual()) {
                        file.write((char *) vertices.vs<-1, 0>().data(), 3 * sizeof(double));
                    }
                    return;
                }
            } else {
                // Сохраняем как простой шестигранник (VTK_HEXAHEDRON)
                file.write((char *) vertices.vs<-1, -1, -1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<+1, -1, -1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<+1, +1, -1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<-1, +1, -1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<-1, -1, +1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<+1, -1, +1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<+1, +1, +1>().data(), 3 * sizeof(double));
                file.write((char *) vertices.vs<-1, +1, +1>().data(), 3 * sizeof(double));
                return;
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент

            // В данном случае подразумеваем, что вершины пронумерованы
            // в соответствии с соглашениями, принятыми в формате VTK

            auto n_vertices = vertices.count();

            file.write((char *) vertices[0].data(), 3 * n_vertices * sizeof(double));
        }
    }

    /// @brief Записать порядок вершин элемента
    void write_connectivity(std::ofstream &file, geom::AmrCell& cell, size_t &counter) {
        auto n = n_points(cell);
        for (int i = 0; i < n; ++i) {
            size_t val = counter++;
            file.write((char *) &val, sizeof(uint64_t));
        }
    }
    
    /// @brief Вариант вызывается при заданых nodes.
    void write_connectivity2(std::ofstream &file, geom::AmrCell& cell, size_t &counter) {
        auto& nodes = cell.nodes;

        if (cell.adaptive) {
            // Адаптивная ячейка
            if (cell.dim < 3) {
                if (hex_only) {
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    const size_t n = 4;
                    std::array<size_t, n> vals = {
                            uint64_t(nodes[SqQuad::iss<-1, -1>()]),
                            uint64_t(nodes[SqQuad::iss<+1, -1>()]),
                            uint64_t(nodes[SqQuad::iss<+1, +1>()]),
                            uint64_t(nodes[SqQuad::iss<-1, +1>()])
                    };

                    file.write((char *) vals.data(), n * sizeof(uint64_t));
                    counter += n;
                    return;
                } else {
                    // Сохраняем как полигон
                    auto &faces = cell.faces;

                    size_t n = 0;
                    std::array<uint64_t, 8> vals;

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

                    file.write((char *) vals.data(), n * sizeof(uint64_t));
                    counter += n;
                    return;
                }
            } else {
                // Сохраняем как простой шестигранник (VTK_HEXAHEDRON)
                const size_t n = 8;
                std::array<uint64_t, n> vals = {
                        uint64_t(nodes[SqCube::iss<-1, -1, -1>()]),
                        uint64_t(nodes[SqCube::iss<+1, -1, -1>()]),
                        uint64_t(nodes[SqCube::iss<+1, +1, -1>()]),
                        uint64_t(nodes[SqCube::iss<-1, +1, -1>()]),
                        uint64_t(nodes[SqCube::iss<-1, -1, +1>()]),
                        uint64_t(nodes[SqCube::iss<+1, -1, +1>()]),
                        uint64_t(nodes[SqCube::iss<+1, +1, +1>()]),
                        uint64_t(nodes[SqCube::iss<-1, +1, +1>()]),

                };
                file.write((char *) vals.data(), n * sizeof(uint64_t));
                counter += n;
                return;
            }
        }
        else {
            // В обратном случае это полигон или обычный трехмерный элемент

            // В данном случае подразумеваем, что вершины пронумерованы
            // в соответствии с соглашениями, принятыми в формате VTK
            
            auto n = n_points(cell);
            for (int i = 0; i < n; ++i) {
                uint64_t val = nodes[i];
                file.write((char *) &val, sizeof(uint64_t));
            }
            counter += n;
        }
    }

    /// @brief Записать порядок вершин элемента
    static void write_connectivity(std::ofstream &file, geom::MovCell& cell, size_t &counter) {
        auto n = cell.nodes.count();
        for (int i = 0; i < n; ++i) {
            size_t val = cell.nodes[i];
            file.write((char *) &val, sizeof(uint64_t));
            ++counter;
        }
    }
};

void write_mesh_header(
        std::ofstream &file, AmrStorage &cells, size_t n_cells,
        const Variables &variables, bool hex_only, const Filter &filter
) {
    Handler handler(hex_only);

    // Количество вершин
    size_t n_points = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell);
        }
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    size_t offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    file << "      </Points>\n";
    offset += 3 * n_points * sizeof(double) + sizeof(uint32_t);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
    offset += n_points * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += n_cells * sizeof(uint8_t) + sizeof(uint32_t);

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
        file << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

        offset += n_cells * field.size() + sizeof(uint32_t);
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
    int n_nodes = nodes.size();
    int n_cells = cells.size();
    int n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += handler.n_points(cell);
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    size_t offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    file << "      </Points>\n";
    offset += 3 * n_nodes * sizeof(double) + sizeof(uint32_t);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
    offset += n_connectivity * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += n_cells * sizeof(uint8_t) + sizeof(uint32_t);

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
        file << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

        offset += n_cells * field.size() + sizeof(uint32_t);
    }
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_mesh_header(
        std::ofstream &file, CellStorage &cells,
        NodeStorage& nodes, const Variables &variables
) {
    int n_nodes = nodes.size();
    int n_cells = cells.size();
    int n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += Handler::n_points(cell);
    }

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    size_t offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    file << "      </Points>\n";
    offset += 3 * n_nodes * sizeof(double) + sizeof(uint32_t);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
    offset += n_connectivity * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += n_cells * sizeof(uint8_t) + sizeof(uint32_t);

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
        file << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

        offset += n_cells * field.size() + sizeof(uint32_t);
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
        file << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

        offset += n_nodes * field.size() + sizeof(uint32_t);
    }
    file << "      </PointData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_mesh_primitives(
        std::ofstream &file, AmrStorage &cells, size_t n_cells,
        const Variables &variables, bool hex_only, const Filter &filter
) {
    Handler handler(hex_only);

    // Количество вершин
    size_t n_points = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell);
        }
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    uint32_t data_size = static_cast<uint32_t>(3 * n_points * sizeof(double));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto &cell: cells) {
        if (filter(cell)) {
            handler.write_points(file, cell);
        }
    }

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(n_points * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    size_t counter = 0;
    for (auto &cell: cells) {
        if (filter(cell)) {
            handler.write_connectivity(file, cell, counter);
        }
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    uint64_t offset = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            offset += handler.n_points(cell);
            file.write((char *) &(offset), sizeof(uint64_t));
        }
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto& cell: cells) {
        if (filter(cell)) {
            uint8_t type = handler.type(cell);
            file.write((char *) &type, sizeof(uint8_t));
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
    int n_nodes = nodes.size();
    int n_cells = cells.size();
    int n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += handler.n_points(cell);
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    uint32_t data_size = static_cast<uint32_t>(3 * n_nodes * sizeof(double));
    file.write((char *) &data_size, sizeof(uint32_t));
    file.write((char *) nodes.data(), data_size);

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(n_connectivity * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    size_t counter = 0;
    for (auto &cell: cells) {
        handler.write_connectivity2(file, cell, counter);
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    uint64_t offset = 0;
    for (auto& cell: cells) {
        offset += handler.n_points(cell);
        file.write((char *) &(offset), sizeof(uint64_t));
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto& cell: cells) {
        uint8_t type = handler.type(cell);
        file.write((char *) &type, sizeof(uint8_t));
    }
}

void write_mesh_primitives(
        std::ofstream &file, CellStorage &cells,
        NodeStorage &nodes, const Variables &variables
) {
    // Количество вершин
    int n_nodes = nodes.size();
    int n_cells = cells.size();
    int n_connectivity = 0;
    for (auto& cell: cells) {
        n_connectivity += Handler::n_points(cell);
    }

    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    uint32_t data_size = static_cast<uint32_t>(3 * n_nodes * sizeof(double));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto &node: nodes) {
        double *data = node.coords.data();
        file.write((char *) data, 3 * sizeof(double));
    }

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(n_connectivity * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    size_t counter = 0;
    for (auto &cell: cells) {
        Handler::write_connectivity(file, cell, counter);
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    uint64_t offset = 0;
    for (auto& cell: cells) {
        offset += Handler::n_points(cell);
        file.write((char *) &(offset), sizeof(uint64_t));
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto& cell: cells) {
        uint8_t type = Handler::type(cell);
        file.write((char *) &type, sizeof(uint8_t));
    }
}

void write_cells_data(
        std::ofstream &file, AmrStorage &cells, size_t n_cells,
        const Variables &variables, const Filter &filter
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        if (!field.is_eu_cell()) {
            continue;
        }

        size_t field_size = field.size();
        uint32_t data_size = n_cells * field_size;

        file.write((char *) &data_size, sizeof(uint32_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto& cell: cells) {
            if (filter(cell)) {
                field.write(cell, temp.data() + counter * field_size);
                ++counter;
            }
        }

        file.write((char *) temp.data(), data_size);
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

        size_t field_size = field.size();
        uint32_t data_size = cells.size() * field_size;

        file.write((char *) &data_size, sizeof(uint32_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto& cell: cells) {
            field.write(cell, temp.data() + counter * field_size);
            ++counter;
        }

        file.write((char *) temp.data(), data_size);
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

        size_t field_size = field.size();
        uint32_t data_size = cells.size() * field_size;

        file.write((char *) &data_size, sizeof(uint32_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto& cell: cells) {
            field.write(cell, temp.data() + counter * field_size);
            ++counter;
        }

        file.write((char *) temp.data(), data_size);
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

        size_t field_size = field.size();
        uint32_t data_size = nodes.size() * field_size;

        file.write((char *) &data_size, sizeof(uint32_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto& node: nodes) {
            field.write(node, temp.data() + counter * field_size);
            ++counter;
        }

        file.write((char *) temp.data(), data_size);
    }
}

/// ===========================================================================
///                            Функции класса
/// ===========================================================================

VtuFile::VtuFile(
        const std::string &filename,
        const Variables &variables, bool hex_only) :
    filename(filename), variables(variables),
    filter(), hex_only(hex_only) {
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
    save(filename, cells, variables, hex_only, filter);
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
    bool hex_only, const Filter &filter
) {
    size_t n_cells = count_cells(cells, filter);
    if (n_cells < 1) {
        return;
    }

    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << filename << "'\n";
        return;
    }

    write_mesh_header(file, cells, n_cells, variables, hex_only, filter);
    write_mesh_primitives(file, cells, n_cells, variables, hex_only, filter);
    write_cells_data(file, cells, n_cells, variables, filter);

    file.close();
}


void VtuFile::save(
        const std::string &filename, AmrStorage &cells,
        const std::vector<geom::Vector3d>& nodes,
        const Variables &variables, bool hex_only) {

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

} // io
} // zephyr