#include <fstream>

#include <zephyr/geom/primitives/amr_cell.h>
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

inline size_t count_cells(Storage &cells, const Filter &filter) {
    if (filter.is_trivial()) {
        return cells.size();
    }
    size_t count = 0;
    for (auto cell: cells) {
        if (filter(cell)) {
            ++count;
        }
    }
    return count;
}

struct Handler {
    int dim;
    bool hex_only;

    Handler(int dim, bool hex_only)
            : dim(dim), hex_only(hex_only) {
    }

    /// @brief VTK тип ячейки
    int type(const geom::AmrCell& cell) {
        auto& faces = cell.faces;
        auto& vertices = cell.vertices;

        // Если последняя вершина (8 или 26) определена, значит все вершины
        // заполнены и перед нами AMR ячейка
        if (dim < 3) {
            if (vertices[8].is_actual()) {
                if (hex_only) {
                    return 9; // VTK_QUAD
                } else {
                    if (faces[Side::LEFT1].is_actual() ||
                        faces[Side::RIGHT1].is_actual() ||
                        faces[Side::BOTTOM1].is_actual() ||
                        faces[Side::TOP1].is_actual()) {
                        return 7; // VTK_POLYGON
                    } else {
                        return 9; // VTK_QUAD
                    }
                }
            }
        } else {
            if (vertices[26].is_actual()) {
                return 12; // VTK_HEXAHEDRON

                // В будущем использовать TRIQUADRIC HEXAHEDRON
                // return 29; // VTK_TRIQUADRIC_HEXAHEDRON
            }
        }

        // В обратном случае это полигон или обычный трехмерный элемент
        int n = vertices.size();

        if (dim < 3) {
            switch (n) {
                case 3:
                    return 5; // VTK_TRIANGLE
                case 4:
                    return 9; // VTK_QUAD
                default:
                    return 7; // VTK_POLYGON
            }
        } else {
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

    /// @brief Количество вершин элемента
    int n_points(geom::AmrCell& cell) {
        auto& faces = cell.faces;
        auto& vertices = cell.vertices;

        // Если последняя вершина (8 или 26) определена, значит все вершины
        // заполнены и перед нами AMR ячейка
        if (dim < 3) {
            if (vertices[8].is_actual()) {
                if (hex_only) {
                    return 4;
                } else {
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
            }
        } else {
            if (vertices[26].is_actual()) {
                return 8;
                // return hex_only ? 8 : 27; // Для VTK_TRIQUADRIC_HEXAHEDRON
            }
        }

        // В обратном случае это полигон или обычный трехмерный элемент
        return vertices.size();
    }

    /// @brief Записать в файл координаты вершин элемента
    void write_points(std::ofstream &file, const geom::AmrCell& cell) {
        auto& faces = cell.faces;
        auto& vertices = cell.vertices;

        // Если последняя вершина (8 или 26) определена, значит все вершины
        // заполнены и перед нами AMR ячейка
        if (dim < 3) {
            if (vertices[8].is_actual()) {
                if (hex_only) {
                    // Сохраняем как простой четырехугольник (VTK_QUAD)
                    file.write((char *) vertices[iww(0, 0)].data(), 3 * sizeof(double));
                    file.write((char *) vertices[iww(2, 0)].data(), 3 * sizeof(double));
                    file.write((char *) vertices[iww(2, 2)].data(), 3 * sizeof(double));
                    file.write((char *) vertices[iww(0, 2)].data(), 3 * sizeof(double));
                    return;
                } else {
                    // Сохраняем как полигон
                    file.write((char *) vertices[iww(0, 0)].data(), 3 * sizeof(double));
                    if (faces[Side::BOTTOM1].is_actual()) {
                        file.write((char *) vertices[iww(1, 0)].data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices[iww(2, 0)].data(), 3 * sizeof(double));
                    if (faces[Side::RIGHT1].is_actual()) {
                        file.write((char *) vertices[iww(2, 1)].data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices[iww(2, 2)].data(), 3 * sizeof(double));
                    if (faces[Side::TOP1].is_actual()) {
                        file.write((char *) vertices[iww(1, 2)].data(), 3 * sizeof(double));
                    }

                    file.write((char *) vertices[iww(0, 2)].data(), 3 * sizeof(double));
                    if (faces[Side::LEFT1].is_actual()) {
                        file.write((char *) vertices[iww(0, 1)].data(), 3 * sizeof(double));
                    }
                    return;
                }
            }
        } else {
            if (vertices[26].is_actual()) {
                // Сохраняем как простой шестигранник (VTK_HEXAHEDRON)
                file.write((char *) vertices[iww(0, 0, 0)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(2, 0, 0)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(2, 2, 0)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(0, 2, 0)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(0, 0, 2)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(2, 0, 2)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(2, 2, 2)].data(), 3 * sizeof(double));
                file.write((char *) vertices[iww(0, 2, 2)].data(), 3 * sizeof(double));
                return;
            }
        }

        // В обратном случае это полигон или обычный трехмерный элемент

        // В данном случае подразумеваем, что вершины пронумерованы
        // в соответствии с соглашениями, принятыми в формате VTK

        auto n_vertices = vertices.size();

        file.write((char *) vertices[0].data(), 3 * n_vertices * sizeof(double));
    }

    /// @brief Записать порядок вершин элемента
    void write_connectivity(std::ofstream &file, geom::AmrCell& cell, size_t &counter) {
        auto n = n_points(cell);
        for (int i = 0; i < n; ++i) {
            size_t val = counter++;
            file.write((char *) &val, sizeof(uint64_t));
        }
    }
};

void write_mesh_header(
        std::ofstream &file, Storage &cells, size_t n_cells,
        const Variables &variables, bool hex_only, const Filter &filter
) {
    size_t dim = 2;
    if (!cells.empty()) {
        dim = cells[0].dim();
    }

    Handler handler(dim, hex_only);

    // Количество вершин
    size_t n_points = 0;
    for (auto cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell.geom());
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

void write_particles_header(
        std::ofstream &file, Storage &cells, size_t n_cells,
        const Variables &variables, const Filter &filter
) {
    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    file << "  <UnstructuredGrid>" << '\n';
    file << "    <Piece NumberOfPoints=\"" << n_cells << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    size_t offset = 0;
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    file << "      </Points>\n";
    offset += 3 * n_cells * sizeof(double) + sizeof(uint32_t);

    // Cells
    file << "      <Cells>" << '\n';
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
         << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
         << offset << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += n_cells * sizeof(uint8_t) + sizeof(uint32_t);

    file << "      </Cells>\n";

    // PointData
    file << "      <PointData>\n";

    for (auto &field: variables.list()) {
        file << "        <DataArray type=\"" << field.type() << "\" Name=\"" << field.name();
        if (!field.is_scalar()) {
            file << "\" NumberOfComponents=\"" << field.n_components();
        }
        file << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

        offset += n_cells * field.size() + sizeof(uint32_t);
    }

    file << "      </PointData>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
}

void write_mesh_primitives(
        std::ofstream &file, Storage &cells, size_t n_cells,
        const Variables &variables, bool hex_only, const Filter &filter
) {
    size_t dim = cells[0].dim();

    Handler handler(dim, hex_only);

    // Количество вершин
    size_t n_points = 0;
    for (auto& cell: cells) {
        if (filter(cell)) {
            n_points += handler.n_points(cell.geom());
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
            handler.write_points(file, cell.geom());
        }
    }

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(n_points * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    size_t counter = 0;
    for (auto &cell: cells) {
        if (filter(cell)) {
            handler.write_connectivity(file, cell.geom(), counter);
        }
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    uint64_t offset = 0;
    for (auto cell: cells) {
        if (filter(cell)) {
            offset += handler.n_points(cell.geom());
            file.write((char *) &(offset), sizeof(uint64_t));
        }
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto cell: cells) {
        if (filter(cell)) {
            uint8_t type = handler.type(cell.geom());
            file.write((char *) &type, sizeof(uint8_t));
        }
    }
}

void write_particles_primitives(
        std::ofstream &file, Storage &cells, size_t n_cells,
        const Variables &variables, const Filter &filter
) {
    // AppendedData
    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    // PointsCoords
    uint32_t data_size = static_cast<uint32_t>(3 * n_cells * sizeof(double));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (auto &cell: cells) {
        if (filter(cell)) {
            file.write((char *) &cell.center(), 3 * sizeof(double));
        }
    }

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    size_t counter = 0;
    for (size_t i = 0; i < n_cells; ++i) {
        file.write((char *) &i, sizeof(uint64_t));
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    for (size_t i = 1; i <= n_cells; ++i) {
        file.write((char *) &i, sizeof(uint64_t));
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    file.write((char *) &data_size, sizeof(uint32_t));

    uint8_t type = 1;
    for (size_t i = 0; i < n_cells; ++i) {
        file.write((char *) &type, sizeof(uint8_t));
    }
}

void write_cells_data(
        std::ofstream &file, Storage &cells, size_t n_cells,
        const Variables &variables, const Filter &filter
) {
    std::vector<char> temp;

    for (auto &field: variables.list()) {
        size_t field_size = field.size();
        uint32_t data_size = n_cells * field_size;

        file.write((char *) &data_size, sizeof(uint32_t));

        temp.resize(data_size);

        size_t counter = 0;
        for (auto cell: cells) {
            if (filter(cell)) {
                field.write(cell, temp.data() + counter * field_size);
                ++counter;
            }
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

void VtuFile::save(Storage &cells) {
    save(filename, cells, variables, hex_only, filter);
}

void VtuFile::save(
    const std::string &filename,
    Storage &cells, const Variables &variables,
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

} // io
} // zephyr