#include <iostream>
#include <sstream>
#include <fstream>

#include <zephyr/io/vtr_file.h>

namespace zephyr::io {

VtrFile::VtrFile(const std::string &filename,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z)
    : m_filename(filename), m_x(x), m_y(y), m_z(z) {

    if (m_x.empty() || m_y.empty()) {
        throw std::runtime_error("Empty nodes coordinates (x or y)");
    }
    m_dimension = m_z.empty() ? 2 : 3;
    if (m_z.empty()) {
        m_z = {0.0};
    }
}

void VtrFile::save() const {
    create_directories(m_filename);

    std::ofstream file(m_filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file '" << m_filename << "'\n";
        return;
    }

    std::stringstream extent;
    extent << "0 " << m_x.size() - 1 << " 0 " << m_y.size() - 1 << " 0 " << m_z.size() - 1;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"" << byteorder() << "\">\n";
    file << "  <RectilinearGrid WholeExtent=\"" << extent.str() << "\">\n";
    file << "    <Piece Extent=\"" << extent.str() << "\">\n";

    offset_t byte_offset = 0;

    file << "      <Coordinates>\n";
    file << "        <DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(m_x);
    file << "        <DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(m_y);
    file << "        <DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(datasize_t) + buffer_size(m_z);
    file << "      </Coordinates>\n";

    file << "      <PointData>\n";
    for (const auto&[name, field]: m_point_data) {
        file << "        <DataArray type=\"" << field.type << "\" Name=\"" << name;
        file << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
        byte_offset += sizeof(datasize_t) + buffer_size(field.data);
    }
    file << "      </PointData>\n";

    file << "      <CellData>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </RectilinearGrid>\n";

    file << "  <AppendedData encoding=\"raw\">\n";
    file << "_";

    write_buffer(file, m_x);
    write_buffer(file, m_y);
    write_buffer(file, m_z);

    for (const auto&[name, field]: m_point_data) {
        write_buffer(file, field.data);
    }

    file << "\n  </AppendedData>\n";
    file << "</VTKFile>\n";

    file.close();
}

} // namespace zephyr::io