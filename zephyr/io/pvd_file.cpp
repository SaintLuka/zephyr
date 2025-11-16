#include <iomanip>
#include <fstream>
#include <filesystem>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/json.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::io {

using utils::mpi;

PvdFile::PvdFile()
    : m_open(false), m_counter(0)
{
}

const std::string default_dir = "output";

PvdFile::PvdFile(const utils::Json& config) : PvdFile() {
    std::string directory;

    if (config["directory"]) {
        directory = config["directory"].as<std::string>();
    } else {
        directory = default_dir;
    }

    if (!config["filename"]) {
        throw std::runtime_error("PvdFile(json) error: Add key 'filename'");
    }
    std::string filename = config["filename"].as<std::string>();

    if (config["polyhedral"]) {
        polyhedral = config["polyhedral"].as<bool>();
    }
    if (config["unique_nodes"]) {
        unique_nodes = config["unique_nodes"].as<bool>();
    }

    open(filename, directory);
}

void PvdFile::open(const char* filename) {
    open(std::string(filename), default_dir, !mpi::single());
}

void PvdFile::open(const char* filename, bool distributed) {
    open(std::string(filename), default_dir, mpi::single() ? false : distributed);
}

void PvdFile::open(const char* filename, const char* directory) {
    open(std::string(filename), std::string(directory), !mpi::single());
}

void PvdFile::open(const char* filename, const char* directory, bool distributed) {
    open(std::string(filename), std::string(directory), distributed);
}

void PvdFile::open(const std::string& filename) {
    open(filename, default_dir, !mpi::single());
}

void PvdFile::open(const std::string& filename, bool distributed) {
    open(filename, default_dir, mpi::single() ? false : distributed);
}

void PvdFile::open(const std::string& filename, const std::string& directory) {
    open(filename, directory, !mpi::single());
}

void PvdFile::open(const std::string& filename, const std::string& _directory, bool distributed) {
    namespace fs = std::filesystem;

    if (m_open) {
        return;
    }

    fs::path directory = fs::current_path();
    if (!_directory.empty()) {
        fs::path dir = _directory;
        if (dir.is_relative()) {
            directory /= dir;
        } else {
            directory = dir;
        }
    }

    // Мастер проверяет наличие директории и создает при необходимости
    if (mpi::master()) {
        if (!fs::exists(directory) || !fs::is_directory(directory)) {
            fs::create_directories(directory);
        }
    }

    m_filename = filename;
    if (filename.size() > 4) {
        if (filename.substr(filename.size() - 4) == ".pvd") {
            m_filename = filename.substr(filename.size() - 4);
        }
    }
    m_fullname = (directory / filename).string();

    // Мастер-процесс пишет заголовок PVD
    m_distributed = mpi::single() ? false : distributed;
    if (m_distributed && !mpi::master()) {
        return;
    }

    /// Откроем файл и запишем заголовок
    std::ofstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Warning: Cannot open .pvd file " << m_fullname << ".pvd\n";
        return;
    }

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    ofs << "    <Collection>" << std::endl;

    m_pos = ofs.tellp();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    m_open = true;
}

void PvdFile::save(mesh::EuMesh& mesh, double timestep) {
    if (unique_nodes) {
        mesh.collect_nodes();
    }
    if (mesh.has_nodes()) {
        VtuFile::save(get_filename(), mesh.locals(), mesh.nodes(), variables, polyhedral);
    }
    else {
        VtuFile::save(get_filename(), mesh.locals(), variables, polyhedral);
    }
    update_pvd(timestep);
}

void PvdFile::save(mesh::AmrCells& elements, double timestep) {
    VtuFile::save(get_filename(), elements, variables, polyhedral);
    update_pvd(timestep);
}

std::string PvdFile::get_filename() const {
    std::string filename = m_fullname + "_" + std::to_string(m_counter);

    if (m_distributed) {
        filename += ".pt" + mpi::srank();
    }

    filename += ".vtu";
    return filename;
}

void PvdFile::update_pvd(double timestep) {
    // Мастер-процесс пишет PVD
    if (m_distributed && !mpi::master()) {
        ++m_counter;
        return;
    }

    if (!m_open) {
        std::string message = "PvdFile::save() error: You need to open PvdFile";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    std::fstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << m_fullname << ".pvd\n";
    }

    ofs.seekg(m_pos, std::ios::beg);

    ofs << std::scientific << std::setprecision(15);

    if (m_distributed) {
        for (int r = 0; r < mpi::size(); ++r) {
            ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"" << r << "\" file=\""
                << m_filename << "_" << m_counter << ".pt" << r << ".vtu" << "\"/>\n";
        }
    }
    else {
        ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << m_filename << "_" << m_counter << ".vtu" << "\"/>\n";
    }

    m_pos = ofs.tellg();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    ++m_counter;
}

} // namespace zephyr::io