#include <iomanip>
#include <fstream>
#include <filesystem>

#include <zephyr/mesh/primitives/mov_node.h>
#include <zephyr/mesh/primitives/mov_cell.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/lagrange/la_mesh.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/json.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>

namespace zephyr::io {

using zephyr::utils::mpi;

inline bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

PvdFile::PvdFile()
    : variables(),
      filter(),
      hex_only(true),
      unique_nodes(false),
      clear_dir(true),
      m_open(false),
      m_counter(0)
{
}

const std::string default_dir = "output";

PvdFile::PvdFile(const utils::Json& config)
    : PvdFile() {

    std::string filename;
    std::string directory;

    if (config["directory"]) {
        directory = config["directory"].as<std::string>();
    } else {
        directory = default_dir;
    }

    if (!config["filename"]) {
        throw std::runtime_error("PvdFile(json) error: Add key 'filename'");
    }
    filename = config["filename"].as<std::string>();

    if (config["hex_only"]) {
        hex_only = config["hex_only"].as<bool>();
    }
    if (config["unique_nodes"]) {
        unique_nodes = config["unique_nodes"].as<bool>();
    }
    if (config["clear_dir"]) {
        clear_dir = config["clear_dir"].as<bool>();
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

    // Такой директории не существует, создать
    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        // Ждем пока все процессы убедятся, что директории
        // не существует, после мастер создает.
        mpi::barrier();
        if (mpi::master()) {
            fs::create_directory(directory);
        }
        mpi::barrier();
    }
    else {
        // Директория существует
        if (clear_dir) {
            mpi::barrier();
            if (mpi::master()) {
                fs::remove_all(directory);
                fs::create_directory(directory);
            }
            mpi::barrier();
        } else {
            bool success = false;
            for (int i = 1; i < 9; ++i) {
                std::string appendix = "_" + std::to_string(i);
                fs::path test_dir = directory.string() + appendix;
                if (!fs::exists(test_dir) || !fs::is_directory(test_dir)) {
                    // Ждем пока все процессы убедятся, что директории
                    // не существует, после мастер создает.
                    mpi::barrier();
                    if (mpi::master()) {
                        directory += appendix;
                        fs::create_directory(directory);
                    }
                    mpi::barrier();
                    success = true;
                    break;
                }
            }
            if (!success) {
                std::string message = "PvdFile open error: only 10 directories with the same name '"
                                      + directory.filename().string() + "'";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
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

    std::string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
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
        save(mesh.locals(), mesh.nodes(), timestep);
    }
    else {
        save(mesh.locals(), timestep);
    }
}

void PvdFile::save(AmrStorage& elements, double timestep) {
    VtuFile::save(get_filename(), elements, variables, hex_only, polyhedral, filter);
    update_pvd(timestep);
}

void PvdFile::save(AmrStorage& elements, const std::vector<geom::Vector3d>& nodes, double timestep) {
    VtuFile::save(get_filename(), elements, nodes, variables, hex_only, polyhedral);
    update_pvd(timestep);
}

void PvdFile::save(zephyr::mesh::LaMesh& mesh, double timestep) {
    save(mesh.locals(), mesh.nodes(), timestep);
}

void PvdFile::save(CellStorage& cells, NodeStorage& nodes, double timestep) {
    VtuFile::save(get_filename(), cells, nodes, variables, filter);
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
