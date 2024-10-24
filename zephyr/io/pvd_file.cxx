#include <iomanip>
#include <fstream>
#include <filesystem>

#include <zephyr/geom/primitives/mov_node.h>
#include <zephyr/geom/primitives/mov_cell.h>
#include <zephyr/geom/primitives/amr_cell.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/mesh/lagrange/la_mesh.h>

namespace zephyr::io {

namespace fs = std::filesystem;

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
Network& dummy_net = ::zephyr::network::network::dummy;
#endif

#ifdef ZEPHYR_ENABLE_MPI
#include "zephyr/utils/mpi.h"
using namespace zephyr::utils;
#endif

inline bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

PvdFile::PvdFile()
    : variables(), filter(), hex_only(true), m_counter(0)
    , m_open(false)
#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    , net(dummy_net)
#endif
{
// #ifdef ZEPHYR_ENABLE_DISTRIBUTED
    // if (network::network::main().size() > 1) {
        // std::cerr << "Warning: MPI is enabled, but single process version of PvdFile constructor is used.\n";
    // }
// #endif
}

PvdFile::PvdFile(const std::string &filename, const std::string &directory)
    : variables(), filter(), hex_only(true), m_counter(0)
    , m_open(false)
#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    , net(dummy_net)
#endif
{
// #ifdef ZEPHYR_ENABLE_DISTRIBUTED
    // if (network::network::main().size() > 1) {
     // std::cerr << "Warning: MPI is enabled, but single process version of PvdFile constructor is used.\n";
    // }
// #endif
    open(filename, directory);
}

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
PvdFile::PvdFile(Network& network)
    : net(network), variables(), filter(), hex_only(false), m_counter(0)
    , m_open(false) {
}

PvdFile::PvdFile(Network& network, const std::string &filename, const std::string &directory)
    : net(network), variables(), filter(), hex_only(false), m_counter(0)
    , m_open(false) {
    open(filename, directory);
}
#endif

void PvdFile::open(const std::string& filename, const std::string& _directory) {
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

    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        fs::create_directories(directory);
    }

    m_filename = filename;
    if (filename.size() > 4) {
        if (filename.substr(filename.size() - 4) == ".pvd") {
            m_filename = filename.substr(filename.size() - 4);
        }
    }
    m_fullname = (directory / filename).string();

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    // Не мастер процесс
    if (&net != &dummy_net && !net.is_master()) {
        m_open = true;
        return;
    }
#endif

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
    if (!m_open) {
        throw std::runtime_error("PvdFile::save() error: You need to open PvdFile");
    }
    VtuFile::save(get_filename(), elements, variables, hex_only, filter);
    update_pvd(timestep);
}

void PvdFile::save(AmrStorage& elements, const std::vector<geom::Vector3d>& nodes, double timestep) {
    if (!m_open) {
        throw std::runtime_error("PvdFile::save() error: You need to open PvdFile");
    }
    VtuFile::save(get_filename(), elements, nodes, variables, hex_only);
    update_pvd(timestep);
}

void PvdFile::save(zephyr::mesh::LaMesh& mesh, double timestep) {
    save(mesh.locals(), mesh.nodes(), timestep);
}

void PvdFile::save(CellStorage& cells, NodeStorage& nodes, double timestep) {
    if (!m_open) {
        throw std::runtime_error("PvdFile::save() error: You need to open PvdFile");
    }
    VtuFile::save(get_filename(), cells, nodes, variables, filter);
#ifdef ZEPHYR_ENABLE_MPI
    if(mpi::master())
#endif
    update_pvd(timestep);
}

std::string PvdFile::get_filename() const {
    std::string filename = m_fullname;

#ifdef ZEPHYR_ENABLE_MPI
    filename += "_" + std::to_string(mpi::rank());
#endif

    filename += "_" + std::to_string(m_counter);

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    if (&net != &dummy_net) {
        filename += ".pt" + std::to_string(net.rank());
    }
#endif
    filename += ".vtu";
    return filename;
}

void PvdFile::update_pvd(double timestep) {
    if (!m_open) {
        throw std::runtime_error("You need to open PvdFile");
    }

    std::fstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << m_fullname << ".pvd\n";
    }

    ofs.seekg(m_pos, std::ios::beg);

    ofs << std::scientific << std::setprecision(15);

#ifdef ZEPHYR_ENABLE_MPI
    for(int r = 0; r < mpi::size(); ++r){
        ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"" << r << "\" file=\""
                << m_filename << "_" << r << "_" << m_counter << ".vtu" << "\"/>\n";
    }
#elif
    ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << m_filename << "_" << m_counter << ".vtu" << "\"/>\n";
#endif

    m_pos = ofs.tellg();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    ++m_counter;
}

} // namespace zephyr::io
