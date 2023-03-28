#include <iomanip>
#include <fstream>

#if __cplusplus >= 201703L
#    include <filesystem>
     namespace tmpfs = ::std::filesystem;
#else
#    include <boost/filesystem/path.hpp>
#    include <boost/filesystem/operations.hpp>
     namespace tmpfs = ::boost::filesystem;
#endif

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>

namespace zephyr { namespace io {

namespace fs = ::tmpfs;

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
Network& dummy_net = ::zephyr::network::network::dummy;
#endif

inline bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

PvdFile::PvdFile()
    : variables(), filter(), hex_only(false), m_counter(0)
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
    : variables(), filter(), hex_only(false), m_counter(0)
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
    ofs << "    <Collection>\n";

    m_pos = ofs.tellp();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    m_open = true;
}

void PvdFile::save(Storage& elements, double timestep) {
    if (!m_open) {
        throw std::runtime_error("PvdFile::save() error: You need to open PvdFile");
    }
    VtuFile::save(get_filename(), elements, variables, hex_only, filter);
    update_pvd(elements, timestep);
}

std::string PvdFile::get_filename() const {
    std::string filename = m_fullname + "_" + std::to_string(m_counter);
#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    if (&net != &dummy_net) {
        filename += ".pt" + std::to_string(net.rank());
    }
#endif
    filename += ".vtu";
    return filename;
}

void PvdFile::update_pvd(Storage& elements, double timestep) {
    if (!m_open) {
        throw std::runtime_error("You need to open PvdFile");
    }

#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    std::VECTOR<size_t> part_sizes;
    if (&net != &dummy_net) {
        // Параллельная запись нескольких файлов
        part_sizes = net.all_gather(elements.size());

        if (!net.is_master()) {
            ++m_counter;
            return;
        }
    }
#endif

    std::fstream ofs;
    ofs.open(m_fullname + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << m_fullname << ".pvd\n";
    }

    ofs.seekg(m_pos, std::ios::beg);

    ofs << std::scientific << std::setprecision(15);
#ifdef ZEPHYR_ENABLE_DISTRIBUTED
    if (&net != &dummy_net) {
        for (size_t rank = 0; rank < net.size(); ++rank) {
            if (part_sizes[rank] > 0) {
                ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"" << rank << "\" file=\""
                    << m_filename << "_" << m_counter << ".pt" << rank << ".vtu" << "\"/>\n";
            }
        }
    } else {
        ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << m_filename << "_" << m_counter << ".vtu" << "\"/>\n";
    }
#else
    ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << m_filename << "_" << m_counter << ".vtu" << "\"/>\n";
#endif

    m_pos = ofs.tellg();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    ++m_counter;
}

} // namespace io
} // namespace zephyr
