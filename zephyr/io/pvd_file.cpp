#include <iomanip>
#include <fstream>
#include <filesystem>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/json.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/vtu_file.h>

#include <zephyr/mesh/cell.h>
#include <zephyr/mesh/mesh.h>

namespace zephyr::io {

using utils::mpi;

PvdFile::PvdFile()
    : open_(false), counter_(0)
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
        options.polyhedral = config["polyhedral"].as<bool>();
    }
    if (config["unique_nodes"]) {
        options.unique_nodes = config["unique_nodes"].as<bool>();
    }

    open(filename, directory);
}

void PvdFile::open(std::string_view filename) {
    open(std::string(filename), default_dir, !mpi::single());
}

void PvdFile::open(std::string_view filename, bool distributed) {
    open(std::string(filename), default_dir, mpi::single() ? false : distributed);
}

void PvdFile::open(std::string_view filename, std::string_view directory) {
    open(std::string(filename), std::string(directory), !mpi::single());
}

void PvdFile::open(std::string_view filename, std::string_view input_dir, bool distributed) {
    namespace fs = std::filesystem;

    if (open_) {
        return;
    }

    fs::path directory = fs::current_path();
    if (!input_dir.empty()) {
        fs::path dir = input_dir;
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

    filename_ = filename;
    if (filename.size() > 4) {
        if (filename.substr(filename.size() - 4) == ".pvd") {
            filename_ = filename.substr(filename.size() - 4);
        }
    }
    fullname_ = (directory / filename).string();

    // Мастер-процесс пишет заголовок PVD
    distributed_ = mpi::single() ? false : distributed;
    if (distributed_ && !mpi::master()) {
        return;
    }

    /// Откроем файл и запишем заголовок
    std::ofstream ofs;
    ofs.open(fullname_ + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Warning: Cannot open .pvd file " << fullname_ << ".pvd\n";
        return;
    }

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" + byteorder() + "\">\n";
    ofs << "    <Collection>" << std::endl;

    pos_ = ofs.tellp();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    open_ = true;
}

void PvdFile::save(mesh::Mesh& mesh, double timestep) {
    VtuFile::save(get_filename(), mesh, variables, options);
    update_pvd(timestep);
}

void PvdFile::save(mesh::RawCells& elements, double timestep) {
    VtuFile::save(get_filename(), elements, variables, options);
    update_pvd(timestep);
}

std::string PvdFile::get_filename() const {
    std::string filename = fullname_ + "_" + std::to_string(counter_);

    if (distributed_) {
        filename += ".pt" + mpi::srank();
    }

    filename += ".vtu";
    return filename;
}

void PvdFile::update_pvd(double timestep) {
    // Мастер-процесс пишет PVD
    if (distributed_ && !mpi::master()) {
        ++counter_;
        return;
    }

    if (!open_) {
        std::string message = "PvdFile::save() error: You need to open PvdFile";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    std::fstream ofs;
    ofs.open(fullname_ + ".pvd");

    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << fullname_ << ".pvd\n";
    }

    ofs.seekg(pos_, std::ios::beg);

    ofs << std::scientific << std::setprecision(15);

    if (distributed_) {
        for (int r = 0; r < mpi::size(); ++r) {
            ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"" << r << "\" file=\""
                << filename_ << "_" << counter_ << ".pt" << r << ".vtu" << "\"/>\n";
        }
    }
    else {
        ofs << "        <DataSet timestep=\"" << timestep << "\" part=\"0\" file=\""
            << filename_ << "_" << counter_ << ".vtu" << "\"/>\n";
    }

    pos_ = ofs.tellg();

    ofs << "    </Collection>\n";
    ofs << "</VTKFile>\n";

    ofs.close();

    ++counter_;
}

} // namespace zephyr::io