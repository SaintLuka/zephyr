#include <zephyr/mesh/generator/generator.h>
#include <zephyr/mesh/generator/cuboid.h>
#include <zephyr/mesh/generator/rectangle.h>
#include <zephyr/mesh/generator/sector.h>

#include <zephyr/mesh/generator/block.h>
#include <zephyr/mesh/generator/curve/curve.h>
#include <zephyr/mesh/generator/block_structured.h>

#include <zephyr/mesh/generator/collection/plane_with_hole.h>

namespace zephyr { namespace mesh { namespace generator {


generator::Generator::Part::Part(int n_cells, int n_parts, int rank)
: n_cells(n_cells), n_parts(n_parts), rank(rank) {
    int bin = n_cells / n_parts;
    from = bin * rank;
    to = bin * (rank + 1);
    if (rank + 1 == n_parts) {
        to = n_cells;
    }
}

int generator::Generator::Part::size() const {
    return to - from;
}

Generator::Generator(const std::string &type)
        : m_type(type), m_size(0) {}

const std::string &Generator::type() const {
    return m_type;
}

int Generator::size() const {
    return m_size;
}

Box Generator::bbox() const {
    std::string message = "Volume is not defined for the mesh generator.";
    std::cerr << message << "\n";
    throw std::runtime_error(message);
}

#ifdef ZEPHYR_ENABLE_YAML
std::unique_ptr<Generator> Generator::create(YAML::Node config) {
    if (!config["type"]) {
        throw std::runtime_error("Mesh config doesn't contain 'type'");
    }

    std::string type = config["type"].as<std::string>();
    if (type == "rectangle") {
        return std::unique_ptr<VRectangle>(new VRectangle(config));
    } else if (type == "cuboid") {
        return std::unique_ptr<Cuboid>(new Cuboid(config));
    } else if (type == "sector") {
        return std::unique_ptr<Sector>(new Sector(config));
    } else if (type.find("collection") != std::string::npos) {
        // Block Structured from collection
        int idx = type.find('.') + 1;
        std::string name = type.substr(idx);

        std::transform(name.begin(), name.end(), name.begin(), ::tolower);

        if (name == "plane-with-hole") {
            return std::unique_ptr<collection::PlaneWithHole>(
                    new collection::PlaneWithHole(config));
        }
        else {
            std::string message = "Create generator error: unknown mesh type '"
                    + name + "' from collection";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }
    }
    else {
        std::string message = "Error: Unknown mesh type '" + type + "'";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
}
#endif

void Generator::check_size() const {
    if (size() < 1) {
        std::string message = "'" + m_type + "' generator error: Please set mesh size";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (size() > 10000000) {
        std::string message =
                "'" + m_type +
                "' generator error: "
                "You are trying to create mesh that contains more than 10 million elements, "
                "I just want to keep your RAM.";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
}

void Generator::check_params() const {

}

void Generator::initialize(Storage &storage, Part part) {
    if (part.rank == 0) {
        storage.resize(part.n_cells);
        initialize(storage);
    }
    else {
        storage.resize(0);
    }
}

} // namespace generator
} // namespace mesh
} // namespace zephyr