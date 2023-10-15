#include <zephyr/geom/box.h>
#include <zephyr/geom/generator/generator.h>

namespace zephyr::geom {

Generator::Generator(const std::string &type)
        : m_name(type) { }

const std::string &Generator::name() const {
    return m_name;
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
            std::string message = "Create generators error: unknown mesh type '"
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
        std::string message = "'" + m_name + "' generator error: Please set mesh size";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (size() > 10000000) {
        std::string message =
                "'" + m_name +
                "' generator error: "
                "You are trying to create mesh that contains more than 10 million elements, "
                "I just want to keep your RAM.";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
}

void Generator::check_params() const {

}

} // namespace zephyr::geom