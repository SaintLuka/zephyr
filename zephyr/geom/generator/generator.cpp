#include <iostream>

#include <zephyr/utils/json.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

namespace zephyr::geom {

Generator::Generator(const std::string &type)
        : m_name(type) { }

Generator::Ptr Generator::create(const Json& config) {
    using namespace generator;

    if (!config["type"]) {
        throw std::runtime_error("EuMesh config doesn't contain 'type'");
    }

    std::string type = config["type"].as<std::string>();
    if (type == "strip") {
        return Strip::create(config);
    }
    else if (type == "rectangle") {
        return Rectangle::create(config);
    } else if (type == "cuboid") {
        return Cuboid::create(config);
    } else if (type == "sector") {
        return Sector::create(config);
    } else if (type.find("collection") != std::string::npos) {
        // Block Structured from collection
        int idx = type.find('.') + 1;
        std::string name = type.substr(idx);

        std::transform(name.begin(), name.end(), name.begin(), ::tolower);

        if (name == "plane-with-hole") {
            return collection::PlaneWithHole::create(config);
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

const std::string &Generator::name() const {
    return m_name;
}

Box Generator::bbox() const {
    std::string message = "Volume is not defined for the mesh generator.";
    std::cerr << message << "\n";
    throw std::runtime_error(message);
}


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