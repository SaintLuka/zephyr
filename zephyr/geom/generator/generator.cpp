#include <iostream>
#include <format>

#include <zephyr/utils/json.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/generator.h>
#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

namespace zephyr::geom {

Generator::Generator(const std::string &name)
    : m_name(name) { }

Generator::Ptr Generator::create(const Json& config) {
    using namespace generator;

    if (!config["type"]) {
        throw std::runtime_error("Generator::create: config doesn't contain key 'type'");
    }

    auto type = config["type"].as<std::string>();
    if (type == "strip") {
        return Strip::create(config);
    }
    if (type == "rectangle") {
        return Rectangle::create(config);
    }
    if (type == "cuboid") {
        return Cuboid::create(config);
    }
    if (type == "sector") {
        return Sector::create(config);
    }
    if (type.find("collection") != std::string::npos) {
        // Block Structured from collection
        auto idx = type.find('.') + 1;
        std::string name = type.substr(idx);

        std::ranges::transform(name, name.begin(), ::tolower);

        if (name == "plane-with-hole") {
            return collection::PlaneWithHole::create(config);
        }
        throw std::runtime_error("Generator::create: unknown mesh type '" + name + "' from collection");
    }
    throw std::runtime_error("Generator::create: unknown mesh type '" + type + "'");
}

const std::string &Generator::name() const {
    return m_name;
}

Box Generator::bbox() const {
    throw std::runtime_error("Generator::bbox: volume is not defined for generator.");
}

void Generator::check_size(size_t size) const {
    if (size < 1) {
        throw std::runtime_error("Generator::check_size: '" + m_name + "' generator error, cells count is zero.");
    }
    if (size > max_grid_size) {
        throw std::runtime_error(std::format("Generator::check_size: attempt to create mesh that contains more than {} elements", max_grid_size));
    }
}

void Generator::set_axial(bool axial) {
    m_axial = axial;
}

void Generator::set_adaptive(bool adaptive) {
    m_adaptive = adaptive;
}

void Generator::set_linear(bool linear) {
    m_linear = linear;
}

} // namespace zephyr::geom