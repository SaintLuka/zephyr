#include <zephyr/mesh/raw/distributor.h>
#include <zephyr/mesh/cell.h>

namespace zephyr::mesh {

Distributor::Distributor() {
    split = [](const Cell &parent, Children &children) {};
    merge = [](const Children &children, Cell &parent) {};
}

Distributor Distributor::empty() {
    return {};
}

Distributor Distributor::simple() {
    Distributor ds;
    ds.split = [](const Cell &parent, Children &children) {
        for (auto child: children) {
            parent.copy_data_to(child);
        }
    };
    ds.merge = [](const Children &children, Cell &parent) {
        children[0].copy_data_to(parent);
    };
    return ds;
}

Distributor Distributor::initializer(std::function<void(Cell&)> func) {
    Distributor distr;
    distr.split = [func](const Cell& parent, Children& children) {
        for (auto child: children) func(child);
    };
    distr.merge = [func](const Children& children, Cell& parent) {
        func(parent);
    };
    return distr;
}

} // namespace zephyr::mesh