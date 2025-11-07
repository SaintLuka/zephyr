#include <zephyr/mesh/euler/distributor.h>
#include <zephyr/mesh/euler/eu_prim.h>

namespace zephyr::mesh {

Distributor::Distributor() {
    split = [](const EuCell &parent, Children &children) {};
    merge = [](const Children &children, EuCell &parent) {};
}

Distributor Distributor::empty() {
    return {};
}

Distributor Distributor::simple() {
    Distributor ds;
    ds.split = [](const EuCell &parent, Children &children) {
        for (auto child: children) {
            parent.copy_data_to(child);
        }
    };
    ds.merge = [](const Children &children, EuCell &parent) {
        children[0].copy_data_to(parent);
    };
    return ds;
}

Distributor Distributor::initializer(std::function<void(EuCell&)> func) {
    Distributor distr;
    distr.split = [func](const EuCell& parent, Children& children) {
        for (auto child: children) func(child);
    };
    distr.merge = [func](const Children& children, EuCell& parent) {
        func(parent);
    };
    return distr;
}

} // namespace zephyr::mesh