#include <zephyr/mesh/euler/distributor.h>

#include <zephyr/mesh/euler/eu_prim.h>

namespace zephyr::mesh {

Distributor::Distributor() {
    split = [](EuCell &parent, Children &children) {};
    merge = [](Children &children, EuCell &parent) {};
}

Distributor Distributor::empty() {
    return Distributor();
}

Distributor Distributor::simple() {
    Distributor ds;
    
    ds.split = [](EuCell &parent, Children &children) {
        for (auto child: children) {
            parent.copy_data_to(child);
        }
    };

    ds.merge = [](Children &children, EuCell &parent) {
        children[0].copy_data_to(parent);
    };

    return ds;
}

} // namespace zephyr::mesh