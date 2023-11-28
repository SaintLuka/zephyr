#include <cstring>

#include <zephyr/mesh/euler/distributor.h>

namespace zephyr::mesh {

Distributor::Distributor() {
    split = [](AmrStorage::Item &parent, Children &children) {};
    merge = [](Children &children, AmrStorage::Item &parent) {};
}

Distributor Distributor::empty() {
    return Distributor();
}

Distributor Distributor::simple() {
    Distributor ds;

    ds.split = [](AmrStorage::Item &parent, Children &children) {
        int ds = children.datasize();
        for (auto &child: children) {
            std::memcpy(child.data(), parent.data(), ds);
        }
    };

    ds.merge = [](Children &children, AmrStorage::Item &parent) {
        std::memcpy(parent.data(), children[0].data(), children.datasize());
    };

    return ds;
}

} // namespace zephyr::mesh