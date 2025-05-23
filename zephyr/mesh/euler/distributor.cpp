#include <cstring>

#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/mesh/euler/distributor.h>

namespace zephyr::mesh {

Distributor::Distributor() {
    split = [](AmrStorage::Item &parent, Children &children) {};
    merge = [](Children &children, AmrStorage::Item &parent) {};

    split_soa = [](QCell &parent, SoaChildren &children) {};
    merge_soa = [](SoaChildren &children, QCell &parent) {};
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
    
    ds.split_soa = [](QCell &parent, SoaChildren &children) {
        for (auto child: children) {
            parent.copy_data_to(child);
        }
    };

    ds.merge_soa = [](SoaChildren &children, QCell &parent) {
        children[0].copy_data_to(parent);
    };

    return ds;
}

} // namespace zephyr::mesh