#include <zephyr/mesh/distributor.h>

namespace zephyr { namespace mesh {

Distributor::Distributor() {
    split2D = [](Storage::Item parent, const std::array<Storage::Item, 4> & children) { };
    split3D = [](Storage::Item parent, const std::array<Storage::Item, 8> & children) { };

    merge2D = [](const std::array<Storage::Item, 4> & children, Storage::Item parent) { };
    merge3D = [](const std::array<Storage::Item, 8> & children, Storage::Item parent) { };
}

Distributor Distributor::empty() {
    return Distributor();
}

Distributor Distributor::simple() {
    Distributor ds;
    ds.split2D = [](Storage::Item parent, const std::array<Storage::Item, 4> & children) {
        int ds = parent.datasize();
        for (auto& child: children) {
            std::memcpy(child.data_ptr(), parent.data_ptr(), ds);
        }
    };
    ds.split3D = [](Storage::Item parent, const std::array<Storage::Item, 8> & children) {
        int ds = parent.datasize();
        for (auto& child: children) {
            std::memcpy(child.data_ptr(), parent.data_ptr(), ds);
        }
    };

    ds.merge2D = [](const std::array<Storage::Item, 4> & children, Storage::Item parent) {
        std::memcpy(parent.data_ptr(), children[0].data_ptr(), parent.datasize());
    };
    ds.merge3D = [](const std::array<Storage::Item, 8> & children, Storage::Item parent) {
        std::memcpy(parent.data_ptr(), children[0].data_ptr(), parent.datasize());
    };

    return ds;
}

} // namespace mesh
} // namespace zephyr