#include <cassert>
#include <charconv>

#include <zephyr/io/variable.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_node.h>

namespace zephyr::io {

using mesh::EuCell;
using mesh::EuNode;

template <typename T>
T& buffer(void* buff, int idx = 0) {
    return static_cast<T*>(buff)[idx];
}

inline bool contain(std::string_view name, const char* substr) {
    return name.find(substr) != std::string_view::npos;
}

inline int get_count(std::string_view sv) {
    auto beg = sv.find('[');
    auto end = sv.find(']', beg);
    if (beg != std::string_view::npos &&
        end != std::string_view::npos &&
        end > beg) {

        std::string_view str_num = sv.substr(beg + 1, end - beg - 1);

        int value = 0;
        auto [ptr, ec] = std::from_chars(str_num.data(), str_num.data() + str_num.size(), value);

        if (ec == std::errc()) {
            return value;
        }
    }
    return -1;
}

Variable::Variable(std::string_view name)
    : name_(name), n_components_(1) {

    // Ищет число в подстроке по типу "faces[4]"
    int n_comp = get_count(name_);

    if (name == "rank") {
        type_ = VtkType::Int32;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int32_t>(out) = cell.rank();
        };
    }
    else if (name == "index") {
        type_ = VtkType::Int32;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int32_t>(out) = cell.index();
        };
    }
    else if (name == "level") {
        type_ = VtkType::Int8;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int8_t>(out) = cell.level();
        };
    }
    else if (name == "next") {
        type_ = VtkType::Int32;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int32_t>(out) = cell.next();
        };
    }
    else if (name == "flag") {
        type_ = VtkType::Int8;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int8_t>(out) = cell.flag();
        };
    }
    else if (name == "b_idx") {
        type_ = VtkType::Int32;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int32_t>(out) = cell.b_idx();
        };
    }
    else if (name == "z_idx") {
        type_ = VtkType::Int32;
        write_ = [](const EuCell& cell, void *out) {
            buffer<int32_t>(out) = cell.z_idx();
        };
    }
    else if (contain(name, "face.rank") && n_comp > 0) {
        name_ = "face.rank";
        type_ = VtkType::Int8;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const int n_faces = std::min(n_comp, cell.face_count());
            for (int i = 0; i < n_faces; ++i) {
                buffer<int8_t>(out, i) = cell.face(i).adj_rank();
            }
            for (int i = n_faces; i < n_comp; ++i) {
                buffer<int8_t>(out, i) = int8_t{-42};
            }
        };
    }
    else if (contain(name, "face.index") && n_comp > 0) {
        name_ = "face.index";
        type_ = VtkType::Int32;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const int n_faces = std::min(n_comp, cell.face_count());
            for (int i = 0; i < n_faces; ++i) {
                buffer<int32_t>(out, i) = cell.face(i).adj_index();
            }
            for (int i = n_faces; i < n_comp; ++i) {
                buffer<int32_t>(out, i) = int32_t{-42};
            }
        };
    }
    else if (contain(name, "face.ghost") && n_comp > 0) {
        name_ = "face.ghost";
        type_ = VtkType::Int32;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const int n_faces = std::min(n_comp, cell.face_count());
            for (int i = 0; i < n_faces; ++i) {
                buffer<int32_t>(out, i) = cell.face(i).adj_ghost();
            }
            for (int i = n_faces; i < n_comp; ++i) {
                buffer<int32_t>(out, i) = int32_t{-42};
            }
        };
    }
    else if (contain(name, "face.boundary") && n_comp > 0) {
        name_ = "face.boundary";
        type_ = VtkType::Int8;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const int n_faces = std::min(n_comp, cell.face_count());
            for (int i = 0; i < n_faces; ++i) {
                buffer<int8_t>(out, i) = static_cast<int8_t>(cell.face(i).flag());
            }
            for (int i = n_faces; i < n_comp; ++i) {
                buffer<int8_t>(out, i) = int8_t{-42};
            }
        };
    }
    else if (contain(name, "face.rotation") && n_comp > 0) {
        name_ = "face.rotation";
        type_ = VtkType::Int8;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const int n_faces = std::min(n_comp, cell.face_count());
            for (int i = 0; i < n_faces; ++i) {
                buffer<int8_t>(out, i) = static_cast<int8_t>(cell.face(i).rotation());
            }
            for (int i = n_faces; i < n_comp; ++i) {
                buffer<int8_t>(out, i) = int8_t{-42};
            }
        };
    }
    else if (name == "coords" || name == "center") {
        type_ = VtkType::Float32;
        n_components_ = 3;
        write_ = [](const EuCell& cell, void *out) {
            buffer<float>(out, 0) = static_cast<float>(cell.center().x());
            buffer<float>(out, 1) = static_cast<float>(cell.center().y());
            buffer<float>(out, 2) = static_cast<float>(cell.center().z());
        };
    }
    else if (contain(name, "vert.rank") && n_comp > 0) {
        name_ = "vert.rank";
        type_ = VtkType::Int8;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const mesh::AmrCells& cells = cell.cells();
            int n_nodes = std::min(n_comp, cell.node_count());
            if (!cells.has_nodes()) n_nodes = 0;

            for (int i = 0; i < n_nodes; ++i) {
                buffer<int8_t>(out, i) = static_cast<int8_t>(cell.node_rank(i));
            }
            for (int i = n_nodes; i < n_comp; ++i) {
                buffer<int8_t>(out, i) = int8_t{-42};
            }
        };
    }
    else if (contain(name, "vert.index") && n_comp > 0) {
        name_ = "vert.index";
        type_ = VtkType::Int32;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const mesh::AmrCells& cells = cell.cells();
            int n_nodes = std::min(n_comp, cell.node_count());
            if (!cells.has_nodes()) n_nodes = 0;

            for (int i = 0; i < n_nodes; ++i) {
                buffer<int32_t>(out, i) = cell.node_index(i);
            }
            for (int i = n_nodes; i < n_comp; ++i) {
                buffer<int32_t>(out, i) = int32_t{-42};
            }
        };
    }
    else if (contain(name, "vert.ghost") && n_comp > 0) {
        name_ = "vert.ghost";
        type_ = VtkType::Int32;
        n_components_ = n_comp;
        write_ = [n_comp](const EuCell& cell, void *out) {
            const mesh::AmrCells& cells = cell.cells();
            int n_nodes = std::min(n_comp, cell.node_count());
            if (!cells.has_nodes()) n_nodes = 0;

            for (int i = 0; i < n_nodes; ++i) {
                buffer<int32_t>(out, i) = cell.node_ghost(i);
            }
            for (int i = n_nodes; i < n_comp; ++i) {
                buffer<int32_t>(out, i) = int32_t{-42};
            }
        };
    }
    else {
        throw std::runtime_error("Unknown variable '" + std::string(name) + "'");
    }
}

bool Variable::cell_data() const {
    return std::holds_alternative<WriteCell<void>>(write_);
}

bool Variable::node_data() const {
    return std::holds_alternative<WriteNode<void>>(write_);
}

void Variable::write(EuCell& cell, void* out) const {
    z_assert(std::get<WriteCell<void>>(write_) != nullptr, "Variable::write: nullptr function");
    std::get<WriteCell<void>>(write_)(cell, out);
}

void Variable::write(EuNode& node, void* out) const {
    z_assert(std::get<WriteNode<void>>(write_) != nullptr, "Variable::write: nullptr function");
    std::get<WriteNode<void>>(write_)(node, out);
}

} // namespace zephyr::io