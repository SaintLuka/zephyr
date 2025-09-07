#include <cstring>

#include <zephyr/io/variable.h>
#include <zephyr/mesh/euler/eu_prim.h>

namespace zephyr::io {

using mesh::EuCell;

Variable::Variable(const char* name)
    : m_name(name) {

    if (!std::strcmp(name, "rank")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            out[0] = cell.rank();
        };
    }
    else if (!std::strcmp(name, "index")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            out[0] = cell.index();
        };
    }
    else if (!std::strcmp(name, "level")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int8_t *>(arg);
            out[0] = static_cast<int8_t>(cell.level());
        };
    }
    else if (!std::strcmp(name, "next")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            out[0] = cell.next();
        };
    }
    else if (!std::strcmp(name, "flag")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int8_t *>(arg);
            out[0] = static_cast<int8_t>(cell.flag());
        };
    }
    else if (!std::strcmp(name, "b_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            out[0] = cell.b_idx();
        };
    }
    else if (!std::strcmp(name, "z_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            out[0] = cell.z_idx();
        };
    }
    else if (!std::strcmp(name, "face2D.rank") || !std::strcmp(name, "face3D.rank")) {
        m_type = VtkType::Int8;
        int n_faces = !std::strcmp(name, "face2D.rank") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = static_cast<int8_t *>(arg);
            for (int i = 0; i < n_faces; ++i) {
                out[i] = static_cast<int8_t>(cell.face(i).adj_rank());
            }
        };
    }
    else if (!std::strcmp(name, "face2D.index") || !std::strcmp(name, "face3D.index")) {
        m_type = VtkType::Int32;
        int n_faces = !std::strcmp(name, "face2D.index") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            for (int i = 0; i < n_faces; ++i) {
                out[i] = cell.face(i).adj_index();
            }
        };
    }
    else if (!std::strcmp(name, "face2D.alien") || !std::strcmp(name, "face3D.alien")) {
        m_type = VtkType::Int32;
        int n_faces = !std::strcmp(name, "face2D.alien") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = static_cast<int32_t *>(arg);
            for (int i = 0; i < n_faces; ++i) {
                out[i] = cell.face(i).adj_alien();
            }
        };
    }
    else if (!std::strcmp(name, "face2D.boundary") || !std::strcmp(name, "face3D.boundary")) {
        m_type = VtkType::Int8;
        int n_faces = !std::strcmp(name, "face2D.boundary") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = static_cast<int8_t *>(arg);
            for (int i = 0; i < n_faces; ++i) {
                out[i] = static_cast<int8_t>(cell.face(i).flag());
            }
        };
    }
    else if (!std::strcmp(name, "coords") || !std::strcmp(name, "center")) {
        m_type = VtkType::Float32;
        m_n_components = 3;

        m_write = [](EuCell& cell, void *arg) {
            auto out = static_cast<float *>(arg);
            out[0] = static_cast<float>(cell.center().x());
            out[1] = static_cast<float>(cell.center().y());
            out[2] = static_cast<float>(cell.center().z());
        };
    }
    else {
        throw std::runtime_error("Unknown variable '" + std::string(name) + "'");
    }
}

Variable::Variable(const std::string& name)
    : Variable(name.c_str()) {

}

void Variable::write(EuCell& cell, void* out) const {
    assert(m_write != nullptr);
    m_write(cell, out);
}

} // namespace zephyr::io