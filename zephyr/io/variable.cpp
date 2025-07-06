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
            auto out = (int32_t *) arg;
            out[0] = cell.rank();
        };
    }
    else if (!std::strcmp(name, "index")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.index();
        };
    }
    else if (!std::strcmp(name, "level")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.level();
        };
    }
    else if (!std::strcmp(name, "next")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.next();
        };
    }
    else if (!std::strcmp(name, "flag")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.flag();
        };
    }
    else if (!std::strcmp(name, "b_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.b_idx();
        };
    }
    else if (!std::strcmp(name, "z_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.z_idx();
        };
    }
    else if (!std::strcmp(name, "face2D.rank") || !std::strcmp(name, "face3D.rank")) {
        m_type = VtkType::Int8;
        int n_faces = !std::strcmp(name, "face2D.rank") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = (int8_t *) arg;
            for (int i = 0; i < n_faces; ++i) {
                out[i] = cell.face(i).neib_rank();
            }
        };
    }
    else if (!std::strcmp(name, "face2D.index") || !std::strcmp(name, "face3D.index")) {
        m_type = VtkType::Int32;
        int n_faces = !std::strcmp(name, "face2D.index") ? 8 : 24;
        m_n_components = n_faces;

        m_write = [n_faces](EuCell& cell, void *arg) {
            auto out = (int32_t *) arg;
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
            auto out = (int32_t *) arg;
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
            auto out = (int8_t *) arg;
            for (int i = 0; i < n_faces; ++i) {
                out[i] = (int8_t) cell.face(i).flag();
            }
        };
    }
    else if (!std::strcmp(name, "coords") || !std::strcmp(name, "center")) {
        m_type = VtkType::Float32;
        m_n_components = 3;

        m_write = [](EuCell& cell, void *arg) {
            auto out = (float *) arg;
            out[0] = float(cell.center().x());
            out[1] = float(cell.center().y());
            out[2] = float(cell.center().z());
        };
    }
    else {
        throw std::runtime_error("Unknown variable '" + std::string(name) + "'");
    }
}

Variable::Variable(const std::string& name)
    : Variable(name.c_str()) {

}

void Variable::write(EuCell& elem, void* out) const {
    assert(m_write != nullptr);
    m_write(elem, out);
}

std::string Variable::name() const {
    return m_name;
}

bool Variable::is_eu_cell() const {
    return m_write != nullptr;
}

VtkType Variable::type() const {
    return m_type;
}

int Variable::n_components() const {
    return m_n_components;
}

bool Variable::is_scalar() const {
    return m_n_components < 2;
}

size_t Variable::size() const {
    return m_n_components * m_type.size();
}

} // namespace zephyr::io