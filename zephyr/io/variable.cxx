#include <cstring>

#include <zephyr/io/variable.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/lagrange/la_mesh.h>

namespace zephyr::io {

using namespace zephyr::mesh;

Variable::Variable(const char* name)
    : m_name(name) {

    if (!std::strcmp(name, "rank")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.rank;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.rank();
        };
    }
    else if (!std::strcmp(name, "index")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.index;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.index();
        };
    }
    else if (!std::strcmp(name, "level")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.level;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.level();
        };
    }
    else if (!std::strcmp(name, "next")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.next;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.next();
        };
    }
    else if (!std::strcmp(name, "flag")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.flag;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.flag();
        };
    }
    else if (!std::strcmp(name, "b_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.b_idx;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.b_idx();
        };
    }
    else if (!std::strcmp(name, "z_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.z_idx;
        };
        m_soa_func = [](QCell& cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.z_idx();
        };
    }
    else if (!std::strcmp(name, "face.rank")) {
        m_type = VtkType::Int8;
        m_n_components = BFaces::max_count;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int8_t *) arg;
            for (int i = 0; i < BFaces::max_count; ++i) {
                out[i] = cell.faces[i].adjacent.rank;
            }
        };
    }
    else if (!std::strcmp(name, "face.alien")) {
        m_type = VtkType::Int32;
        m_n_components = BFaces::max_count;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            for (int i = 0; i < BFaces::max_count; ++i) {
                out[i] = cell.faces[i].adjacent.alien;
            }
        };
    }
    else if (!std::strcmp(name, "face.index")) {
        m_type = VtkType::Int32;
        m_n_components = BFaces::max_count;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int32_t *) arg;
            for (int i = 0; i < BFaces::max_count; ++i) {
                out[i] = cell.faces[i].adjacent.index;
            }
        };
    }
    else if (!std::strcmp(name, "face.boundary")) {
        m_type = VtkType::Int8;
        m_n_components = BFaces::max_count;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (int8_t *) arg;
            for (int i = 0; i < BFaces::max_count; ++i) {
                out[i] = int(cell.faces[i].boundary);
            }
        };
    }
    else if (!std::strcmp(name, "coords") || !std::strcmp(name, "center")) {
        m_type = VtkType::Float32;
        m_n_components = 3;

        m_amr_func = [](AmrStorage::Item& cell, void *arg) {
            auto out = (float *) arg;
            out[0] = float(cell.center.x());
            out[1] = float(cell.center.y());
            out[2] = float(cell.center.z());
        };
    }
    else {
        throw std::runtime_error("Unknown variable '" + std::string(name) + "'");
    }
}

Variable::Variable(const std::string& name)
    : Variable(name.c_str()) {

}

void Variable::write(AmrStorage::Item& elem, void* out) const {
    assert(m_amr_func != nullptr);
    m_amr_func(elem, out);
}

void Variable::write(QCell& elem, void* out) const {
    assert(m_soa_func != nullptr);
    m_soa_func(elem, out);
}

void Variable::write(CellStorage::Item& elem, void* out) const {
    assert(m_cell_func != nullptr);
    m_cell_func(elem, out);
}

void Variable::write(NodeStorage::Item& elem, void* out) const {
    assert(m_node_func != nullptr);
    m_node_func(elem, out);
}

std::string Variable::name() const {
    return m_name;
}

bool Variable::is_eu_cell() const {
    return m_amr_func != nullptr;
}

bool Variable::is_soa_cell() const {
    return m_soa_func != nullptr;
}

bool Variable::is_lag_cell() const {
    return m_cell_func != nullptr;
}

bool Variable::is_node() const {
    return m_node_func != nullptr;
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