#include <cstring>

#include <zephyr/io/variable.h>

namespace zephyr { namespace io {

Variable::Variable(const char* name)
    : m_name(name) {

    if (!std::strcmp(name, "rank")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.rank();
        };
    }
    else if (!std::strcmp(name, "index")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.index();
        };
    }
    else if (!std::strcmp(name, "level")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.level();
        };
    }
    else if (!std::strcmp(name, "next")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.next();
        };
    }
    else if (!std::strcmp(name, "flag")) {
        m_type = VtkType::Int8;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int8_t *) arg;
            out[0] = cell.flag();
        };
    }
    else if (!std::strcmp(name, "b_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.b_idx();
        };
    }
    else if (!std::strcmp(name, "z_idx")) {
        m_type = VtkType::Int32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            out[0] = cell.z_idx();
        };
    }
    else if (!std::strcmp(name, "face.rank")) {
        m_type = VtkType::UInt8;
        m_n_components = geom::Faces::max_size;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int8_t *) arg;
            for (int i = 0; i < geom::Faces::max_size; ++i) {
                out[i] = cell.geom().faces[i].adjacent.rank;
            }
        };
    }
    else if (!std::strcmp(name, "face.index")) {
        m_type = VtkType::Int32;
        m_n_components = geom::Faces::max_size;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int32_t *) arg;
            for (int i = 0; i < geom::Faces::max_size; ++i) {
                out[i] = cell.geom().faces[i].adjacent.index;
            }
        };
    }
    else if (!std::strcmp(name, "face.boundary")) {
        m_type = VtkType::Int8;
        m_n_components = geom::Faces::max_size;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (int8_t *) arg;
            for (int i = 0; i < geom::Faces::max_size; ++i) {
                out[i] = cell.geom().faces[i].boundary;
            }
        };
    }
    else if (!std::strcmp(name, "coords") || !std::strcmp(name, "center")) {
        m_type = VtkType::Float32;
        m_n_components = 3;

        m_function = [](Storage::Item cell, void *arg) {
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

void Variable::write(Storage::Item elem, void* out) const {
    m_function(elem, out);
}

std::string Variable::name() const {
    return m_name;
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

} // namespace io
} // namespace zephyr