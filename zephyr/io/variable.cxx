#include <cstring>

#include <zephyr/io/variable.h>

namespace zephyr { namespace io {

Variable::Variable(const char* name)
    : m_name(name) {

    if (!std::strcmp(name, "base_id")) {
        m_type = VtkType::UInt32;
        m_n_components = 1;

        m_function = [](Storage::Item cell, void *arg) {
            auto out = (uint32_t *) arg;
            out[0] = cell.b_idx();
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
        throw std::runtime_error("Unknown variable");
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