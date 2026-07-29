#include <cstring>

#include <zephyr/io/variables.h>

namespace zephyr::io {

Variables::Variables(const char* name) {
    append(name);
}

Variables::Variables(const std::string& name) {
    append(name);
}

Variables::Variables(std::initializer_list<const char*> names) {
    for (const char *name: names) {
        append(name);
    }
}

Variables::Variables(std::initializer_list<std::string> names) {
    for (const std::string& name: names) {
        append(name);
    }
}

Variables::Variables(const std::vector<const char*>& names) {
    for (const char *name: names) {
        append(name);
    }
}

Variables::Variables(const std::vector<std::string>& names) {
    for (const std::string& name: names) {
        append(name);
    }
}

void Variables::append(const char* name) {
    if (!std::strcmp(name, "faces")) {
        // Здесь добавляются сложные типы данных
        m_list.emplace_back("face.rank");
        m_list.emplace_back("face.alien");
        m_list.emplace_back("face.index");
        m_list.emplace_back("face.boundary");
    }
    else if (!std::strcmp(name, "faces2D")) {
        // Для EuMesh
        m_list.emplace_back("face2D.rank");
        m_list.emplace_back("face2D.index");
        m_list.emplace_back("face2D.alien");
        m_list.emplace_back("face2D.boundary");
        m_list.emplace_back("face2D.rotation");
    }
    else if (!std::strcmp(name, "faces3D")) {
        // Для EuMesh
        m_list.emplace_back("face3D.rank");
        m_list.emplace_back("face3D.index");
        m_list.emplace_back("face3D.alien");
        m_list.emplace_back("face3D.boundary");
        m_list.emplace_back("face3D.rotation");
    }
    else {
        if (contains_name(name)) return;
        m_list.emplace_back(name);
    }
}

bool Variables::contains_name(const char *name) const {
    const auto it = std::ranges::find_if(m_list,
        [&name](const Variable& var) -> bool {
           return var.name() == name;
        });
    if (it != m_list.end()) {
        std::cerr << "Attempt to add a variable with an existing name '" << name << "'\n";
        return true;
    }
    return false;
}

bool Variables::contains_name(const std::string& name) const {
    return contains_name(name.c_str());
}

void Variables::append(const std::string& name) {
    append(name.c_str());
}

void Variables::append(std::initializer_list<const char *> names) {
    for (auto& name: names) {
        append(name);
    }
}

void Variables::append(std::initializer_list<std::string> names) {
    for (auto& name: names) {
        append(name);
    }
}

void Variables::append(const std::vector<const char *> &names) {
    for (auto& name: names) {
        append(name);
    }
}

void Variables::append(const std::vector<std::string> &names) {
    for (auto& name: names) {
        append(name);
    }
}

void Variables::append(const Variables &variables) {
    for (auto& desc: variables.list()) {
        if (contains_name(desc.name())) return;
        m_list.emplace_back(desc);
    }
}

void Variables::reset() {
    m_list.clear();
}

const Variable& Variables::operator[](int i) const {
    return m_list[i];
}

size_t Variables::size() const {
    return m_list.size();
}

const std::vector<Variable>& Variables::list() const {
    return m_list;
}

} // namespace zephyr::io