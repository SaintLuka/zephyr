#include <cstring>

#include <zephyr/mesh/primitives/mov_node.h>
#include <zephyr/mesh/primitives/mov_cell.h>
#include <zephyr/mesh/primitives/amr_cell.h>

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
    else {
        m_list.emplace_back(name);
    }
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