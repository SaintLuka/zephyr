#include <cstring>
#include <charconv>

#include <zephyr/io/variables.h>

namespace zephyr::io {

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

Variables::Variables(std::string_view name) {
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

void Variables::append(std::string_view name) {
    int n_comp = get_count(name);
    if (contain(name, "faces")) {
        if (name == "faces2D") n_comp = 8;
        if (name == "faces3D") n_comp = 24;
        if (n_comp > 0) {
            std::string count = "[" + std::to_string(n_comp) + "]";

            // Здесь добавляются сложные типы данных
            list_.emplace_back("face.rank" + count);
            list_.emplace_back("face.ghost" + count);
            list_.emplace_back("face.index" + count);
            list_.emplace_back("face.boundary" + count);
            list_.emplace_back("face.rotation" + count);
        }
    }
    else if (contain(name, "verts")) {
        if (name == "verts2D") n_comp = 9;
        if (name == "verts3D") n_comp = 27;
        if (n_comp > 0) {
            std::string count = "[" + std::to_string(n_comp) + "]";

            // Здесь добавляются сложные типы данных
            list_.emplace_back("vert.rank" + count);
            list_.emplace_back("vert.ghost" + count);
            list_.emplace_back("vert.index" + count);
        }
    }
    else {
        list_.emplace_back(name);
    }
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
        list_.emplace_back(desc);
    }
}

void Variables::reset() {
    list_.clear();
}

const Variable& Variables::operator[](int i) const {
    return list_[i];
}

size_t Variables::size() const {
    return list_.size();
}

const std::vector<Variable>& Variables::list() const {
    return list_;
}

} // namespace zephyr::io