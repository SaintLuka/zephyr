#include <cstring>

#include <zephyr/utils/storage.h>

namespace zephyr::utils {

// Название базового типа
template <typename T>
std::string basic_type_name() {
    if constexpr (std::is_same_v<T, short>) {
        return "short";
    } else if constexpr (std::is_same_v<T, unsigned short>) {
        return "unsigned short";
    } else if constexpr (std::is_same_v<T, int>) {
        return "int";
    } else if constexpr (std::is_same_v<T, unsigned int>) {
        return "unsigned int";
    } else if constexpr (std::is_same_v<T, long>) {
        return "long";
    } else if constexpr (std::is_same_v<T, unsigned long>) {
        return "unsigned long";
    } else if constexpr (std::is_same_v<T, long long>) {
        return "long long";
    } else if constexpr (std::is_same_v<T, unsigned long long>) {
        return "unsigned long long";
    } else if constexpr (std::is_same_v<T, float>) {
        return "float";
    } else if constexpr (std::is_same_v<T, double>) {
        return "double";
    } else if constexpr (std::is_same_v<T, long double>) {
        return "long double";
    } else if constexpr (std::is_same_v<T, geom::Vector3d>) {
        return "Vector3d";
    } else {
        return "Unknown basic type";
    }
}

// Выполняет resize для полей данных, нужно для копии сигнатуры SoaStorage
template <int I = 0, typename... Ts>
void resize_by_example(const std::tuple<Ts...>& old_values,
                       std::tuple<Ts...>& new_values) {
    std::get<I>(new_values).resize(std::get<I>(old_values).size());
    if constexpr (I < sizeof...(Ts) - 1) {
        resize_by_example<I + 1>(old_values, new_values);
    }
}

SoaStorage SoaStorage::same() const {
    SoaStorage new_s;
    new_s.m_names = m_names;
    resize_by_example(m_values, new_s.m_values);
    return new_s;
}

template <size_t I = 0, typename... Ts>
void resize_base(std::tuple<Ts...>& values, size_t new_size) {
    if constexpr (I < sizeof...(Ts)) {
        for (auto &vec: std::get<I>(values)) {
            vec.resize(new_size);
        }
        resize_base<I + 1>(values, new_size);
    }
}

void SoaStorage::resize(size_t new_size) {
    m_size = new_size;
    resize_base(m_values, new_size);
}

template <size_t I = 0, typename... Ts>
void reserve_base(std::tuple<Ts...>& values, size_t new_size) {
    if constexpr (I < sizeof...(Ts)) {
        for (auto &vec: std::get<I>(values)) {
            vec.reserve(new_size);
        }
        reserve_base<I + 1>(values, new_size);
    }
}

void SoaStorage::reserve(size_t new_size) {
    reserve_base(m_values, new_size);
}

template <size_t I = 0, typename... Ts>
void shrink_base(std::tuple<Ts...>& values) {
    if constexpr (I < sizeof...(Ts)) {
        for (auto &vec: std::get<I>(values)) {
            vec.shrink_to_fit();
        }
        shrink_base<I + 1>(values);
    }
}

void SoaStorage::shrink_to_fit() {
    shrink_base(m_values);
}

// Рекурсивная печать, возвращает количество напечатаных полей
template <int I = 0, typename... Ts>
int print_base(const std::array<soa::names_t, sizeof...(Ts)>& names, const std::tuple<Ts...>& values) {
    // База рекурсии
    if constexpr (I < sizeof...(Ts)) {
        int count = 0;
        if (!std::get<I>(values).empty()) {
            // разворачиваем двойной массив, чтобы узнать тип аргумента
            using vec_of_vec_T = std::tuple_element_t<I, std::tuple<Ts...>>;
            using vec_of_T = typename vec_of_vec_T::value_type;
            using elem_type = typename vec_of_T::value_type;

            // Запись примитивных величин (есть вывод cout)
            if constexpr (soa::is_basic_type<elem_type>()) {
                int K = std::get<I>(values).size();
                std::cout << "  " << basic_type_name<elem_type>() << " arrays:\n";
                for (int k = 0; k < K; ++k) {
                    auto &vec = std::get<I>(values).at(k);
                    auto &name = names[I].at(k);

                    std::cout << "    '" << name << "': [";
                    if (vec.empty()) {
                        std::cout << "]\n";
                    } else {
                        if constexpr (I == soa::type_index<geom::Vector3d>()) {
                            for (size_t i = 0; i < vec.size() - 1; ++i) {
                                std::cout << "{" << vec[i].transpose() << "}, ";
                            }
                            std::cout << "{" << vec.back().transpose() << "}]\n";
                        } else {
                            for (size_t i = 0; i < vec.size() - 1; ++i) {
                                std::cout << vec[i] << ", ";
                            }
                            std::cout << vec.back() << "]\n";
                        }
                    }
                }
            }
            else {
                std::cout << "  " << sizeof(soa::storable_types::type<I>) << " byte type arrays:\n    ";
                for (int j = 0; j < names[I].size() - 1; ++j) {
                    std::cout << "'" << names[I][j] << "', ";
                }
                std::cout << "'" << names[I].back() << "'\n";
            }

            ++count;
        }
        return print_base<I + 1>(names, values) + count;
    }
    return 0;
}

void SoaStorage::print() const {
    if (empty()) {
        std::cout << "SoaStorage is empty\n";
    }
    else {
        std::cout << "SoaStorage contains " << size() << " elements\n";
    }
    //int count = 0;
    int count = print_base(m_names, m_values);

    if (count < 1) {
        std::cout << "  Has no fields yet\n";
    }
    std::cout << "\n";
}


void SoaStorage::copy_data(size_t from, size_t to) {
    copy_data(from, this, to);
}

// Скопировать для одного базового типа значения
template <typename T>
static void copy_one(const soa::vec_of_vec<T>& src, size_t from,
                     soa::vec_of_vec<T>& dst, size_t to) {
    if (!src.empty() && src.size() == dst.size()) {
        for (size_t k = 0; k < src.size(); ++k) {
            dst[k][to] = src[k][from];
        }
    }
}

// Скопировать для всех базовых типов (рекурсивно)
template <int I = 0, typename... Ts>
static void copy_values(const std::tuple<Ts...>& src, size_t from,
                        std::tuple<Ts...>& dst, size_t to) {
    copy_one(std::get<I>(src), from, std::get<I>(dst), to);
    if constexpr (I < sizeof...(Ts) - 1) {
        copy_values<I + 1>(src, from, dst, to);
    }
}

void SoaStorage::copy_data(size_t from, SoaStorage* dst, size_t to) const {
    copy_values(m_values, from, dst->m_values, to);
}

} // namespace zephyr::utils