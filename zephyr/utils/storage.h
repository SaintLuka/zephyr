#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>

#include <zephyr/geom/vector.h>

namespace zephyr::utils {

/// @brief Простая структура для доступа к полям SoaStorage
template<typename T>
struct Storable {
public:
    int idx = -1;  ///< Смещение в массиве структур

    /// @brief Со смещением для векторных полей
    Storable<T> operator[](int shift) const { return {idx + shift}; }

    /// @brief MPI-тэг, используется в пересылках
    int tag() const;
};

/// @class SoaStorage storage.h
/// @brief Класс для хранения разнородных данных в форме SoA (Structure of Arrays),
///
/// Все данные хранятся в массивах одинаковой длины, памятью управляют простые
/// стандартные std::vector. Примитивные типы (int, size_t, float, double,
/// Vector3d и так далее) хранятся в виде массивов своего типа, к примеру,
/// std::vector<int>, std::vector<double>. Сложные пользовательские типы
/// добавляются в хранилище в виде байтовых массивов std::vector<std::byte>.
/// Пользовательские типы должны быть кратны 4 байтам и ограничены по размеру.
class SoaStorage final {
public:
    /// @brief Пустое хранилище заданного размера
    explicit SoaStorage(size_t size = 0) : m_size(size) { }

    /// @brief Создать пустое хранилище с таким же размером данных
    SoaStorage same() const;

    /// @brief Размер массивов данных равен нулю? (даже если данных нет)
    inline bool empty() const { return m_size == 0; }

    /// @brief Получить размер массивов данных (даже если данных нет)
    inline size_t size() const { return m_size; }

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void resize(size_t new_size);

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void reserve(size_t new_size);

    /// @brief Изменить размеры буфферов под размеры данных
    void shrink_to_fit();

    /// @brief Добавить новое СКАЛЯРНОЕ поле
    /// @param name Имя для поля данных, должно быть уникальным среди полей
    /// данного типа, поля типа double и int могут иметь одинаковое имя, но
    /// зачем так делать?
    template<typename T>
    Storable<T> add(const std::string &name);

    /// @brief Добавить новое ВЕКТОРНОЕ поле
    /// @param name Имя для поля данных, должно быть уникальным среди полей
    /// данного типа, поля типа double и int могут иметь одинаковое имя, но
    /// зачем так делать?
    template<typename T>
    Storable<T> add(const std::string &name, int count);

    /// @brief Содержит поле с именем name типа T?
    template<typename T>
    bool contain(const std::string &name) const {
        const auto &row = names<T>();
        return std::find(row.cbegin(), row.cend(), name) != row.cend();
    }

    /// @brief Вытащить поле с именем name для типа T
    template<typename T>
    Storable<T> storable(const std::string &name) const {
        const auto &row = names<T>();
        auto it = std::find(row.cbegin(), row.cend(), name);
        if (it != row.cend()) {
            return {int(std::distance(row.cbegin(), it))};
        }
        throw std::runtime_error("Has no storable name '" + name + "'");
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    const T *data(Storable<T> type) const {
        assert_storable_type<T>();
        if constexpr (is_basic_type<T>()) {
            return basic_values<T>().at(type.idx).data();
        } else {
            return reinterpret_cast<const T *>(m_values2[custom_type_index<T>()].at(type.idx).data());
        }
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    T *data(Storable<T> type) {
        assert_storable_type<T>();
        if constexpr (is_basic_type<T>()) {
            return basic_values<T>().at(type.idx).data();
        } else {
            return reinterpret_cast<T *>(m_values2[custom_type_index<T>()].at(type.idx).data());
        }
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template <typename T>
    T* operator()(Storable<T> type) { return data(type); }

    /// @brief Извлечь поле по переменной Storable<T>
    template <typename T>
    const T* operator()(Storable<T> type) const { return data(type); }

    /// @brief Скопировать все данные по индексу from
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, size_t to);

    /// @brief Скопировать все данные по индексу from
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, SoaStorage* dst, size_t to) const;

    /// @brief Поменять два массива данных местами
    template <typename T>
    void swap(const Storable<T>& type1, const Storable<T>& type2) {
        if constexpr (is_basic_type<T>()) {
            std::swap(basic_values<T>()[type1.idx], basic_values<T>()[type2.idx]);
        }
        else {
            std::swap(custom_values<T>()[type1.idx], custom_values<T>()[type2.idx]);
        }
    }

    /// @brief Вывести в консоль массивы данных, массивы пользовательских
    /// типов не выводятся, только имена и размеры типов данных
    void print() const;

public:
    // ========================================================================
    //          Описание базовых и пользовательских типов + функции
    // ========================================================================

    // Кортеж с базовыми типами, доступными для хранения в SoaStorage
    using basic_types_t = std::tuple<
            short, unsigned short, int, unsigned int,
            long, unsigned long, long long, unsigned long long,
            float, double, long double, geom::Vector3d>;

    // Число типов в basic_types_t
    static constexpr int n_basic_types = std::tuple_size<basic_types_t>::value;

    // Индекс типа в кортеже basic_types_t или -1, если тип не найден
    template <typename T>
    static constexpr int basic_type_index() {
        int index = -1;
        [&index] <size_t... Is>(std::index_sequence<Is...>) {
                ((std::is_same_v<T, std::tuple_element_t<Is, basic_types_t>> ? (index = Is) : 0), ...);
        }
        (std::make_index_sequence<n_basic_types>());
        return index;
    }

    // Название базового типа
    template <typename T>
    static std::string basic_type_name();

    // Название базового типа по индексу в кортеже
    template<int I>
    static std::string basic_type_name();

    // Является ли тип базовым? (присутствует в кортеже basic_types_t)
    template <typename T>
    static constexpr bool is_basic_type() { return basic_type_index<T>() >= 0; }

    // Допустимое число пользовательских типов разного размера. Допускаются
    // пользовательские типы, кратные 4 байтам: 4, 8, 12, 16 байт и так далее.
    // Использование типов длиннее 64 байт считается сомнительным.
    static constexpr int n_custom_types = 20;

    // Допустимое выравнивание пользовательского типа (кратно 4 байтам)
    static constexpr int type_alignment = 4;

    // Максимальный размер пользовательского типа в байтах (включительно)
    static constexpr int custom_type_max_sizeof = type_alignment * n_custom_types;

    // Индекс пользовательского типа данных, возвращает -1, если тип нельзя
    // поместить в SoaStorage по некоторым причинам (выравнивание / размер).
    template <typename T>
    static constexpr int custom_type_index() {
        if constexpr (is_basic_type<T>()) {
            return -1;
        }
        if constexpr (sizeof(T) % type_alignment != 0) {
            return -1;
        }
        if constexpr (sizeof(T) > custom_type_max_sizeof) {
            return -1;
        }
        return sizeof(T) / type_alignment - 1;
    }

    // Размер пользовательского типа по индексу
    static constexpr int custom_type_sizeof(int idx) {
        assert(0 <= idx && idx < custom_type_max_sizeof &&
                "Индекс пользовательского типа out_of_range");
        return type_alignment * (idx + 1);
    }

    // Является допустимым пользовательским типом?
    template <typename T>
    static constexpr bool is_custom_type() { return custom_type_index<T>() >= 0; }

    // Является базовым типом
    template <typename T>
    static constexpr void assert_basic_type() {
        static_assert(is_basic_type<T>(), "В данном контексте допускается только базовый тип");
    }

    // Является пользовательским типом
    template <typename T>
    static constexpr void assert_custom_type() {
        static_assert(!is_basic_type<T>() && is_custom_type<T>(),
                      "В данном контексте допускается только пользовательский тип");
    }

    // Тип можно ли сохранить в SoaStorage
    template <typename T>
    static constexpr void assert_storable_type() {
        if constexpr (!is_basic_type<T>()) {
            static_assert(sizeof(T) % type_alignment == 0, "Пользовательский тип должен быть кратен 4 байтам");
            static_assert(sizeof(T) <= custom_type_max_sizeof, "Слишком большой пользовательский тип");
        }
    }

    // ========================================================================
    //                         Некоторые типы данных
    // ========================================================================

    // Объявление типов
    using names_t = std::vector<std::string>;

    // Двойной массив
    template<typename T>
    using vec_of_vec = std::vector<std::vector<T>>;

    template<typename... Ts>
    static auto to_vec_of_vec_tuple(const std::tuple<Ts...> &)
    -> std::tuple<vec_of_vec<Ts>...> { return {}; }

    using vec_of_vec_types = decltype(to_vec_of_vec_tuple(basic_types_t{}));

    // ========================================================================
    //                        Приватные поля класса
    // ========================================================================

    // Длина используемых массивов
    size_t m_size = 0;

    // Имена для массивов базовых типов
    std::array<names_t, n_basic_types> m_names1 = {};
    
    // Массивы массивов в кортеже для базовых типов
    vec_of_vec_types m_values1 = {};
    
    // Имена для массивов пользовательских типов
    std::array<names_t, n_custom_types> m_names2 = {};

    // Массивы массивов для пользовательских типов
    std::array<vec_of_vec<std::byte>, n_custom_types> m_values2 = {};

    // ========================================================================
    //                       Приватные методы класса
    // ========================================================================

    // Приватный доступ к именам по типу
    template<typename T>
    inline std::vector<std::string> &names() {
        assert_storable_type<T>();

        if constexpr (is_basic_type<T>()) {
            return m_names1[basic_type_index<T>()];
        } else {
            return m_names2[custom_type_index<T>()];
        }
    }

    // Приватный доступ к именам по типу
    template<typename T>
    inline const names_t &names() const {
        assert_storable_type<T>();

        if constexpr (is_basic_type<T>()) {
            return m_names1[basic_type_index<T>()];
        } else {
            return m_names2[custom_type_index<T>()];
        }
    }

    // Приватный доступ к данным по типу (для базовых типов)
    template<typename T>
    inline vec_of_vec<T> &basic_values() {
        assert_basic_type<T>();
        return std::get<vec_of_vec<T>>(m_values1);
    }

    // Приватный доступ к данным по типу (для базовых типов)
    template<typename T>
    inline const vec_of_vec<T> &basic_values() const {
        assert_basic_type<T>();
        return std::get<vec_of_vec<T>>(m_values1);
    }

    // Приватный доступ к данным по типу (для пользовательских типов)
    template<typename T>
    inline vec_of_vec<std::byte> &custom_values() {
        assert_custom_type<T>();
        return m_values2[custom_type_index<T>()];
    }

    // Приватный доступ к данным по типу (для пользовательских типов)
    template<typename T>
    inline const vec_of_vec<std::byte> &custom_values() const {
        assert_custom_type<T>();
        return m_values2[custom_type_index<T>()];
    }


    // Количество уже добавленных полей типа T
    template<typename T>
    int type_count() const {
        if constexpr (is_basic_type<T>()) {
            return std::get<vec_of_vec<T>>(m_values1).size();
        } else {
            return m_values2[custom_type_index<T>()].size();
        }
    }

    // Рекурсивная печать базовых типов, возвращает количество напечатаных полей
    template<int I = 0>
    std::enable_if_t<I < n_basic_types, int> print_basic() const {
        int count = 0;
        if (!std::get<I>(m_values1).empty()) {
            int K = std::get<I>(m_values1).size();
            std::cout << "  " << basic_type_name<I>() << " arrays:\n";
            for (int k = 0; k < K; ++k) {
                auto &vec = std::get<I>(m_values1).at(k);
                auto &name = m_names1[I].at(k);

                std::cout << "    '" << name << "': [";
                if (vec.empty()) {
                    std::cout << "]\n";
                }
                else {
                    if constexpr (I == basic_type_index<geom::Vector3d>()) {
                        for (int i = 0; i < vec.size() - 1; ++i) {
                            std::cout << "{" << vec[i].transpose() << "}, ";
                        }
                        std::cout << "{" << vec.back().transpose() << "}]\n";
                    }
                    else {
                        for (int i = 0; i < vec.size() - 1; ++i) {
                            std::cout << vec[i] << ", ";
                        }
                        std::cout << vec.back() << "]\n";
                    }
                }
            }
            ++count;
        }
        return print_basic<I + 1>() + count;
    }

    // Рекурсивная печать базовых типов, возвращает количество напечатаных полей
    template<int I = 0>
    std::enable_if_t<I >= n_basic_types, int> print_basic() const { return 0;}

    // Рекурсивный resize для базовых типов
    template<int I = 0>
    std::enable_if_t<I < n_basic_types> resize_basic(size_t new_size) {
        for (auto &vec: std::get<I>(m_values1)) {
            vec.resize(new_size);
        }
        resize_basic<I + 1>(new_size);
    }

    // Рекурсивный resize для базовых типов
    template<int I = 0>
    std::enable_if_t<I >= n_basic_types> resize_basic(size_t new_size) { }

    // Рекурсивный reserve для базовых типов
    template<int I = 0>
    std::enable_if_t<I < n_basic_types> reserve_basic(size_t new_size) {
        for (auto &vec: std::get<I>(m_values1)) {
            vec.reserve(new_size);
        }
        reserve_basic<I + 1>(new_size);
    }

    // Рекурсивный reserve для базовых типов
    template<int I = 0>
    std::enable_if_t<I >= n_basic_types> reserve_basic(size_t new_size) { }

    // Рекурсивный shrink_to_fit для базовых типов
    template<int I = 0>
    std::enable_if_t<I < n_basic_types> shrink_basic() {
        for (auto &vec: std::get<I>(m_values1)) {
            vec.shrink_to_fit();
        }
        shrink_basic<I + 1>();
    }

    // Рекурсивный reserve для базовых типов
    template<int I = 0>
    std::enable_if_t<I >= n_basic_types> shrink_basic() { }

    // Скопировать для одного базового типа значения
    template <typename T>
    static void copy_one_basic(const vec_of_vec<T>& src, size_t from,
                                     vec_of_vec<T>& dst, size_t to) {
        if (!src.empty() && src.size() == dst.size()) {
            for (size_t k = 0; k < src.size(); ++k) {
                dst[k][to] = src[k][from];
            }
        }
    }

    // Скопировать для всех базовых типов (рекурсивно)
    template <int I = 0>
    static void copy_basics(const vec_of_vec_types& src, size_t from,
                                  vec_of_vec_types& dst, size_t to) {
        copy_one_basic(std::get<I>(src), from, std::get<I>(dst), to);
        if constexpr (I < n_basic_types - 1) {
            copy_basics<I + 1>(src, from, dst, to);
        }
    }

    // Выполняет resize для полей данных, нужно для копии сигнатуры SoaStorage
    template <int I = 0>
    void resize_basic_values(SoaStorage& new_storage) const {
        std::get<I>(new_storage.m_values1).resize(std::get<I>(m_values1).size());
        if constexpr (I < n_basic_types - 1) {
            resize_basic_values<I + 1>(new_storage);
        }
    }
};

// Название базового типа
template <typename T>
std::string SoaStorage::basic_type_name() {
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

// Название хранимого типа по индексу в кортеже
template<int I>
std::string SoaStorage::basic_type_name() {
    return basic_type_name<std::tuple_element_t<I, basic_types_t>>();
}

template<typename T>
Storable<T> SoaStorage::add(const std::string &name) {
    assert_storable_type<T>();

    int c = type_count<T>();
    if constexpr (is_basic_type<T>()) {
        basic_values<T>().push_back(std::vector<T>(m_size));
    } else {
        // Пользовательские типы разворачиваются в массив байтов
        custom_values<T>().push_back(std::vector<std::byte>(m_size * sizeof(T)));
    }

    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("SoaStorage::add() error: SoaStorage contains name '" + name + "' already");
    }
    names<T>().push_back(name);

    return {c};
}

template<typename T>
Storable<T> SoaStorage::add(const std::string &name, int count) {
    assert_storable_type<T>();

    int c = type_count<T>();
    if constexpr (is_basic_type<T>()) {
        // Добавить count раз
        for (int i = 0; i < count; ++i) {
            basic_values<T>().push_back(std::vector<T>(m_size));
        }
    } else {
        // Пользовательские типы разворачиваются в массив байтов
        for (int i = 0; i < count; ++i) {
            custom_values<T>().push_back(std::vector<std::byte>(m_size * sizeof(T)));
        }
    }

    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("SoaStorage::add() error: SoaStorage contains name '" + name + "' already");
    }
    for (int i = 0; i < count; ++i) {
        names<T>().push_back(name);
    }

    return {c};
}

template <typename T>
int Storable<T>::tag() const {
    if constexpr (SoaStorage::is_basic_type<T>()) {
        // 1024 * type_index + idx
        return (SoaStorage::basic_type_index<T>() << 10) + idx;
    }
    else {
        // 1024 * 1024 * type_index + idx
        return (SoaStorage::custom_type_index<T>() << 20) + idx;
    }
}

} // namespace zephyr::utils