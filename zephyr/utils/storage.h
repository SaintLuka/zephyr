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

// Появляется в C++20
template <class T>
struct remove_cvref { using type = std::remove_cv_t<std::remove_reference_t<T>>; };

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;


/// @brief Добавлю различные полезные типы и функции
namespace soa {

/// @brief Двойной массив
template<typename T>
using vec_of_vec = std::vector<std::vector<T>>;


/// @brief Список типов (почти кортеж, но ничего не хранит)
template<typename... Ts>
struct type_list {

    /// @brief Число типов в списке
    static constexpr int size() { return sizeof...(Ts); }

    /// @brief Проверить наличие типа в списке
    template<typename T>
    static constexpr bool contain() {
        return (std::is_same_v<T, Ts> || ...);
    }

    /// @brief Индекс типа в списке параметров, или -1, если не найден
    template<typename T>
    static constexpr int index_of() {
        int idx = 0;
        bool found = ((std::is_same_v<T, Ts> ? (++idx, true) : (++idx, false)) || ...);
        return found ? (idx - 1) : -1;
    }

    /// @brief Получить тип из списка по индексу
    template<int I>
    using type = std::tuple_element_t<I, std::tuple<Ts...>>;

    /// @brief Кортеж, который для каждого параметра хранит двойной вектор
    using tuple_of_vec = std::tuple<vec_of_vec<Ts>...>;
};


template<typename List1, typename List2>
struct type_list_cat;

template<typename... Ts1, typename... Ts2>
struct type_list_cat<type_list<Ts1...>, type_list<Ts2...>> {
    using type = type_list<Ts1..., Ts2...>;
};

/// @brief Конкатенация двух списков типов
template<typename List1, typename List2>
using type_list_cat_t = typename type_list_cat<List1, List2>::type;


// ========================================================================
//          Описание базовых и пользовательских типов + функции
// ========================================================================

/// @brief Базовые типы, доступные для хранения в SoaStorage
using basic_types = type_list<
        short, unsigned short, int, unsigned int,
        long, unsigned long, long long, unsigned long long,
        float, double, long double, geom::Vector3d>;

/// @brief Допустимое число пользовательских типов разного размера. Допускаются
/// пользовательские типы, кратные 4 байтам: 4, 8, 12, 16 байт и так далее.
/// Использование типов длиннее 64 байт считается сомнительным.
static constexpr int n_custom_types = 20;

/// @brief Допустимое выравнивание пользовательского типа (кратно 4 байтам)
static constexpr int type_alignment = 4;

/// @brief Максимальный размер пользовательского типа в байтах (включительно)
static constexpr int custom_type_max_sizeof = type_alignment * n_custom_types;

// Пользовательские типы будем сохранять в массивах:
//    std::array<std::byte, 4>,
//    std::array<std::byte, 8>,
//    ...
//    std::array<std::byte, 4 * n_custom_types>
template<size_t N, typename = std::make_index_sequence<N>>
struct generate_custom_types;

template<size_t N, size_t... Is>
struct generate_custom_types<N, std::index_sequence<Is...>> {
    using type = type_list<std::array<std::byte, 4 * (Is + 1)>...>;
};

/// @brief Пользовательские типы (после редукции), доступные для хранения
/// в SoaStorage. Набор типов:
//    std::array<std::byte, 4>,
//    std::array<std::byte, 8>,
//    ...
//    std::array<std::byte, 4 * n_custom_types>
using custom_types = generate_custom_types<n_custom_types>::type;

/// @brief Является ли тип базовым? (присутствует в кортеже basic_types)
template<typename T>
static constexpr bool is_basic_type() { return basic_types::contain<T>(); }

/// @brief Является ли тип пользовательским?
template<typename T>
static constexpr bool is_custom_type() {
    // Не базовый тип
    if constexpr (is_basic_type<T>()) {
        return false;
    }
    // Размер кратен 4 байтам
    if constexpr (sizeof(T) % type_alignment != 0) {
        return false;
    }
    // Меньше определенного размера
    return sizeof(T) <= custom_type_max_sizeof;
}

/// @brief Редуцировать тип до такого, который можно хранить в SoaStorage
/// Получается void, если редукция невозможна
template<typename T>
using reduce_t = std::conditional_t<
        is_basic_type<T>(), T,
        std::conditional_t<
                is_custom_type<T>(),
                std::array<std::byte, sizeof(T)>,
                void>>;

/// @brief Тип может храниться в SoaStorage?
template<typename T>
static constexpr bool is_storable_type() {
    return !std::is_same_v<reduce_t<T>, void> && (is_basic_type<T>() || is_custom_type<T>());
}

/// @brief Статическая проверка, что тип можно хранить в SoaStorage
template<typename T>
static constexpr void assert_storable_type() {
    static_assert(is_storable_type<T>(), "Using of not storable type");
}

// Полный список допустимых типов
using storable_types = type_list_cat_t<basic_types, custom_types>;

// Кортеж векторов для всех допустимых типов
using tuple_of_values = storable_types::tuple_of_vec;

// Объявление типов
using names_t = std::vector<std::string>;

// Массив имен для всех допустимых типов
using array_of_names = std::array<names_t, storable_types::size()>;

/// @brief Индекс типа в списке
template<typename T>
static constexpr int type_index() {
    assert_storable_type<T>();
    return storable_types::index_of<reduce_t<T>>();
}

} // namespace soa


/// @brief Простая структура для доступа к полям SoaStorage
template<typename T>
struct Storable {
public:
    using type = T;

    int idx = -1;  ///< Смещение в массиве структур

    /// @brief Со смещением для векторных полей
    Storable<T> operator[](int shift) const { return {idx + shift}; }

    /// @brief MPI-тэг, используется в пересылках
    int tag() const { return (soa::type_index<T>() << 10) + idx; }

    // Проверка хранимых типов
    static_assert(soa::is_storable_type<T>(), "Storable<T>: using of not storable type T");
};


template <typename >
struct is_storable : std::false_type { };

template <typename T>
struct is_storable<Storable<T>> : std::true_type { };

/// @brief Все типы в parameter pack являются Storable<T>
template <typename... Args>
constexpr bool is_all_storable() {
    return (is_storable<std::decay_t<Args>>::value && ...);
}

/// @brief Статическая проверка, что все типы являются Storable<T>
template<typename... Args>
static constexpr void assert_all_storable() {
    static_assert(is_all_storable<Args...>(), "All types must be storable");
}


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
    bool empty() const { return m_size == 0; }

    /// @brief Получить размер массивов данных (даже если данных нет)
    size_t size() const { return m_size; }

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void resize(size_t new_size);

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void reserve(size_t new_size);

    /// @brief Изменить размеры буфферов под размеры данных
    void shrink_to_fit();

    /// @brief Добавить одно или несколько СКАЛЯРНЫХ полей
    /// @param names Имена должна быть уникальными
    /// @details Принимает произвольное число аргументов (> 0) типов, которые
    /// можно конвертировать в std::string, для всех выполняет add_one
    /// @return Единственный Storable<T> при единственном аргументе, кортеж из
    /// Storable<T> при нескольких аргументах (для structure binding).
    template<typename T, typename... Names, typename = std::enable_if_t<
            (sizeof...(Names) > 0) && (std::is_convertible_v<Names, std::string> && ...)>>
    auto add(Names&&... names) {
        soa::assert_storable_type<T>();
        if constexpr (sizeof...(Names) == 1) {
            return add_one<T>(std::forward<Names>(names)...);
        } else {
            // Короче жесть, тут сложно было добиться соблюдения порядка
            // с круглыми скобками std::tuple() или std::make_tuple() не работают.
            return std::tuple{add_one<T>(std::string(std::forward<Names>(names)))...};
        }
    }

    /// @brief Добавить новое ВЕКТОРНОЕ поле
    /// @param name Имя поля данных должно быть уникальным
    template<typename T>
    Storable<T> add(const std::string &name, int count);

    /// @brief Содержит поле с именем name типа T?
    template<typename T>
    bool contain(const std::string &name) const {
        const auto &row = get_names<T>();
        return std::find(row.cbegin(), row.cend(), name) != row.cend();
    }

    /// @brief Вытащить поле с именем name для типа T
    template<typename T>
    Storable<T> storable(const std::string &name) const {
        soa::assert_storable_type<T>();
        const auto &row = get_names<T>();
        auto it = std::find(row.cbegin(), row.cend(), name);
        if (it != row.cend()) {
            return {int(std::distance(row.cbegin(), it))};
        }
        throw std::runtime_error("Has no storable name '" + name + "'");
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    T *data(Storable<T> type) {
        soa::assert_storable_type<T>();
        if constexpr (soa::is_basic_type<T>()) {
            return get_values<T>()[type.idx].data();
        } else {
            return reinterpret_cast<T *>(get_values<T>()[type.idx].data()->data());
        }
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    const T *data(Storable<T> type) const {
        soa::assert_storable_type<T>();
        if constexpr (soa::is_basic_type<T>()) {
            return get_values<T>()[type.idx].data();
        } else {
            return reinterpret_cast<const T *>(get_values<T>()[type.idx].data()->data());
        }
    }

    /// @brief Кортеж указателей на данные (может быть пустым)
    template<typename... Args>
    auto data_tuple(const std::tuple<Args...>& vars) {
        assert_all_storable<Args...>();
        return std::apply([this](const auto&... args) {
            return std::tuple{data(args)...};
        }, vars);
    }

    /// @brief Кортеж указателей на данные (может быть пустым)
    template<typename... Args>
    auto data_tuple_1(Args&&... vars) {
        assert_all_storable<Args...>();
        return std::tuple{data(std::forward<Args>(vars))...};
    }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    auto operator()(Storable<T> var) { return data(var); }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    auto operator()(Storable<T> var) const { return data(var); }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    auto operator[](Storable<T> var) { return data(var); }

    /// @brief Извлечь поле по переменной Storable<T>
    template<typename T>
    auto operator[](Storable<T> var) const { return data(var); }

    /// @brief Скопировать все данные по индексу from
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, size_t to);

    /// @brief Скопировать все данные по индексу from
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, SoaStorage* dst, size_t to) const;

    /// @brief Поменять два массива данных местами
    template <typename T>
    void swap(const Storable<T>& type1, const Storable<T>& type2) {
        std::swap(get_values<T>()[type1.idx], get_values<T>()[type2.idx]);
    }

    /// @brief Вывести в консоль массивы данных, массивы пользовательских
    /// типов не выводятся, только имена и размеры типов данных
    void print() const;

    /// @brief Для списка переменных типа Storable<T1>, Storable<T2>, ...
    /// находятся места в существующем хранилище, возможно, с заменой уже
    /// существующих данных. Возвращается кортеж с теми же типами
    /// Storable<T1>, Storable<T2>, ... которые указывают на нужные
    /// массивы данных в хранлище.
    /// @param vars Переменные типа Storable<T1>, Storable<T2>, ...
    /// @return Кортеж переменных типа <Storable<T1>, Storable<T2>, ...>, который
    /// соответствует входным аргументам.
    template <typename... Args>
    std::tuple<remove_cvref_t<Args>...> add_replace(Args&&... vars);

public:

    // ========================================================================
    //                        Приватные поля класса
    // ========================================================================

    /// @brief Длина всех используемых массивов данных
    size_t m_size = 0;

    /// @brief Имена полей данных
    soa::array_of_names  m_names  = {};

    /// @brief Массивы данных
    soa::tuple_of_values m_values = {};

    // ========================================================================
    //                       Приватные методы класса
    // ========================================================================

    /// @brief Имена переменных для заданного типа
    template <typename T>
    soa::names_t &get_names() {
        soa::assert_storable_type<T>();
        return m_names[soa::type_index<T>()];
    }

    /// @brief Имена переменных для заданного типа
    template <typename T>
    const soa::names_t &get_names() const {
        soa::assert_storable_type<T>();
        return m_names[soa::type_index<T>()];
    }

    // Приватный доступ к данным по типу (для базовых типов)
    template<typename T>
    soa::vec_of_vec<soa::reduce_t<T>> &get_values() {
        soa::assert_storable_type<T>();
        return std::get<soa::type_index<T>()>(m_values);
    }

    // Приватный доступ к данным по типу (для базовых типов)
    template<typename T>
    const soa::vec_of_vec<soa::reduce_t<T>> &get_values() const {
        soa::assert_storable_type<T>();
        return std::get<soa::type_index<T>()>(m_values);
    }

    // Количество уже добавленных полей типа T
    template<typename T>
    int type_count() const {
        soa::assert_storable_type<T>();
        return get_values<T>().size();
    }

    /// @brief Добавить новое СКАЛЯРНОЕ поле
    /// @param name Имя поля данных должно быть уникальным
    template<typename T>
    Storable<T> add_one(const std::string &name);

    /// @brief Добавить поле данных типа T после индекса idx.
    /// Если после idx уже существует поле данных, то используется оно само,
    /// фактически происходит замена.  Если поля данных нет, то оно
    /// добавляется в конец. Полю данных назначается имя "DEF_{idx}"
    template<typename T>
    Storable<T> paste_one(int idx);
};

template<typename T>
Storable<T> SoaStorage::add_one(const std::string &name) {
    soa::assert_storable_type<T>();
    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("SoaStorage::add() error: SoaStorage contains name '" + name + "' already");
    }

    int c = type_count<T>();
    get_names<T>().push_back(name);
    get_values<T>().push_back(std::vector<soa::reduce_t<T>>(m_size));

    return {c};
}

template<typename T>
Storable<T> SoaStorage::add(const std::string &name, int count) {
    soa::assert_storable_type<T>();
    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("SoaStorage::add() error: SoaStorage contains name '" + name + "' already");
    }

    int c = type_count<T>();

    // Добавить count раз
    for (int i = 0; i < count; ++i) {
        get_names<T>().push_back(name);
        get_values<T>().push_back(std::vector<soa::reduce_t<T>>(m_size));
    }

    return {c};
}

template<typename T>
Storable<T> SoaStorage::paste_one(int idx) {
    soa::assert_storable_type<T>();

    int c = type_count<T>();
    if (idx < c) {
        get_names<T>()[idx] = "DEF_" + std::to_string(idx);
        return {idx};
    }
    else {
        get_names<T>().push_back("DEF_" + std::to_string(idx));
        get_values<T>().push_back(std::vector<soa::reduce_t<T>>(m_size));
        return {c};
    }
}

template <typename... Args>
std::tuple<remove_cvref_t<Args>...> SoaStorage::add_replace(Args&&... vars) {
    // Проверяет, что все типы
    assert_all_storable<Args...>();

    // Результирующий кортеж с типами Storable<T1>, Storable<T2>, ...
    auto result = std::tuple{std::forward<Args>(vars)...};

    // Обнуляем число полей, при добавлении поля выполняется count++
    std::array<int, soa::storable_types::size()> count = {};

    // Для каждого Storable<T> в vars выполняется paste_one
    // Было бы неплохо немного упростить запись, но пока так.
    [&]<size_t... Is> (std::index_sequence<Is...>) {
        ((std::get<Is>(result) =
                paste_one<typename std::tuple_element_t<Is, decltype(result)>::type>(
                        count[soa::type_index<typename std::tuple_element_t<Is, decltype(result)>::type>()]++)
                        ), ... );
    }(std::index_sequence_for<Args...>());

    return result;
}

} // namespace zephyr::utils