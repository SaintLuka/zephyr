#pragma once

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <tuple>
#include <type_traits>

#include <zephyr/geom/vector.h>

namespace zephyr::mesh {

// Появляется в C++20
template <class T>
struct remove_cvref { using type = std::remove_cv_t<std::remove_reference_t<T>>; };

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;


/// @brief Вспомогательные функции и типы для SoA
namespace soa {

/// @brief Массив массивов
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
    using get_type = std::tuple_element_t<I, std::tuple<Ts...>>;

    /// @brief Кортеж, каждому типу ставится в соответстиве массив массивов
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

/// @brief Базовые типы, доступные для хранения в Storage
using basic_types = type_list<
        short, unsigned short, int, unsigned int,
        long, unsigned long, long long, unsigned long long,
        float, double, long double, geom::Vector3d>;

/// @brief Допустимый базовый тип? (присутствует в кортеже basic_types)
template<typename T>
constexpr bool is_basic_type() { return basic_types::contain<T>(); }

/// @brief Допустимое число пользовательских типов разного размера. Допускаются
/// пользовательские типы, кратные 4 байтам: 4, 8, 12, 16 байт и так далее.
/// Использование типов длиннее 64 байт считается сомнительным.
constexpr int n_custom_types = 20;

/// @brief Допустимое выравнивание пользовательского типа (кратно 4 байтам)
constexpr int type_alignment = 4;

/// @brief Максимальный размер пользовательского типа в байтах (включительно)
constexpr int custom_type_max_sizeof = type_alignment * n_custom_types;

// Пользовательские типы хранятся в массивах:
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
/// в Storage. Набор типов:
//    std::array<std::byte, 4>,
//    std::array<std::byte, 8>,
//    ...
//    std::array<std::byte, 4 * n_custom_types>
using custom_types = generate_custom_types<n_custom_types>::type;

/// @brief Допустимый пользовательский тип?
template<typename T>
constexpr bool is_custom_type() {
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

/// @brief Редуцировать тип до такого, который можно непосредственно хранить
/// в Storage. Равен void, если редукция невозможна.
template<typename T>
using reduce_t = std::conditional_t<
        is_basic_type<T>(), T,
        std::conditional_t<
                is_custom_type<T>(),
                std::array<std::byte, sizeof(T)>,
                void>>;

/// @brief Тип может храниться в Storage?
template<typename T>
constexpr bool is_storable_type() {
    return !std::is_same_v<reduce_t<T>, void> && (is_basic_type<T>() || is_custom_type<T>());
}

/// @brief Статическая проверка, что тип можно хранить в Storage
template <typename T>
constexpr void assert_storable_type() {
    static_assert(is_storable_type<T>(), "Not storable type");
}

/// @brief Полный список допустимых типов
using storable_types = type_list_cat_t<basic_types, custom_types>;

/// @brief Кортеж векторов для всех допустимых типов
using tuple_of_values = storable_types::tuple_of_vec;

/// @brief Вектор имен переменных
using names_t = std::vector<std::string>;

/// @brief Массив имен для всех допустимых типов
using array_of_names = std::array<names_t, storable_types::size()>;

/// @brief Индекс типа в полном списке типов
template <typename T>
constexpr int type_index() {
    assert_storable_type<T>();
    return storable_types::index_of<reduce_t<T>>();
}

} // namespace zephyr::mesh::soa


/// @brief Структура для доступ к полям Storage
template<typename T>
struct Storable {
    /// @brief Смещение в массиве структур
    int idx = -1;

    /// @brief Дополнительное смещение для векторных полей
    Storable operator[](int shift) const { return {idx + shift}; }

    /// @brief MPI-тэг, используется в пересылках
    int tag() const { return (soa::type_index<T>() << 10) + idx; }


    using type = T; ///< Хранимый тип

    // Проверка хранимых типов
    static_assert(soa::is_storable_type<T>(), "Not storable type");
};

namespace soa {

template<typename>
struct is_storable_one : std::false_type { };

template<typename T>
struct is_storable_one<Storable<T>> : std::true_type { };

/// @brief Все типы в parameter pack являются `Storable<T>`
template<typename... Args>
constexpr bool is_storable() {
    return (is_storable_one<std::decay_t<Args>>::value && ...);
}

/// @brief Статическая проверка, что все типы в паке являются Storable<T>
template<typename... Args>
constexpr void assert_storable() {
    static_assert(is_storable<Args...>(), "All types must be storable");
}

} // namespace soa


/// @brief Класс для хранения данных в форме SoA (Structure of Arrays).
///
/// Все данные хранятся в массивах **одинаковой** длины, памятью управляют
/// стандартные `std::vector`. Примитивные типы (`int`, `size_t`, `float`,
/// `double`, `geom::Vector3d` и так далее) хранятся в виде массивов своего
/// типа, к примеру, `std::vector<int>` или `std::vector<double>`.
/// Пользовательские типы добавляются в хранилище в виде байтовых массивов
/// константного размера `std::array<std::byte, N>`. К примеру, структура
/// вида `struct Vec3i { int x, y, z; }` будет храниться в виде массива
/// `std::vector<std::array<std::byte, 12>>`. Пользовательские типы должны
/// быть кратны 4 байтам и ограничены по размеру.
///
/// Пример:
/// @code
///     Storage s(10);                   // Хранилище на 10 элементов, пока без полей
///     auto rho = s.add<double>("rho"); // Добавить поле типа double с именем "rho"
///                                      // rho имеет тип Storable<double>
///     double* rho_data = s.data(rho);  // Указатель на данные
///     s[rho][5] = 42;                  // Аналогично s.data(rho)[5] = 42;
///
///     struct Foo { double x; char z; } // Не кратно 4 байтам
///     struct Bar { double x; int  y; } // Подходящее выравнивание
///
///     auto foo = s.add<Foo>("bar");    // Error: Foo is not storable
///     auto bar = s.add<Bar>("foo");    // OK!
/// @endcode
class Storage final {
public:
    /// @brief Пустое хранилище заданного размера
    explicit Storage(size_t size = 0) : m_size(size) { }

    /// @brief Создать пустое хранилище с таким же размером данных
    Storage same() const;

    /// @brief Размер массивов данных равен нулю? (даже если данных нет)
    bool empty() const { return m_size == 0; }

    /// @brief Получить размер массивов данных (даже если данных нет)
    size_t size() const { return m_size; }

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void resize(size_t new_size);

    /// @brief Изменить размер массивов данных (даже если данных нет)
    void reserve(size_t new_size);

    /// @brief Изменить размеры буферов под размеры данных
    void shrink_to_fit();

    /// @brief Добавить одно или несколько **скалярных** полей
    /// @param names Имена полей данных, должны быть уникальными
    /// @tparam Names Произвольное число аргументов (> 0), которые можно
    /// конвертировать в `std::string`, для всех выполнится `add_one<T>`.
    /// @return Единственный `Storable<T>` при одном аргументе, кортеж из
    /// `Storable<T>` для нескольких аргументов (для structure binding).
    template<typename T, typename... Names, typename = std::enable_if_t<
            (sizeof...(Names) > 0) && (std::is_convertible_v<Names, std::string> && ...)>>
    auto add(Names&&... names) {
        soa::assert_storable_type<T>();
        if constexpr (sizeof...(Names) == 1) {
            return add_one<T>(std::forward<Names>(names)...);
        } else {
            // Тут сложно было добиться соблюдения порядка, варианты
            // с круглыми скобками std::tuple() или std::make_tuple() не работают.
            return std::tuple{add_one<T>(std::string(std::forward<Names>(names)))...};
        }
    }

    /// @brief Добавить новое **векторное** поле
    /// @param name Имя поля данных должно быть уникальным
    template<typename T>
    Storable<T> add(const std::string &name, int count);

    /// @brief Содержит поле с именем name типа T?
    template<typename T>
    bool contain(const std::string &name) const {
        const auto &row = get_names<T>();
        return std::find(row.cbegin(), row.cend(), name) != row.cend();
    }

    /// @brief Извлечь поле с именем name для типа T
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

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    T *data(Storable<T> type) {
        soa::assert_storable_type<T>();
        if constexpr (soa::is_basic_type<T>()) {
            return get_values<T>()[type.idx].data();
        } else {
            return reinterpret_cast<T *>(get_values<T>()[type.idx].data()->data());
        }
    }

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    const T *data(Storable<T> type) const {
        soa::assert_storable_type<T>();
        if constexpr (soa::is_basic_type<T>()) {
            return get_values<T>()[type.idx].data();
        } else {
            return reinterpret_cast<const T *>(get_values<T>()[type.idx].data()->data());
        }
    }

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    auto operator()(Storable<T> var) { return data(var); }

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    auto operator()(Storable<T> var) const { return data(var); }

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    auto operator[](Storable<T> var) { return data(var); }

    /// @brief Получить указатель на данные по переменной `Storable<T>`
    template<typename T>
    auto operator[](Storable<T> var) const { return data(var); }

    /// @brief Получить кортеж указателей на данные (может быть пустым)
    /// @tparam Args Произвольное число аргументов, каждый из которых имеет
    /// тип `Storable<T>`
    template<typename... Args>
    auto data_tuple(const std::tuple<Args...>& vars) {
        soa::assert_storable<Args...>();
        return std::apply([this](const auto&... args) {
            return std::tuple{data(args)...};
        }, vars);
    }

    /// @brief Скопировать данные с индекса from на индекс to
    void copy_data(size_t from, size_t to);

    /// @brief Скопировать все данные с индекса from текущего хранилища
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, Storage* dst, size_t to) const;

    /// @brief Поменять два массива данных местами
    template <typename T>
    void swap(const Storable<T>& type1, const Storable<T>& type2) {
        std::swap(get_values<T>()[type1.idx], get_values<T>()[type2.idx]);
    }

    /// @brief Вывести в консоль массивы данных. Массивы пользовательских
    /// типов не выводятся, только имена и размеры типов данных
    void print() const;

    /// @brief Для списка переменных типа `Storable<T1>, Storable<T2>, ...`
    /// находятся места в существующем хранилище, возможно, с заменой уже
    /// существующих данных. Возвращается кортеж с теми же типами
    /// `Storable<T1>, Storable<T2>, ...` которые указывают на нужные
    /// массивы данных в хранилище.
    /// @param vars Переменные типа `Storable<T1>, Storable<T2>, ...`
    /// @return Кортеж переменных типа `<Storable<T1>, Storable<T2>, ...>`, который
    /// соответствует входным аргументам.
    template <typename... Args>
    std::tuple<remove_cvref_t<Args>...> add_replace(Args&&... vars);

public:

    /// @brief Длина всех используемых массивов данных
    size_t m_size = 0;

    /// @brief Имена полей данных
    soa::array_of_names  m_names  = {};

    /// @brief Массивы данных
    soa::tuple_of_values m_values = {};

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

    /// @brief Приватный доступ к данным по типу
    template<typename T>
    soa::vec_of_vec<soa::reduce_t<T>> &get_values() {
        soa::assert_storable_type<T>();
        return std::get<soa::type_index<T>()>(m_values);
    }

    /// @brief Приватный доступ к данным по типу
    template<typename T>
    const soa::vec_of_vec<soa::reduce_t<T>> &get_values() const {
        soa::assert_storable_type<T>();
        return std::get<soa::type_index<T>()>(m_values);
    }

    /// @brief Количество добавленных полей типа T
    template<typename T>
    int type_count() const {
        soa::assert_storable_type<T>();
        return get_values<T>().size();
    }

    /// @brief Добавить новое **скалярное* поле
    /// @param name Имя поля данных должно быть уникальным
    template<typename T>
    Storable<T> add_one(const std::string &name);

    /// @brief Добавить поле данных типа T после индекса idx. Если после idx
    /// уже существует поле данных, то используется оно само, фактически
    /// происходит замена.  Если поля данных нет, то оно добавляется в конец.
    /// Полю данных назначается имя "DEF_{idx}".
    template<typename T>
    Storable<T> paste_one(int idx);
};

template<typename T>
Storable<T> Storage::add_one(const std::string &name) {
    soa::assert_storable_type<T>();
    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("Storage::add() error: Storage contains name '" + name + "' already");
    }

    int c = type_count<T>();
    get_names<T>().push_back(name);
    get_values<T>().push_back(std::vector<soa::reduce_t<T>>(m_size));

    return {c};
}

template<typename T>
Storable<T> Storage::add(const std::string &name, int count) {
    soa::assert_storable_type<T>();
    // Проверка на наличие массива данных с таким именем
    if (contain<T>(name)) {
        throw std::runtime_error("Storage::add() error: Storage contains name '" + name + "' already");
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
Storable<T> Storage::paste_one(int idx) {
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
std::tuple<remove_cvref_t<Args>...> Storage::add_replace(Args&&... vars) {
    // Проверяет, что все типы
    soa::assert_storable<Args...>();

    // Результирующий кортеж с типами Storable<T1>, Storable<T2>, ...
    auto result = std::tuple{std::forward<Args>(vars)...};

    // Обнуляем число полей, при добавлении поля выполняется count++
    std::array<int, soa::storable_types::size()> count = {};

    // Для каждого Storable<T> в vars выполняется paste_one.
    // Было бы неплохо немного упростить запись, но пока так.
    [&]<size_t... Is> (std::index_sequence<Is...>) {
        ((std::get<Is>(result) =
                paste_one<typename std::tuple_element_t<Is, decltype(result)>::type>(
                        count[soa::type_index<typename std::tuple_element_t<Is, decltype(result)>::type>()]++)
                        ), ... );
    }(std::index_sequence_for<Args...>());

    return result;
}

} // namespace zephyr::mesh