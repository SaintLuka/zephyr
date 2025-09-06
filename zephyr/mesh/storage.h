#pragma once

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <tuple>
#include <type_traits>

#include <zephyr/utils/buffer.h>

namespace zephyr::mesh {

// Появляется в C++20
template <class T>
struct remove_cvref { using type = std::remove_cv_t<std::remove_reference_t<T>>; };

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;


/// @brief Структура для доступа к полям Storage
template<typename T>
class Storable {
public:
    using type = T; ///< Хранимый тип

    Storable() = default;

    bool operator!=(Storable<T> other) const { return index != other.index; }

    bool operator==(Storable<T> other) const { return index == other.index; }

    /// @brief MPI-тэг, используется в пересылках
    int tag() const { return 123000 + index; }

private:
    /// @brief Только Storage может создавать
    friend class Storage;

    /// @brief Создать из целочисленного типа, может только Storage
    template <typename V, typename = std::enable_if_t<std::is_integral_v<V>>>
    explicit Storable(V index) : index(static_cast<int>(index)) { }

    /// @brief Смещение в массиве структур
    int index = -1;
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
/// Для каждой переменной хранит имя и байтовый буфер, внутри буферов также
/// есть информация о размере элементов и MPI-тип для пересылок.
class Storage final {
    using Buffer = utils::Buffer;

public:
    /// @brief Пустое хранилище заданного размера
    explicit Storage(size_t size = 0) : m_size(size) { }

    /// @brief Создать пустое хранилище с такими же массивами данных
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
        static_assert(utils::is_simple_type_v<T>, "T must be simple type");
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
    Storable<T[]> add(const std::string &name, int count);

    /// @brief Содержит поле с именем name?
    bool contain(const std::string &name) const {
        for (const auto& ex_name: m_names) {
            if (ex_name == name) return true;
        }
        return false;
    }

    /// @brief Содержит и поддерживает тип?
    template<typename T>
    bool contain(const std::string &name) const {
        for (int idx = 0; idx < m_data.size(); ++idx) {
            if (m_names[idx] == name) {
                return m_data[idx].support<T>();
            }
        }
        return false;
    }


    /// @brief Извлечь поле с именем name для типа T
    template<typename T>
    Storable<T> find(const std::string &name) const {
        for (int idx = 0; idx < m_data.size(); ++idx) {
            if (m_names[idx] == name) {
                if (m_data[idx].support<T>()) { return Storable<T>(idx); }
                throw std::runtime_error("Storable \"" + name + "\" does not support type T");
            }
        }
        throw std::runtime_error("Don't have a Storable named " + name);
    }

    /// @brief Значение или span (предполагается единственный массив)
    template<typename T>
    auto get_val(size_t idx) { return m_data[0].get_val<T>(idx); }

    /// @brief Значение или span (предполагается единственный массив)
    template<typename T>
    auto get_val(size_t idx) const { return m_data[0].get_val<T>(idx); }


    /// @brief Получить ссылку на значение
    template<typename T>
    T& get_val(Storable<T> var, size_t idx) {
        return m_data[var.index].template get_val<T>(idx);
    }

    /// @brief Получить константную ссылку на значение
    template<typename T>
    const T& get_val(Storable<T> var, size_t idx) const {
        return m_data[var.index].template get_val<T>(idx);
    }

    /// @brief Получить span на изменяемый участок
    template<typename T>
    utils::span<T> get_val(Storable<T[]> var, size_t idx) {
        return m_data[var.index].template get_val<T[]>(idx);
    }

    /// @brief Получить span на константный участок
    template<typename T>
    utils::span<const T> get_val(Storable<T[]> var, size_t idx) const {
        return m_data[var.index].template get_val<T[]>(idx);
    }


    /// @brief Получить буфер, соответствующий `Storable<T>`
    template<typename T>
    Buffer& operator[](Storable<T> var) { return m_data[var.index]; }

    /// @brief Получить буфер, соответствующий `Storable<T>`
    template<typename T>
    const Buffer& operator[](Storable<T> var) const { return m_data[var.index]; }

    /// @brief Получить массив указателей на данные (может быть пустым)
    /// @tparam Args Произвольное число аргументов, каждый из которых имеет
    /// тип `Storable<T>`
    template<typename... Args>
    std::array<Buffer*, sizeof...(Args)>
    operator[](const std::tuple<Args...>& vars) {
        soa::assert_storable<Args...>();
        return std::apply([this](const auto&... args) {
            return std::array<Buffer*, sizeof...(Args)>{(m_data.data() + args.index)...};
        }, vars);
    }

    /// @brief Скопировать данные с индекса from на индекс to
    void copy_data(size_t from, size_t to);

    /// @brief Скопировать все данные с индекса from текущего хранилища
    /// в хранилище dst по индексу to.
    void copy_data(size_t from, Storage* dst, size_t to) const;

    /// @brief Поменять два массива данных местами
    template <typename T>
    void swap(Storable<T> var1, Storable<T> var2) {
        std::swap(m_data[var1.index], m_data[var2.index]);
    }

    void print_names() const;

    /// @brief Вывести в консоль массивы данных. Массивы пользовательских
    /// типов не выводятся, только имена и размеры типов данных
    template<typename T>
    void print(Storable<T> var) const;

    /// @brief Для списка переменных типа `Storable<T1>, Storable<T2>, ...`
    /// находятся места в существующем хранилище, возможно, с заменой уже
    /// существующих данных. Возвращается кортеж с теми же типами
    /// `Storable<T1>, Storable<T2>, ...` которые указывают на нужные
    /// массивы данных в хранилище.
    /// @param vars Переменные типа `Storable<T1>, Storable<T2>, ...`
    /// @return Кортеж переменных типа `<Storable<T1>, Storable<T2>, ...>`, который
    /// соответствует входным аргументам.
    template <typename... Vars>
    std::tuple<remove_cvref_t<Vars>...> add_replace(
        const Storage& s, const std::tuple<Vars...>& vars);

public:
    size_t                   m_size;  ///< Длина всех массивов данных
    std::vector<Buffer>      m_data;  ///< Массивы данных
    std::vector<std::string> m_names; ///< Массив имён данных

    /// @brief Добавить новое **скалярное* поле
    /// @param name Имя поля данных должно быть уникальным
    template<typename T>
    Storable<T> add_one(const std::string &name);
};

// ============================================================================
//                        Реализации шаблонных функций
// ============================================================================

template<typename T>
Storable<T> Storage::add_one(const std::string &name) {
    static_assert(utils::is_simple_type_v<T>, "T must be simple type");

    // Проверка на наличие массива данных с таким именем
    if (contain(name)) {
        throw std::runtime_error("Storable \"" + name + "\" already exists");
    }

    Storable<T> var{m_data.size()};
    m_data.emplace_back(Buffer::make<T>(m_size));
    m_names.emplace_back(name);
    return var;
}

template<typename T>
Storable<T[]> Storage::add(const std::string &name, int count) {
    static_assert(utils::is_simple_type_v<T>, "T must be simple type");
    // Проверка на наличие массива данных с таким именем
    if (contain(name)) {
        throw std::runtime_error("Storable \"" + name + "\" already exists");
    }

    Storable<T[]> var{m_data.size()};
    m_data.emplace_back(Buffer::make<T[]>(count, m_size));
    m_names.emplace_back(name);
    return var;
}

inline void Storage::print_names() const {
    if (m_names.empty()) { return; }
    std::cout << "Buffer names: ";
    for (size_t i = 0; i < m_names.size() - 1; ++i) {
        std::cout << m_names[i] << ", ";
    }
    std::cout << m_names.back() << "\n";
}

template<typename T>
void Storage::print(Storable<T> var) const {
    if (m_data.size() <= var.index) {
        throw std::out_of_range("Storable #" + std::to_string(var.index) + " does not exists");
    }
    if (!m_data[var.index].template support<T>()) {
        throw std::runtime_error("Buffer doesn't support type T");
    }
    m_data[var.index].template print<T>(m_names[var.index]);
}

template <typename... Vars>
std::tuple<remove_cvref_t<Vars>...> Storage::add_replace(
    const Storage& s, const std::tuple<Vars...>& vars) {
    soa::assert_storable<Vars...>();

    // Число переменных
    constexpr int n_vars = sizeof...(Vars);

    // Копируем имена
    m_names.resize(n_vars, "");
    [&]<size_t... I>(std::index_sequence<I...>) {
        ((m_names[I] = s.m_names[std::get<I>(vars).index]), ...);
    }(std::make_index_sequence<n_vars>{});

    // Заменяем буферы на новые
    m_data.resize(n_vars, Buffer::dummy());
    [&]<size_t... I>(std::index_sequence<I...>) {
        (m_data[I].replace(s.m_data[I], m_size), ...);
    }(std::make_index_sequence<n_vars>{});

    // Новые Storable имеют последовательные индексы
    std::tuple<Vars...> res = vars;
    [&]<size_t... I>(std::index_sequence<I...>) {
        ((std::get<I>(res).index = I), ...);
    }(std::make_index_sequence<n_vars>{});

    return res;
}

} // namespace zephyr::mesh