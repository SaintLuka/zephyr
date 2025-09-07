#pragma once

#include <string>
#include <zephyr/geom/vector.h>
#include <zephyr/utils/span.h>

#ifdef ZEPHYR_MPI
#include <zephyr/utils/mpi.h>
#endif

namespace zephyr::utils {

/// @brief Достаточно простой тип, который допускается хранить в буфере
template<typename T>
struct is_simple_type {
    static constexpr bool value =
        !std::is_pointer_v<T> &&
        !std::is_array_v<T> &&
        !std::is_polymorphic_v<T> &&
         std::is_trivially_destructible_v<T>;
};

template<typename T>
static constexpr bool is_simple_type_v = is_simple_type<T>::value;

/// @brief Буфер в виде массива байтов для хранения элементов произвольного
/// размера. Элементы буфера могут быть *скалярными* или *векторными*.
/// Под векторными понимаются типы, для которых std::is_array_v == true,
/// то есть массивы T[] и T[N]. Таким образом, в буфере можно хранить
/// тройки int'ов, к примеру. При этом тип std::array<int, 3> при хранении
/// рассматривается как скалярный.
class Buffer final {
public:
    /// @{ @name Конструкторы / фабричные функции

    /// @brief Создать буфер для скалярного типа
    /// @size Число элементов
    /// @code
    ///   Buffer buf = Buffer::make<int>("val", 1000);
    /// @endcode
    template<class T>
    static std::enable_if_t<!std::is_array_v<T>, Buffer>
    make(size_t size = 0) {
        static_assert(is_simple_type_v<T>, "Type must be a simple type");
        return Buffer{sizeof(T), size, 1};
    }

    /// @brief Создать буфер для статического векторного типа
    /// @param size Число элементов
    /// @code
    ///   Buffer arr = Buffer::make<int[3]>("vec", 1000);
    /// @endcode
    template<class T>
    static std::enable_if_t<(std::is_array_v<T> && std::extent_v<T> > 0), Buffer>
    make(size_t size = 0) {
        static_assert(is_simple_type_v<std::remove_extent_t<T>>, "Type must be a simple type");
        return Buffer{sizeof(T), size, std::extent_v<T>};
    }

    /// @brief Создать буфер для динамического векторного типа
    /// @param count Число компонент элемента
    /// @param size Число элементов
    /// @code
    ///   Buffer buf = Buffer::make<int[]>("vec", 3, 1000);
    /// @endcode
    template<class T>
    static std::enable_if_t<(std::is_array_v<T> && std::extent_v<T> == 0), Buffer>
    make(size_t count, size_t size = 0) {
        static_assert(is_simple_type_v<std::remove_extent_t<T>>, "Type must be a simple type");
        return Buffer{count * sizeof(std::remove_extent_t<T>), size, count};
    }

    /// @brief Пустой буфер с типом int
    static Buffer dummy() { return Buffer::make<int>(); }

    /// @brief Создать аналогичный буфер
    Buffer same(size_t size = 0) const { return Buffer{m_element_size, size, m_count}; }

    /// @}

    /// @{ @name Общая информация

    /// @brief Число компонент векторного типа
    size_t count() const noexcept { return m_count; }

    /// @brief Скалярные элементы?
    bool is_scalar() const noexcept { return m_count == 1; }

    /// @brief Векторные элементы
    bool is_vector() const noexcept { return m_count > 1; }

    /// @brief Буфер поддерживает скалярный тип?
    template <typename T>
    std::enable_if_t<!std::is_array_v<T>, bool>
    support() const { return m_count == 1 && m_element_size == sizeof(T); }

    /// @brief Буфер поддерживает статический векторный тип?
    template <typename T>
    std::enable_if_t<(std::is_array_v<T> && std::extent_v<T> > 0), bool>
    support() const { return m_element_size == sizeof(std::remove_extent_t<T>) * m_count && std::extent_v<T> == m_count; }

    /// @brief Буфер поддерживает динамический векторный тип?
    template <typename T>
    std::enable_if_t<(std::is_array_v<T> && std::extent_v<T> == 0), bool>
    support() const { return m_element_size == sizeof(std::remove_extent_t<T>) * m_count; }

#ifdef ZEPHYR_MPI
    /// @brief MPI-тип для пересылок, по умолчанию contiguous из байтов
    MPI_Datatype dtype() const { return m_dtype; }
#endif

    /// @}

    /// @{ @name Стандартные функции динамических массивов

    /// @brief Пустой буфер?
    bool empty() const noexcept { return m_data.empty(); }

    /// @brief Число элементов буфера
    size_t size() const noexcept { return m_data.size() / m_element_size; }

    /// @brief Предельное число элементов без переаллокации
    size_t capacity() const noexcept { return m_data.capacity() / m_element_size; }

    /// @brief Очистить буфер
    void clear() noexcept { m_data.clear(); }

    /// @brief Изменить число элементов буфера
    void resize(size_t new_size) { m_data.resize(new_size * m_element_size); }

    /// @brief Расширить допустимые размеры
    void reserve(size_t new_capacity) { m_data.reserve(new_capacity * m_element_size); }

    /// @brief Уменьшить размер буфера до актуальных данных
    void shrink_to_fit() { m_data.shrink_to_fit(); }

    /// @brief Размер элемента буфера в байтах
    size_t element_size() const noexcept { return m_element_size; }

    /// @brief Размер всего буфера в байтах
    size_t byte_size() const noexcept { return m_data.size(); }

    /// @}

    /// @{ @name Работа с данными

    /// @brief Указатель на начало
    std::byte* data() { return m_data.data(); }

    /// @brief Константный указатель на начало
    const std::byte* data() const { return m_data.data(); }

    /// @brief Указатель на начало данных по индексу
    std::byte* get_ptr(size_t index) { return m_data.data() + index * m_element_size; }

    /// @brief Константный указатель на начало данных по индексу
    const std::byte* get_ptr(size_t index) const { return m_data.data() + index * m_element_size; }

    /// @brief Доступ к элементу скалярного типа
    template<typename T>
    std::enable_if_t<!std::is_array_v<T>, T&>
    get_val(size_t index) { return *reinterpret_cast<T*>(get_ptr(index)); }

    /// @brief Доступ к элементу скалярного типа
    template<typename T>
    std::enable_if_t<!std::is_array_v<T>, const T&>
    get_val(size_t index) const { return *reinterpret_cast<const T*>(get_ptr(index)); }

    /// @brief Доступ к элементу векторного типа
    template<typename T>
    std::enable_if_t<std::is_array_v<T>, span<std::remove_extent_t<T>>>
    get_val(size_t index) {
        using V = std::remove_extent_t<T>;
        return span<V>{reinterpret_cast<V*>(get_ptr(index)), m_count};
    }

    /// @brief Доступ к элементу векторного типа
    template<typename T>
    std::enable_if_t<std::is_array_v<T>, span<const std::remove_extent_t<T>>>
    get_val(size_t index) const {
        using V = std::remove_extent_t<T>;
        return span<const V>{reinterpret_cast<const V*>(get_ptr(index)), m_count}; }

    /// @brief Скопировать данные с индекса from на индекс to
    void copy_data(size_t from, size_t to) {
        std::memcpy(get_ptr(to), get_ptr(from), m_element_size);
    }

    /// @brief Скопировать все данные с индекса from текущего буфера
    /// в буфер dst по индексу to.
    void copy_data(size_t from, Buffer& dst, size_t to) const {
        std::memcpy(dst.get_ptr(to), get_ptr(from), m_element_size);
    }

    /// @brief Заменить текущий буфер на буфер с другим element_size и count
    void replace(const Buffer& other, size_t size) {
        m_element_size = other.m_element_size;
        m_count = other.m_count;
#ifdef ZEPHYR_MPI
        m_dtype = other.m_dtype;
#endif
        resize(size);
    }

    /// @}

    /// @brief Вывести массив данных, для типа должен быть определен operator<<
    template<typename T>
    void print(const std::string& name = "") const;

protected:
    /// @brief Буфер с заданным размером элементов, заданной длиной и кратностью
    explicit Buffer(size_t element_size, size_t size, size_t count)
        : m_element_size(element_size),
          m_data(size * m_element_size),
          m_count(count) {
#ifdef ZEPHYR_MPI
        m_dtype = mpi::datatype::contiguous(element_size, MPI_BYTE);
#endif
    }

    size_t m_element_size;          ///< Размер элемента в байтах
    std::vector<std::byte> m_data;  ///< Байтовый буфер
    size_t m_count = 1;             ///< Число компонент (для векторного типа)

#ifdef ZEPHYR_MPI
    mpi::datatype m_dtype;          ///< MPI_Type для обменов
#endif
};


template<typename T>
void Buffer::print(const std::string& name) const {
    std::cout << "Buffer ";
    if (!name.empty()) {
        std::cout << "\"" << name << "\" ";
    }
    std::cout << "(size: " << size();
    if (is_vector()) { std::cout << ", count: " << m_count; }
    std::cout << "):";
    if (empty()) {
        std::cout << "\n  [ ]\n";
        return;
    }

    using V = std::conditional_t<std::is_array_v<T>, std::remove_extent_t<T>, T>;
    std::cout << "\n";
    if (is_scalar()) {
        std::cout << "  [ ";
        if constexpr (!std::is_same_v<T, geom::Vector3d>) {
            for (size_t i = 0; i < size() - 1; ++i) { std::cout << get_val<T>(i) << ", "; }
            std::cout << get_val<V>(size() - 1) << " ]\n";
        }
        else {
            for (size_t i = 0; i < size() - 1; ++i) { std::cout << get_val<T>(i).transpose() << ", "; }
            std::cout << get_val<V>(size() - 1).transpose() << " ]\n";
        }
    }
    else {
        std::cout << "  [ ";
        for (size_t i = 0; i < size() - 1; ++i) {
            std::cout << get_val<V[]>(i) << ", ";
        }
        std::cout << get_val<V[]>(size() - 1) << " ]\n";
    }
}

} // namespace zephyr::utils