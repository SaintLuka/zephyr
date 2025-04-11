#pragma once

#include <vector>
#include <cstring>
#include <iostream>

namespace zephyr::mesh {

using Byte = unsigned char;

/// @brief Вспомогательный тип для точечного извлечения данных из Storage
template <typename T>
struct VarExtra {
    int offset;
};

/// @class Storage storage.h
/// @brief Хранилище для расчетных элементов. Каждый элемент хранилища
/// содержит геометрию + данные элемента (расчетные величины).
/// Используется три типа геометрии: AmrCell, MovCell и MovNode.
/// @tparam Geom Геометрический тип данных
template <class Geom>
class Storage {
public:
    /// @brief Пустое и бессмысленное хранилище
    Storage() :
        m_size(0),
        m_itemsize(sizeof(Geom)),
        m_data() {
    }

    /// @brief Хранилище без указания типа для хранения, содержит только геометрию
    explicit Storage(size_t size) :
        m_size(size),
        m_itemsize(sizeof(Geom)),
        m_data(m_size * m_itemsize) {
    }

    /// @brief Хранилище с указанием типа для хранения
    template <class U>
    Storage(size_t size, const U& type) :
        m_size(size),
        m_itemsize(sizeof(Geom) + sizeof(U)),
        m_data(m_size * m_itemsize) {
    }

    /// @brief Хранилище с указанием размера буффера
    /// @param dynamic Фиктивный параметр
    /// @param datasize Размер данных элемента для хранения в байтах
    Storage(size_t size, bool dynamic, size_t datasize) :
        m_size(size),
        m_itemsize(sizeof(Geom) + datasize),
        m_data(m_size * m_itemsize) {
    }

    /// @brief Пустое ли хранилище
    inline bool empty() const {
        return m_size < 1;
    }

    /// @brief Количество элементов хранилища
    inline size_t size() const {
        return m_size;
    }

    /// @brief  Вместимость текущего буфера
    inline size_t capacity() const {
        return m_data.capacity() / m_itemsize;
    }

    /// @brief Размер данных элемента хранилища в байтах
    inline size_t datasize() const {
        return m_itemsize - sizeof(Geom);
    }

    /// @brief Размер элемента хранилища в байтах (> datasize)
    inline size_t itemsize() const {
        return m_itemsize;
    }

    /// @brief Изменить размер хранилища
    inline void resize(size_t new_size) {
        m_size = new_size;
        m_data.resize(m_size * size_t(m_itemsize));
    }

    /// @brief Создать пустое хранилище с таким же размером данных
    inline Storage same() const { return Storage(0, true, datasize()); }

    /// @brief A non-binding request to reduce capacity() to size()
    inline void shrink_to_fit() { m_data.shrink_to_fit(); }

    /// @brief Скопировать элемент хранилища с индексом from в элемент
    /// хранилища с индексом to (геометрия и данные).
    inline void move_item(size_t from, size_t to) {
        std::memcpy(
                m_data.data() + to * m_itemsize,
                m_data.data() + from * m_itemsize,
                m_itemsize);
    }


    /// @brief Представление элемента данных хранилища.
    /// @details Все конструкторы запрещены, тип может быть получен
    /// только в виде ссылки на память из существующего хранилища.
    class Item : public Geom {
    public:
        /// @brief Любое создание запрещено
        Item() = delete;

        /// @brief Копирование запрещено, использовать только
        /// по ссылкам
        Item(const Item&) = delete;

        /// @brief Перемещение запрещено
        Item(Item&&) = delete;

        /// @brief Существующему элементу хранилища можно присводить
        /// геометрию, это позволяет создать хранилище с некоторым числом
        /// элементов, а затем вручную установить геометрию ячеек.
        inline Item& operator=(const Geom& g) {
            Geom::operator=(g);
            return *this;
        }

        /// @brief Указатель на начало элемента
        inline Byte* ptr() {
            return reinterpret_cast<Byte *>(this);
        }

        /// @brief Указатель на начало элемента
        inline const Byte* ptr() const {
            return reinterpret_cast<const Byte *>(this);
        }

        /// @brief Указатель на данные элемента
        /// (указатель на память сразу за экземпляром класса)
        inline Byte* data() {
            return reinterpret_cast<Byte *>(this) + sizeof(Geom);
        }

        /// @brief Указатель на данные элемента
        /// (указатель на память сразу за экземпляром класса)
        inline const Byte* data() const {
            return reinterpret_cast<const Byte *>(this) + sizeof(Geom);
        }

        /// @brief Разыменовать данные элемента как тип U.
        template<class U>
        inline const U& data(const U &) const {
            return *reinterpret_cast<const U *>(data());
        }

        /// @brief Разыменовать данные элемента как тип U.
        template<class U>
        inline U& data(const U &) {
            return *reinterpret_cast<U *>(data());
        }

        /// @brief Разыменовать данные элемента как тип U.
        template<class U>
        inline const U& operator()(const U &) const {
            return *reinterpret_cast<const U *>(data());
        }

        /// @brief Разыменовать данные элемента как тип U.
        template<class U>
        inline U& operator()(const U &) {
            return *reinterpret_cast<U *>(data());
        }

        /// @brief Извлечь данные элемента
        template<class T>
        inline const T& operator()(const VarExtra<T>& var) const {
            return *reinterpret_cast<const T *>(data() + var.offset);
        }

        /// @brief Извлечь данные элемента
        template<class T>
        inline T& operator()(const VarExtra<T>& var) {
            return *reinterpret_cast<T *>(data() + var.offset);
        }
    };


    /// @brief Итератор по элементам хранилища
    class Iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Item&;
        using pointer    = Item*;
        using reference  = Item&;

        /// @brief Создание нулевого итератора из nullptr
        inline Iterator(std::nullptr_t null)
                : m_ptr(nullptr), m_itemsize(0) { }

        /// @brief Создание нормального итератора
        inline Iterator(Byte *ptr, size_t itemsize)
                : m_ptr(ptr), m_itemsize(itemsize) { }

        /// @brief Проверка на совпадение с нулевым указателем
        inline operator bool() const {
            return m_ptr;
        }

        /// @brief Размер данных в байтах
        inline size_t datasize() const {
            return m_itemsize - int(sizeof(Geom));
        }

        /// @brief -> Перенаправляется на Storage<Geom>::Item*
        inline Item *operator->() {
            return reinterpret_cast<Item *>(m_ptr);
        }

        /// @brief -> Перенаправляется на Storage<Geom>::Item*
        inline const Item *operator->() const {
            return reinterpret_cast<const Item *>(m_ptr);
        }

        /// @brief Разыменование в Storage<Geom>::Item&
        inline Item &operator*() {
            return *reinterpret_cast<Item *>(m_ptr);
        }

        /// @brief Разыменование в Storage<Geom>::Item&
        inline const Item &operator*() const {
            return *reinterpret_cast<const Item *>(m_ptr);
        }

        /// @brief Префиксный инкремент (++it)
        inline Iterator &operator++() {
            m_ptr += m_itemsize;
            return *this;
        }

        /// @brief Префиксный декремент (++it)
        inline Iterator &operator--() {
            m_ptr -= m_itemsize;
            return *this;
        }

        /// @brief Сдвиг итератора (it += step)
        inline Iterator &operator+=(size_t step) {
            m_ptr += step * m_itemsize;
            return *this;
        }

        /// @brief Смещенный итератор (it2 = it + step)
        inline Iterator operator+(size_t step) {
            return {m_ptr + step * m_itemsize, m_itemsize};
        }

        /// @brief Расстояние между парой итераторов
        inline size_t operator-(const Iterator &it) const {
            return (m_ptr - it.m_ptr) / m_itemsize;
        }

        /// @brief Разыменование со сдвигом
        inline Item& operator[](size_t idx) const {
            return *reinterpret_cast<Item *>(m_ptr + idx * m_itemsize);
        }

        // Решил не мелочиться и реализовать все операции сравнения

        inline bool operator<(const Iterator &it) const { return m_ptr < it.m_ptr; }

        inline bool operator>(const Iterator &it) const { return m_ptr > it.m_ptr; }

        inline bool operator<=(const Iterator &it) const { return m_ptr <= it.m_ptr; }

        inline bool operator>=(const Iterator &it) const { return m_ptr >= it.m_ptr; }

        inline bool operator==(const Iterator &it) const { return m_ptr == it.m_ptr; }

        inline bool operator!=(const Iterator &it) const { return m_ptr != it.m_ptr; }

    private:
        Byte *m_ptr;        ///< Ссылка на начало данных
        size_t m_itemsize;  ///< Размер элемента в байтах
    };

    /// @brief Получить указатель на данные по индексу
    inline const Byte* data_at(size_t idx) const {
        return m_data.data() + m_itemsize * idx + sizeof(Geom);
    }

    /// @brief Получить итератор по индексу
    inline Iterator iterator(size_t idx) {
        return {m_data.data() + m_itemsize * idx, m_itemsize};
    }

    /// @brief Получить элемент по индексу
    inline Item& item(size_t idx) {
        return *reinterpret_cast<Item*>(m_data.data() + m_itemsize * idx);
    }

    /// @brief Получить элемент по индексу
    inline const Item& item(size_t idx) const {
        return *reinterpret_cast<const Item*>(m_data.data() + m_itemsize * idx);
    }

    /// @brief Ссылка на элемент хранилища
    inline Item& operator[](size_t idx) {
        return item(idx);
    }

    /// @brief Константная ссылка на элемент хранилища
    inline const Item& operator[](size_t idx) const {
        return item(idx);
    }

    /// @brief Итератор на начало хранилища
    Storage::Iterator begin() {
        return iterator(0);
    }

    /// @brief Итератор за концом хранилища
    Storage::Iterator end() {
        return iterator(m_size);
    }

private:
    size_t m_size;
    size_t m_itemsize;
    std::vector<Byte> m_data;
};

} // namespace zephyr::mesh