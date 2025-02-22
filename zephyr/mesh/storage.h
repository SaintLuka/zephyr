#pragma once

#include <vector>
#include <cstring>
#include <iostream>

namespace zephyr::mesh {

using Byte = unsigned char;

/// @brief Вспомогательный тип для
/// точечного извлечения данных из Storage
template <typename T>
struct VarExtra {
    int offset;
};

/// @class Хранилище для расчетных элементов. Каждый элемент хранилища
/// содержит геометрию + данные элемента (расчетные величины).
/// Используется три типа геометрии: AmrCell, BCell и BNode.
/// @tparam Geom Геометрический тип данных
template <class Geom>
class Storage {
public:

    /// @brief Пустое и бессмысленное хранилище
    Storage() {
        m_size = 0;
        m_itemsize = int(sizeof(Geom));
        m_data = {};
    };

    /// @brief Хранилище без указания типа для хранения,
    /// содержит только геометрию
    Storage(int size) {
        m_size = std::max(0, size);
        m_itemsize = int(sizeof(Geom));
        m_data.resize(m_size * m_itemsize);
    };

    /// @brief Создать хранилище размера size для хранения типа U
    /// @tparam U Тип данных, который предполагается располагать в хранилище,
    /// используется только для выделения достаточного количества памяти.
    /// @warning Не устанавливает значения в хранилище
    template<class U>
    explicit Storage(const U &, int size) {
        m_size = std::max(0, size);
        m_itemsize = int(sizeof(Geom) + sizeof(U));
        m_data.resize(m_size * m_itemsize);
    }

    /// @brief Пустое ли хранилище
    inline bool empty() const {
        return m_size < 1;
    }

    /// @brief Количество элементов хранилища
    inline int size() const {
        return m_size;
    }

    /// @brief Размер данных элемента хранилища в байтах
    inline int datasize() const {
        return m_itemsize - int(sizeof(Geom));
    }

    /// @brief Размер элемента хранилища в байтах (> datasize)
    inline int itemsize() const {
        return m_itemsize;
    }

    /// @brief Изменить размер хранилища
    inline void resize(int new_size) {
        m_size = std::max(0, new_size);
        m_data.resize(m_size * m_itemsize);
    }

    /// @brief Скопировать элемент хранилища с индексом from в элемент
    /// хранилища с индексом to (геометрия и данные).
    inline void move_item(int from, int to) {
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
        inline Byte * data() {
            return reinterpret_cast<Byte *>(this) + sizeof(Geom);
        }

        /// @brief Указатель на данные элемента
        /// (указатель на память сразу за экземпляром класса)
        inline const Byte * data() const {
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
        inline const T& operator()(const VarExtra<T> & var) const {
            return *reinterpret_cast<const T *>(data() + var.offset);
        }

        /// @brief Извлечь данные элемента
        template<class T>
        inline T& operator()(const VarExtra<T> & var) {
            return *reinterpret_cast<T *>(data() + var.offset);
        }
    };


    /// @brief Итератор по элементам хранилища
    class Iterator {
    public:
        /// @brief Создание нулевого итератора из nullptr
        inline Iterator(std::nullptr_t null)
                : m_ptr(nullptr), m_itemsize(0) { }

        /// @brief Создание нормального итератора
        inline Iterator(Byte *ptr, int itemsize)
                : m_ptr(ptr), m_itemsize(itemsize) { }

        /// @brief Проверка на совпадение с нулевым указателем
        inline operator bool() const {
            return m_ptr;
        }

        /// @brief Размер данных в байтах
        inline int datasize() const {
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

        /// @brief Сдвиг итератора (it += step)
        inline Iterator &operator+=(int step) {
            m_ptr += step * m_itemsize;
            return *this;
        }

        /// @brief Смещенный итератор (it2 = it + step)
        inline Iterator operator+(int step) {
            return {m_ptr + step * m_itemsize, m_itemsize};
        }

        /// @brief Расстояние между парой итераторов
        inline int operator-(const Iterator &it) const {
            return int(m_ptr - it.m_ptr) / m_itemsize;
        }

        // Решил не мелочиться и реализовать все операции сравнения

        inline bool operator<(const Iterator &it) const { return m_ptr < it.m_ptr; }

        inline bool operator>(const Iterator &it) const { return m_ptr > it.m_ptr; }

        inline bool operator<=(const Iterator &it) const { return m_ptr <= it.m_ptr; }

        inline bool operator>=(const Iterator &it) const { return m_ptr >= it.m_ptr; }

        inline bool operator==(const Iterator &it) const { return m_ptr == it.m_ptr; }

        inline bool operator!=(const Iterator &it) const { return m_ptr != it.m_ptr; }

    private:
        Byte *m_ptr;     ///< Ссылка на начало данных
        int m_itemsize;  ///< Размер элемента в байтах
    };

    /// @brief Получить указатель на данные по индексу
    inline const Byte* data_at(int idx) const {
        return m_data.data() + m_itemsize * idx + sizeof(Geom);
    }

    /// @brief Получить итератор по индексу
    inline Iterator iterator(int idx) {
        return {m_data.data() + m_itemsize * idx, m_itemsize};
    }

    /// @brief Получить элемент по индексу
    inline Item& item(int idx) {
        return *reinterpret_cast<Item*>(m_data.data() + m_itemsize * idx);
    }

    /// @brief Получить элемент по индексу
    inline const Item& item(int idx) const {
        return *reinterpret_cast<const Item*>(m_data.data() + m_itemsize * idx);
    }

    /// @brief Ссылка на элемент хранилища
    inline Item& operator[](int idx) {
        return item(idx);
    }

    /// @brief Константная ссылка на элемент хранилища
    inline const Item& operator[](int idx) const {
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
    int m_size;
    int m_itemsize;
    std::vector<Byte> m_data;
};

} // namespace zephyr::mesh