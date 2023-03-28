#pragma once

#include <iostream>
#include <array>
#include <vector>

#include <zephyr/geom/cell.h>

namespace zephyr { namespace mesh {

using zephyr::geom::Cell;
using zephyr::geom::Vector3d;

using Byte = unsigned char;

class Storage {
public:

    // Запрещаем создание без типа
    Storage() = delete;

    // Запрещаем любое копирование

    Storage(Storage &) = delete;

    Storage(Storage &&) = delete;

    Storage &operator=(Storage &) = delete;

    Storage &operator=(Storage &&) = delete;

    /// @brief Создать хранилище размера size для хранения типа T
    /// @warning Не устанавливает значения в хранилище
    template<class T>
    explicit Storage(const T &, int size = 0) { init(sizeof(T), size); }



    /// @brief Пустое ли хранилище
    bool empty() const;

    /// @brief Количество элементов хранилища
    int size() const;

    /// @brief Размер данных элемента хранилища в байтах
    int datasize() const;

    /// @brief Размер всего элемента хранилища в байтах
    int itemsize() const;

    /// @brief Изменить размер хранилища
    void resize(int new_size);


    /// @brief Итератор по хранилищу
    struct iterator {
    public:

        iterator(Byte *ptr, int itemsize) : m_ptr(ptr), m_itemsize(itemsize) {}

        /// @brief Присвоение таких странных операторов имеет неочевидное поведение
        iterator &operator=(const iterator &) = delete;

        virtual const iterator &operator*() { return *this; }

        virtual const iterator &operator++() { m_ptr += m_itemsize; return *this; }

        bool operator<(const iterator &it) const { return m_ptr < it.m_ptr; }

        bool operator>(const iterator &it) const { return m_ptr > it.m_ptr; }

        bool operator<=(const iterator &it) const { return m_ptr <= it.m_ptr; }

        bool operator>=(const iterator &it) const { return m_ptr >= it.m_ptr; }

        bool operator==(const iterator &it) const { return m_ptr == it.m_ptr; }

        bool operator!=(const iterator &it) const { return m_ptr != it.m_ptr; }

        /// @brief Ссылка на геометрию
        inline Cell &geom() { return *((Cell *) m_ptr); }

        /// @brief Ссылка на геометрию
        inline const Cell &geom() const { return *((const Cell *) m_ptr); }

        /// @brief Размерность ячейки
        inline int dim() const { return geom().dim; }

        /// @brief Индекс среди базовых ячеек
        inline int base_id() const { return geom().amr.base_id; }

        /// @brief Индекс ячейки на z-кривой
        inline int z() const { return geom().amr.z; }

        /// @brief Индекс новой ячейки (в алгоритмах)
        inline int next() const { return geom().amr.next; }

        /// @brief Уровень адаптации ячейки (0 для базовой)
        inline short level() const { return geom().amr.level; }

        /// @brief Желаемый флаг адаптации
        inline short flag() const { return geom().amr.flag; }

        /// @brief Барицентр ячейки
        inline const Vector3d& coords() const { return geom().coords; }

        /// @brief Барицентр ячейки
        inline const Vector3d& center() const { return geom().coords; }

        /// @brief Барицентр ячейки
        inline const Vector3d& centroid() const { return geom().coords; }

        /// @brief Линейный размер ячейки
        inline double size() const { return geom().size; }

        /// @brief Площадь (в 2D) или объем (в 3D) ячейки
        inline double volume() const { return geom().volume(); }

        /// @brief Ссылка на данные
        template<class T>
        inline T &data() {
            return *((T *) (m_ptr + sizeof(Cell)));
        }

        /// @brief Ссылка на данные
        template<class T>
        inline T &operator()(const T&) {
            return *((T *) (m_ptr + sizeof(Cell)));
        }

    protected:
        Byte *m_ptr;      ///< Ссылка на начало данных
        int m_itemsize;   ///< Размер элемента в байтах
    };

    iterator operator[](int i);

    iterator begin();

    iterator end();



private:

    /// @brief Проинициализировать хранилище
    void init(int datasize, int size);

    int m_size;
    int m_itemsize;
    std::vector<Byte> m_data;
};

} // mesh
} // zephyr