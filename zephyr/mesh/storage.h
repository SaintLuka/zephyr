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


    /// @brief Элемент хранилища
    struct Item {
    public:

        Item(Byte *ptr, int itemsize) : m_ptr(ptr), m_itemsize(itemsize) {}

        virtual Item &operator*() { return *this; }

        virtual Item &operator++() { m_ptr += m_itemsize; return *this; }

        bool operator<(const Item &it) const { return m_ptr < it.m_ptr; }

        bool operator>(const Item &it) const { return m_ptr > it.m_ptr; }

        bool operator<=(const Item &it) const { return m_ptr <= it.m_ptr; }

        bool operator>=(const Item &it) const { return m_ptr >= it.m_ptr; }

        bool operator==(const Item &it) const { return m_ptr == it.m_ptr; }

        bool operator!=(const Item &it) const { return m_ptr != it.m_ptr; }

        int operator-(const Item &it) const { return int(it.m_ptr - m_ptr) / m_itemsize; }

        /// @brief Ссылка на геометрию
        inline Byte* geom_ptr() { return m_ptr; }

        /// @brief Ссылка на геометрию
        inline Cell &geom() { return *((Cell *) m_ptr); }

        /// @brief Ссылка на геометрию
        inline const Cell &geom() const { return *((Cell *) m_ptr); }

        /// @brief Приведение к чистой геометрии
        inline operator Cell&() { return *((Cell *) m_ptr); }

        /// @brief Приведение к чистой геометрии
        inline operator const Cell &() const { return *((Cell *) m_ptr); }

        /// @brief Вывести полную информацию о ячейке
        inline void print_info() const { geom().print_info(); }

        /// @brief Вывести информацию о ячейке в виде python скрипта
        inline void visualize() const { geom().print_info(); }

        /// @brief Размерность ячейки
        inline const int& dim() const { return geom().dim; }

        /// @brief Индекс среди базовых ячеек
        inline const int& b_idx() const { return geom().b_idx; }

        /// @brief Индекс ячейки на z-кривой
        inline const int& z_idx() const { return geom().z_idx; }

        /// @brief Индекс новой ячейки (в алгоритмах)
        inline const int& next() const { return geom().next; }

        /// @brief Уровень адаптации ячейки (0 для базовой)
        inline const int& level() const { return geom().level; }

        /// @brief Желаемый флаг адаптации
        inline const int& flag() const { return geom().flag; }

        /// @brief Желаемый флаг адаптации
        inline void set_flag(int flag) { geom().flag = flag; }

        /// @brief Ранг процесса, который обрабатывает ячейку
        inline const int& rank() const { return geom().rank; }

        /// @brief Пометить ячейку как неопределенную (вне сетки)
        inline void set_undefined() { geom().set_undefined(); }

        /// @brief Является ли ячейка актуальной?
        inline bool is_actual() const { return geom().is_actual(); }

        /// @brief Является ли ячейка неопределенной?
        inline bool is_undefined() const { return geom().is_undefined(); }

        /// @brief Барицентр ячейки
        inline const Vector3d& center() const { return geom().coords; }

        /// @brief Линейный размер ячейки
        inline const double& size() const { return geom().size; }

        /// @brief Площадь (в 2D) или объем (в 3D) ячейки
        inline double volume() const { return geom().volume(); }

        /// @brief Список вершин
        inline const geom::Vertices& vertices() const { return geom().vertices; };

        /// @brief Конкретная вершина
        inline geom::Vertex& vertices(int i) { return geom().vertices[i]; };

        /// @brief Конкретная вершина
        inline const geom::Vertex& vertices(int i) const { return geom().vertices[i]; };

        /// @brief Список граней
        inline const geom::Faces& faces() const { return geom().faces; };

        /// @brief Конкретная грань
        inline geom::Face& faces(int i) { return geom().faces[i]; };

        /// @brief Конкретная грань
        inline const geom::Face& faces(int i) const { return geom().faces[i]; };


        /// @brief Размер данных в байтах
        inline int datasize() const { return m_itemsize - int(sizeof(Cell)); }

        /// @brief Ссылка на начало данных
        inline Byte* data_ptr() const { return m_ptr + sizeof(Cell); }

        /// @brief Копирует данные из одной ячейки в другую.
        /// @param src Указатель на источник данных.
        /// @param dest Указатель на область, куда надо записать данные.
        inline void copy_to(const Storage::Item& dest) const {
            std::memcpy(m_ptr, dest.m_ptr, m_itemsize);
        }

        /// @brief Ссылка на данные
        template<class T>
        inline T &data() { return *((T *) (m_ptr + sizeof(Cell))); }

        /// @brief Ссылка на данные
        template<class T>
        inline T &operator()(const T&) { return *((T *) (m_ptr + sizeof(Cell))); }

    protected:
        Byte *m_ptr;      ///< Ссылка на начало данных
        int m_itemsize;   ///< Размер элемента в байтах
    };

    Item operator[](int i);

    Item begin();

    Item end();

private:

    /// @brief Проинициализировать хранилище
    void init(int datasize, int size);

    int m_size;
    int m_itemsize;
    std::vector<Byte> m_data;
};

} // mesh
} // zephyr