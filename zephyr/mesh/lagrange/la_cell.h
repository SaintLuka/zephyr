#pragma once

#include <zephyr/geom/primitives/bface.h>

#include <zephyr/mesh/lagrange/mov_storage.h>
#include <zephyr/mesh/lagrange/la_face.h>


namespace zephyr::mesh {

using zephyr::geom::Adjacent;
using zephyr::geom::Vector3d;

/// @brief Ячейка (итератор для доступа к элементам класса LaMesh)
/// @details Копирует функции CellStorage::Item, а также добавляет
/// специфические функции для ячейки: volume, center, faces.
class LaCell {
public:
    // Конструкторы

    /// @brief Конструктор копирования
    LaCell(const LaCell &cell) = default;

    /// @brief Конструктор перемещения
    LaCell(LaCell &&cell) = default;

    /// @brief Основной конструктор. Взять ячейку из хранилища по индексу.
    /// Ячейка может быть из locals или из aliens. Если индекс меньше
    /// locals.size, тогда ячейка берется из хранилища locals, иначе
    /// выбирается ячейка aliens[idx - locals.size]. Если индекс слишком
    /// большой (больше locals.size + aliens.size), тогда возвращается
    /// итератор на конец хранилища locals.
    /// @param idx Индекс ячейки.
    LaCell(CellStorage &locals, int idx);

    /// @brief Вспомогательный конструктор. Взять смежную ячейку из locals
    /// или aliens хранилища по индексу смежности с грани.
    /// @param adj Индекс смежности ячейки.
    LaCell(CellStorage &locals, const Adjacent &adj);


    // Функции итератора

    /// @brief Разыменование итератора (для цикла for)
    inline LaCell &operator*() {
        return *this;
    }

    /// @brief Следующая ячейка
    inline LaCell &operator++() {
        ++m_it;
        return *this;
    }

    /// @brief Ячейка через step
    inline LaCell &operator+=(int step) {
        m_it += step;
        return *this;
    }

    /// @brief Ячейка через step
    inline LaCell operator+(int step) {
        LaCell res(*this);
        res += step;
        return res;
    }

    /// @brief Расстояние между двумя ячейками
    int operator-(const LaCell& cell) const {
        return m_it - cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator<(const LaCell& cell) {
        return m_it < cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator!=(const LaCell& cell) {
        return m_it != cell.m_it;
    }


    // Получение данных, обход граней и соседей

    /// @brief Указатель на данные ячейки
    inline Byte* data() {
        return m_it->data();
    }

    /// @brief Указатель на данные ячейки
    inline const Byte* data() const {
        return m_it->data();
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    U& operator()(const U &) {
        return *reinterpret_cast<U *>(m_it->data());
    }

    /// @brief Ссылка на данные ячейки
    template<class U>
    const U& operator()(const U &) const {
        return *reinterpret_cast<const U *>(m_it->data());
    }

    /// @brief Итератор по граням
    /// @param dir Выбрать грани по некоторым направлениям
    LaFaces faces(Direction dir = Direction::ANY) const;

    /// @brief Получить смежную ячейку через грань
    LaCell neib(const BFace &face) const;

    /// @brief Получить смежную ячейку через грань
    LaCell neib(const LaFace &face) const;

    /// @brief Получить смежную ячейку через грань
    LaCell neib(int face_idx) const;

    /// @brief Для данной ячейки вернуть ячейку из local CellStorage
    /// @param idx Индекс ячейки в локальном хранилище
    LaCell locals(int idx);


    // Далее простые get/set функции
    /*
    /// @brief Размерность ячейки
    inline const int& dim() const { return m_it->dim; }

    /// @brief Индекс новой ячейки (в алгоритмах)
    inline const int& next() const { return m_it->next; }

    /// @brief Барицентр ячейки
    inline const Vector3d& center() const { return m_it->center; }

    /// @brief Линейный размер ячейки
    inline const double& size() const { return m_it->size; }

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    inline double volume() const { return m_it->volume(); }

    /// @brief Вершина ячейки по индексу
    inline const Vector3d& vs(int idx) const {
        return m_it->vertices[idx];
    }
     */

protected:
    friend class LaFaces;

    CellStorage::Iterator m_it;  ///< Итератор по хранилищу
    CellStorage &m_locals;       ///< Ссылка на хранилище текущего процесса
};

} // namespace zephyr::mesh