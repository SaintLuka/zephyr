#pragma once

#include <zephyr/geom/primitives/bface.h>

#include <zephyr/mesh/euler/amr_storage.h>
#include <zephyr/mesh/euler/eu_face.h>


namespace zephyr::mesh {

using zephyr::geom::Adjacent;
using zephyr::geom::Vector3d;

/// @brief Ячейка (итератор для доступа к элементам класса EuMesh)
/// @details Копирует функции AmrStorage::Item, а также добавляет
/// специфические функции для ячейки: volume, center, faces.
class EuCell {
public:
    // Конструкторы

    /// @brief Конструктор копирования
    EuCell(const EuCell &cell) = default;

    /// @brief Конструктор перемещения
    EuCell(EuCell &&cell) = default;

    /// @brief Основной конструктор. Взять ячейку из хранилища по индексу.
    /// Ячейка может быть из locals или из aliens. Если индекс меньше
    /// locals.size, тогда ячейка берется из хранилища locals, иначе
    /// выбирается ячейка aliens[idx - locals.size]. Если индекс слишком
    /// большой (больше locals.size + aliens.size), тогда возвращается
    /// итератор на конец хранилища locals.
    /// @param idx Индекс ячейки.
    EuCell(AmrStorage &locals, AmrStorage &aliens, int idx);

    /// @brief Вспомогательный конструктор. Взять смежную ячейку из locals
    /// или aliens хранилища по индексу смежности с грани.
    /// @param adj Индекс смежности ячейки.
    EuCell(AmrStorage &locals, AmrStorage &aliens, const Adjacent &adj);


    // Функции итератора

    /// @brief Разыменование итератора (для цикла for)
    inline EuCell &operator*() {
        return *this;
    }

    /// @brief Следующая ячейка
    inline EuCell &operator++() {
        ++m_it;
        return *this;
    }

    /// @brief Ячейка через step
    inline EuCell &operator+=(int step) {
        m_it += step;
        return *this;
    }

    /// @brief Ячейка через step
    inline EuCell operator+(int step) {
        EuCell res(*this);
        res += step;
        return res;
    }

    /// @brief Расстояние между двумя ячейками
    int operator-(const EuCell& cell) const {
        return m_it - cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator<(const EuCell& cell) {
        return m_it < cell.m_it;
    }

    /// @brief Оператор сравнения
    bool operator!=(const EuCell& cell) {
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
    EuFaces faces(Direction dir = Direction::ANY) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(const BFace &face) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(const EuFace &face) const;

    /// @brief Получить смежную ячейку через грань
    EuCell neib(int face_idx) const;

    /// @brief Для данной ячейки вернуть ячейку из local AmrStorage
    /// @param idx Индекс ячейки в локальном хранилище
    EuCell locals(int idx);


    // Далее простые get/set функции

    /// @brief Размерность ячейки
    inline const int& dim() const { return m_it->dim; }

    /// @brief Индекс среди базовых ячеек
    inline const int& b_idx() const { return m_it->b_idx; }

    /// @brief Индекс ячейки на z-кривой
    inline const int& z_idx() const { return m_it->z_idx; }

    /// @brief Индекс новой ячейки (в алгоритмах)
    inline const int& next() const { return m_it->next; }

    /// @brief Уровень адаптации ячейки (0 для базовой)
    inline const int& level() const { return m_it->level; }

    /// @brief Желаемый флаг адаптации
    inline const int& flag() const { return m_it->flag; }

    /// @brief Желаемый флаг адаптации
    inline void set_flag(int flag) { m_it->flag = flag; }

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

    /// @brief Скорпировать вершины в полигон (двумерные ячейки)
    /// Для нелинейных AMR-ячеек возвращает до 8 граней.
    inline geom::PolygonS<8> polygon() const {
        return m_it->polygon();
    }


    // Функции для дебага

    /// @brief Вывести полную информаци о ячейке и её соседях
    void print_neibs_info() const;

    /// @brief Вывести полную информацию о ячейке
    inline void print_info() const { m_it->print_info(); }

    /// @brief Вывести информацию о ячейке в виде python скрипта
    inline void visualize() const { m_it->print_info(); }

protected:
    friend class EuFaces;

    AmrStorage::Iterator m_it;  ///< Итератор по хранилищу
    AmrStorage &m_locals;       ///< Ссылка на хранилище текущего процесса
    AmrStorage &m_aliens;       ///< Ссылка на хранилище соседей
};

} // namespace zephyr::mesh