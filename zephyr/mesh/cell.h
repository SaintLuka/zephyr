#pragma once

#include <zephyr/geom/primitives/basic_face.h>

#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/face.h>


namespace zephyr::mesh {

using zephyr::geom::Adjacent;
using zephyr::geom::Vector3d;

/// @brief Ячейка (итератор для доступа к элементам класса Mesh)
/// @details Копирует функции AmrStorage::Iterator, а также добавляет
/// специфические функции для ячейки: volume, center, faces.
class ICell {
public:
    /// @brief Конструктор копирования
    ICell(const ICell &cell) = default;

    /// @brief Конструктор перемещения
    ICell(ICell &&cell) = default;

    /// @brief Основной конструктор. Взять ячейку из хранилища по индексу.
    /// Ячейка может быть из locals или из aliens. Если индекс меньше
    /// locals.size, тогда ячейка берется из хранилища locals, иначе
    /// выбирается ячейка aliens[idx - locals.size]. Если индекс слишком
    /// большой (больше locals.size + aliens.size), тогда возвращается
    /// итератор на конец хранилища locals.
    /// @param idx Индекс ячейки.
    ICell(AmrStorage &locals, AmrStorage &aliens, int idx);

    /// @brief Вспомогательный конструктор. Взять смежную ячейку из locals
    /// или aliens хранилища по индексу смежности с грани.
    /// @param adj Индекс смежности ячейки.
    ICell(AmrStorage &locals, AmrStorage &aliens, const Adjacent &adj);

    /// @brief Для данной ячейки вернуть ячейку из local AmrStorage
    /// @param idx Индекс ячейки в локальном хранилище
    ICell locals(int idx);

    /// @brief Разыменование итератора (для цикла for)
    ICell &operator*() {
        return *this;
    }

    /// @brief Следующая ячейка
    ICell &operator++() {
        ++m_it;
        return *this;
    }

    /// @brief Следующая ячейка
    ICell &operator+=(int step) {
        m_it += step;
        return *this;
    }

    /// @brief Следующая ячейка
    ICell operator+(int step) {
        ICell res(*this);
        res += step;
        return res;
    }

    int operator-(const ICell& cell) const {
        return cell.m_it - m_it;
    }

    bool operator<(const ICell& cell) {
        return m_it < cell.m_it;
    }

    bool operator!=(const ICell& cell) {
        return m_it != cell.m_it;
    }

    template<class U>
    const U& operator()(const U &) const {
        return *reinterpret_cast<const U *>(m_it->data());
    }

    template<class U>
    U& operator()(const U &) {
        return *reinterpret_cast<U *>(m_it->data());
    }

    geom::AmrCell& geom() {
        return m_it->geom();
    }

    const geom::AmrCell& geom() const {
        return m_it->geom();
    }

    /// @brief Массив граней
    IFaces faces(Direction dir = Direction::ANY) const;

    /// @brief Получить смежную ячейку через грань
    ICell neib(const AmrFace &face) const;

    /// @brief Получить смежную ячейку через грань
    ICell neighbor(const AmrFace &face) const;

    /// @brief Получить указатель на смежную ячейку через грань
    ICell neib(const IFace &face) const;

    /// @brief Получить смежную ячейку через грань
    ICell neighbor(const IFace &face) const;

    /// @brief Получить указатель на смежную ячейку через грань
    ICell neib(int face_idx) const;

    /// @brief Получить смежную ячейку через грань
    ICell neighbor(int face_idx) const;

    /// @brief Вывести полную информаци о ячейке и её соседях
    void print_neibs_info() const;


    /// @brief Вывести полную информацию о ячейке
    inline void print_info() const { m_it->print_info(); }

    /// @brief Вывести информацию о ячейке в виде python скрипта
    inline void visualize() const { m_it->print_info(); }

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

private:
    AmrStorage::Iterator m_it;  ///< Итератор хранилища
    AmrStorage &m_locals;       ///< Ссылка на хранилище текущего процесса
    AmrStorage &m_aliens;       ///< Ссылка на хранилище соседей
};

} // namespace zephyr::mesh