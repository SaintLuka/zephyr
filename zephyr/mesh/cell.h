#pragma once

#include <zephyr/geom/face.h>

#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/face.h>


namespace zephyr { namespace mesh {

using zephyr::geom::Adjacent;
using zephyr::geom::Vector3d;

/// @brief Ячейка (итератор для доступа к элементам класса Mesh)
/// @details Копирует функции Storage::Item, а также добавляет
/// специфические функции для ячейки: volume, center, faces.
class ICell : public Storage::Item {
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
    ICell(Storage &locals, Storage &aliens, int idx);

    /// @brief Вспомогательный конструктор. Взять смежную ячейку из locals
    /// или aliens хранилища по индексу смежности с грани.
    /// @param adj Индекс смежности ячейки.
    ICell(Storage &locals, Storage &aliens, const Adjacent &adj);

    /// @brief Разыменование итератора (для цикла for)
    ICell &operator*() final {
        return *this;
    }

    /// @brief Следующая ячейка
    ICell &operator++() final {
        m_ptr += m_itemsize;
        return *this;
    }

    /// @brief Следующая ячейка
    ICell &operator+=(int step) final {
        m_ptr += step * m_itemsize;
        return *this;
    }

    /// @brief Следующая ячейка
    ICell operator+(int step) {
        ICell res(*this);
        res += step;
        return res;
    }

    /// @brief Массив граней
    IFaces faces() const;

    /// @brief Получить смежную ячейку через грань
    ICell neib(const Face &face) const;

    /// @brief Получить смежную ячейку через грань
    ICell neighbor(const Face &face) const;

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

private:
    Storage &locals;  ///< Ссылка на хранилище текущего процесса
    Storage &aliens;  ///< Ссылка на хранилище соседей
};

} // namespace mesh
} // namespace zephyr