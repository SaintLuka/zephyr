#pragma once

#include <zephyr/mesh/storage.h>
#include <zephyr/mesh/cell.h>

namespace zephyr { namespace mesh {

/// @brief Интерфейс диапазона внутренних ячеек сетки
class Range {
public:

    /// @brief Конструктор класса для работы с распределенными сетками
    /// @param locals Ссылка на локальное хранилище
    /// @param aliens Хранилище для ячеек с удаленных процессов
    /// @param from Индекс начала диапазона в массиве locals
    /// @param to Индекс после конца диапазона в массиве locals
    /// @param id Номер для данного диапазона ячеек
    Range(Storage &cells, Storage &aliens, int from, int to, int id = 0);


    operator Storage&() {
        return m_locals;
    }

    /// @brief Количество ячеек в диапазоне
    int size() const;

    /// @brief Индекс диапазона
    int id() const;

    /// @brief Возращает итератор в начало диапазона
    ICell begin() const;

    /// @brief Возвращает итератор на конец диапазона
    ICell end() const;

private:

    int m_id;
    int m_begin;
    int m_end;

    Storage &m_locals;
    Storage &m_aliens;
};

} // mesh
} // zephyr