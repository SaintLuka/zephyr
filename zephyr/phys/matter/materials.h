#pragma once

#include <vector>

#include <zephyr/phys/matter/material.h>
#include "mixture_pt.h"

namespace zephyr::phys {

/// @brief Простой контейнер для хранения массива материалов
/// Почти ничем не отличается от std::vector кроме пары доп функций.
class Materials {
public:

    /// @brief Конструктор по умолчанию (пустой список)
    Materials() = default;

    /// @brief Создать из списка инициалиации
    Materials(const std::initializer_list<Material>& il)
            : m_materials(il) { }

    /// @brief Число компонент
    inline int size() const { return m_materials.size(); }

    /// @brief Единственный материал?
    inline bool single() const { return m_materials.size() == 1; }

    /// @brief Удалить материалы
    void clear();

    /// @brief Добавить материал в список
    void append(Material& mat);

    /// @brief Добавить материал в список
    void operator+=(Material& mat);

    /// @brief Добавить материал в список
    void append(Eos::Ref eos);

    /// @brief Добавить материал в список
    void operator+=(Eos::Ref eos);

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    inline Material& operator[](int idx) {
        return m_materials[idx];
    }

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    inline const Material& operator[](int idx) const {
        return m_materials[idx];
    }

    /// @brief Уравнение состояния смеси с PT-замыканием
    MixturePT mixture_PT() const;

protected:
    /// @brief Массив материалов
    std::vector<Material> m_materials;
};

} // namespace zephyr::phys