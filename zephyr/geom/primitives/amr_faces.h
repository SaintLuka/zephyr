#pragma once

#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/amr_face.h>

namespace zephyr::geom {

/// @brief Список граней ячейки
class AmrFaces {
public:
    /// @brief Максимальное число граней
    static const int max_count = 24;

    /// @brief Конструктор по умолчанию
    AmrFaces() = default;

    /// @brief Конструктор инициализирует индексы граней
    AmrFaces(int dim);

    /// @brief Число актуальных граней в списке
    int count() const;

    /// @brief Начало списка
    inline std::array<AmrFace, max_count>::iterator begin() {
        return m_list.begin();
    }

    /// @brief Начало списка
    inline std::array<AmrFace, max_count>::const_iterator begin() const {
        return m_list.begin();
    }

    /// @brief Конец списка
    inline std::array<AmrFace, max_count>::iterator end() {
        return m_list.end();
    }

    /// @brief Конец списка
    inline std::array<AmrFace, max_count>::const_iterator end() const {
        return m_list.end();
    }

    /// @brief Доступ к грани по индексу
    inline AmrFace &operator[](int idx) {
        return m_list[idx];
    }

    /// @brief Доступ к грани по индексу
    inline const AmrFace &operator[](int idx) const {
        return m_list[idx];
    }

private:
    /// @brief Массив граней ячейки
    std::array<AmrFace, max_count> m_list;
};

} // namespace zephyr::geom