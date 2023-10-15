#pragma once

#include <zephyr/geom/primitives/amr_face.h>

namespace zephyr { namespace geom {

/// @brief Именованные индексы граней ячейки
enum Side : int {
    LEFT    = 0,  L = 0,
    RIGHT   = 1,  R = 1,
    BOTTOM  = 2,  B = 2,
    TOP     = 3,  T = 3,
    BACK    = 4,  X = 4,
    FRONT   = 5,  F = 5,
    LEFT0   = 0,  LEFT1   =  6,  LEFT2   = 12,  LEFT3   = 18,
    RIGHT0  = 1,  RIGHT1  =  7,  RIGHT2  = 13,  RIGHT3  = 19,
    BOTTOM0 = 2,  BOTTOM1 =  8,  BOTTOM2 = 14,  BOTTOM3 = 20,
    TOP0    = 3,  TOP1    =  9,  TOP2    = 15,  TOP3    = 21,
    BACK0   = 4,  BACK1   = 10,  BACK2   = 16,  BACK3   = 22,
    FRONT0  = 5,  FRONT1  = 11,  FRONT2  = 17,  FRONT3  = 23
};

/// @brief Грань в название
inline std::string side_to_string(Side s) {
    switch (s) {
        case Side::LEFT:   return "left";
        case Side::RIGHT:  return "right";
        case Side::BOTTOM: return "bottom";
        case Side::TOP:    return "top";
        case Side::BACK:   return "back";
        case Side::FRONT:  return "front";
        default: return "other";
    }
}

/// @brief Индекс грани в строку
inline std::string side_to_string(int s) {
    return side_to_string(Side(s));
}

/// @brief Требуемый порядок граней
static_assert(Side::LEFT   == 0, "Left   != 0");
static_assert(Side::RIGHT  == 1, "Right  != 1");
static_assert(Side::BOTTOM == 2, "Bottom != 2");
static_assert(Side::TOP    == 3, "Top    != 3");
static_assert(Side::BACK   == 4, "Back   != 4");
static_assert(Side::FRONT  == 5, "Front  != 5");

inline constexpr short undef_index() {
    return -1;
}

template <int dim, Side side>
inline std::array<short, 4> face_indices();

template <>
inline std::array<short, 4> face_indices<2, Side::LEFT>() {
    return {iv(0, 0), iv(0, 1), undef_index(), undef_index()};
}

template <>
inline std::array<short, 4> face_indices<2, Side::RIGHT>() {
    return {iv(1, 0), iv(1, 1), undef_index(), undef_index()};
}

template <>
inline std::array<short, 4> face_indices<2, Side::BOTTOM>() {
    return  {iv(0, 0), iv(1, 0), undef_index(), undef_index()};
}

template <>
inline std::array<short, 4> face_indices<2, Side::TOP>() {
    return {iv(0, 1), iv(1, 1), undef_index(), undef_index()};
}

template <>
inline std::array<short, 4> face_indices<3, Side::LEFT>() {
    return {iv(0, 0, 0), iv(0, 1, 0), iv(0, 0, 1), iv(0, 1, 1)};
}

template <>
inline std::array<short, 4> face_indices<3, Side::RIGHT>() {
    return {iv(1, 0, 0), iv(1, 1, 0), iv(1, 0, 1), iv(1, 1, 1)};
}

template <>
inline std::array<short, 4> face_indices<3, Side::BOTTOM>() {
    return {iv(0, 0, 0), iv(1, 0, 0), iv(0, 0, 1), iv(1, 0, 1)};
}

template <>
inline std::array<short, 4> face_indices<3, Side::TOP>() {
    return {iv(0, 1, 0), iv(1, 1, 0), iv(0, 1, 1), iv(1, 1, 1)};
}

template <>
inline std::array<short, 4> face_indices<3, Side::BACK>() {
    return {iv(0, 0, 0), iv(1, 0, 0), iv(0, 1, 0), iv(1, 1, 0)};
}

template <>
inline std::array<short, 4> face_indices<3, Side::FRONT>() {
    return {iv(0, 0, 1), iv(1, 0, 1), iv(0, 1, 1), iv(1, 1, 1)};
}

/// @brief Список граней ячейки
struct AmrFaces {

    /// @brief Максимальное число граней
    static const int max_size = 24;

    /// @brief Конструктор по умолчанию
    AmrFaces() = default;

    /// @brief Число актуальных граней в списке
    int size() const {
        int count = 0;
        for (int i = 0; i < max_size; ++i) {
            if (m_list[i].is_actual()) {
                ++count;
            }
        }
        return count;
    }

    /// @brief Массив граней
    /// @details Следует исключить
    std::array<AmrFace, max_size>& list() {
        return m_list;
    }

    /// @brief Установить неопределенные грани
    void set_undefined() {
        for (AmrFace &f: m_list) {
            f.set_undefined();
        }
    }

    /// @brief Начало списка
    inline std::array<AmrFace, max_size>::iterator begin() {
        return m_list.begin();
    }

    /// @brief Начало списка
    inline std::array<AmrFace, max_size>::const_iterator begin() const {
        return m_list.begin();
    }

    /// @brief Конец списка
    inline std::array<AmrFace, max_size>::iterator end() {
        return m_list.end();
    }

    /// @brief Конец списка
    inline std::array<AmrFace, max_size>::const_iterator end() const {
        return m_list.end();
    }

    /// @brief Доступ к грани по индексу
    inline AmrFace &operator[](int i) {
        return m_list[i];
    }

    /// @brief Доступ к грани по индексу
    inline const AmrFace &operator[](int i) const {
        return m_list[i];
    }

    template <int dim, Side side>
    inline static std::array<short, 4> indices() {
        return face_indices<dim, side>();
    }

private:
    /// @brief Массив граней ячейки
    std::array<AmrFace, max_size> m_list;
};

} // geom
} // zephyr