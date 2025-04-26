#pragma once

#include <array>
#include <string>

namespace zephyr::mesh {

template <int dim, typename = void>
struct Side { };

/// @brief Именованные индексы граней 2D AMR-ячейки.
/// Константы формируются следующим образом:
/// LEFT[0] = LEFT = L,
/// LEFT[1] = LEFT + 4
/// Для остальных аналогично.
template <int dim>
struct Side<dim, std::enable_if_t<dim == 2>> {
    // Названия граней
    static constexpr Side<dim> LEFT   = 0;
    static constexpr Side<dim> RIGHT  = 1;
    static constexpr Side<dim> BOTTOM = 2;
    static constexpr Side<dim> TOP    = 3;

    // Короткие синонимы
    static constexpr Side<dim> L = LEFT;
    static constexpr Side<dim> R = RIGHT;
    static constexpr Side<dim> B = BOTTOM;
    static constexpr Side<dim> T = TOP;

    int val;

    constexpr Side(int val_) : val(val_) { }

    constexpr operator int() const { return val; }

    constexpr Side<dim> operator[](int idx) const {
        return Side<dim>(val + 4 * idx);
    }

    /*
    LEFT[0   = 0,  LEFT[1]   =  4,
    RIGHT[0]  = 1,  RIGHT[1]  =  5,
    BOTTOM[0] = 2,  BOTTOM[1] =  6,
    TOP[0]    = 3,  TOP[1]    =  7
    */
};

/// @brief Именованные индексы граней 3D AMR-ячейки.
/// Константы формируются следующим образом:
/// LEFT[0] = LEFT = L,
/// LEFT[1] = LEFT + 6,
/// LEFT[2] = LEFT + 12,
/// LEFT[3] = LEFT + 18.
/// Для остальных аналогично.
/// Сокращенное название для грани BACK -- X.
template <int dim>
struct Side<dim, std::enable_if_t<dim == 3>> {
    // Полные имена
    static constexpr Side<dim> LEFT   = 0;
    static constexpr Side<dim> RIGHT  = 1;
    static constexpr Side<dim> BOTTOM = 2;
    static constexpr Side<dim> TOP    = 3;
    static constexpr Side<dim> BACK   = 4;
    static constexpr Side<dim> FRONT  = 5;

    // Сокращенные
    static constexpr Side<dim> L = 0;
    static constexpr Side<dim> R = 1;
    static constexpr Side<dim> B = 2;
    static constexpr Side<dim> T = 3;
    static constexpr Side<dim> X = 4;
    static constexpr Side<dim> F = 5;

    // Значение смещения
    int val;

    // Неявное преобразование из int
    constexpr Side(int val_) : val(val_) { }

    // Неявное преобразование в int
    constexpr operator int() const { return val; }

    // Подграни для сложных граней
    constexpr Side<dim> operator[](int idx) const {
        return Side<dim>(val + 6 * idx);
    }

    /*
    LEFT[0]   = 0,  LEFT[1]   =  6,  LEFT[2]   = 12,  LEFT[3]   = 18,
    RIGHT[0]  = 1,  RIGHT[1]  =  7,  RIGHT[2]  = 13,  RIGHT[3]  = 19,
    BOTTOM[0] = 2,  BOTTOM[1] =  8,  BOTTOM[2] = 14,  BOTTOM[3] = 20,
    TOP[0]    = 3,  TOP[1]    =  9,  TOP[2]    = 15,  TOP[3]    = 21,
    BACK[0]   = 4,  BACK1   = 10,  BACK2   = 16,  BACK3   = 22,
    FRONT0  = 5,  FRONT1  = 11,  FRONT2  = 17,  FRONT3  = 23
    */
};

using Side2D = Side<2>;
using Side3D = Side<3>;

/// @brief Грань в название, часто используется в дебаге
inline std::string side_to_string(Side3D s) {
    std::string num = " face (" + std::to_string(s / 6) + ")";
    switch (s % 6) {
        case Side3D::LEFT:   return "left" + num;
        case Side3D::RIGHT:  return "right" + num;
        case Side3D::BOTTOM: return "bottom" + num;
        case Side3D::TOP:    return "top" + num;
        case Side3D::BACK:   return "back" + num;
        case Side3D::FRONT:  return "front" + num;
        default: return "other" + num;
    }
}

/// @brief Индекс грани в строку
inline std::string side_to_string(int s) {
    return side_to_string(Side3D(s));
}

/// @brief Требуемый порядок граней
static_assert(Side3D::LEFT   == 0, "Left   != 0");
static_assert(Side3D::RIGHT  == 1, "Right  != 1");
static_assert(Side3D::BOTTOM == 2, "Bottom != 2");
static_assert(Side3D::TOP    == 3, "Top    != 3");
static_assert(Side3D::BACK   == 4, "Back   != 4");
static_assert(Side3D::FRONT  == 5, "Front  != 5");

/// @brief Передаем нормальный вектор, функция выбирает сторону
template <int i, int j, int k = 0>
inline constexpr Side3D side_by_dir();

template <>
inline constexpr Side3D side_by_dir<-1, 0, 0>() { return Side3D::LEFT; }

template <>
inline constexpr Side3D side_by_dir<+1, 0, 0>() { return Side3D::RIGHT; }

template <>
inline constexpr Side3D side_by_dir<0, -1, 0>() { return Side3D::BOTTOM; }

template <>
inline constexpr Side3D side_by_dir<0, +1, 0>() { return Side3D::TOP; }

template <>
inline constexpr Side3D side_by_dir<0, 0, -1>() { return Side3D::BACK; }

template <>
inline constexpr Side3D side_by_dir<0, 0, +1>() { return Side3D::FRONT; }


/// @brief Для AMR-ячейки вершины и грани упорядочены определеным образом,
/// поэтому индексы вершин на гранях заведомо известны.
namespace face_indices {

/// @brief Неопределенный индекс
inline constexpr int undef() {
    return -1;
}

/// @brief Индексы вершин простой грани в AMR-ячейке
/// sf -- Simple Face
template<int dim, int side>
inline constexpr std::array<int, 4> sf() {
    return {undef(), undef(), undef(), undef()};
}

/// @brief Индексы вершин части составной грани AMR-ячейки
/// cf -- Complex Face
template<int dim, int side>
inline constexpr std::array<int, 4> cf() {
    return {undef(), undef(), undef(), undef()};
}

/// @brief Индексация вершин в кубической AMR-ячейке
/// Повторяет индексацию функции geom::SqCube::iss
template<int i, int j, int k = -1>
inline int constexpr iss() {
    static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1,
                  "Available indices: {-1, 0, 1}");
    return 9 * (k + 1) + 3 * (j + 1) + i + 1;
}

template<>
inline constexpr std::array<int, 4> sf<2, Side3D::LEFT>() {
    return {iss<-1, -1>(), iss<-1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side3D::RIGHT>() {
    return {iss<+1, -1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side3D::BOTTOM>() {
    return {iss<-1, -1>(), iss<+1, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side3D::TOP>() {
    return {iss<-1, +1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::LEFT>() {
    return {iss<-1, -1, -1>(), iss<-1, +1, -1>(), iss<-1, -1, +1>(), iss<-1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::RIGHT>() {
    return {iss<+1, -1, -1>(), iss<+1, +1, -1>(), iss<+1, -1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::BOTTOM>() {
    return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, -1, +1>(), iss<+1, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::TOP>() {
    return {iss<-1, +1, -1>(), iss<+1, +1, -1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::BACK>() {
    return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, +1, -1>(), iss<+1, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side3D::FRONT>() {
    return {iss<-1, -1, +1>(), iss<+1, -1, +1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::LEFT>() {
    return {iss<-1, -1>(), iss<-1, 0>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::LEFT[1]>() {
    return {iss<-1, 0>(), iss<-1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::RIGHT[0]>() {
    return {iss<+1, -1>(), iss<+1, 0>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::RIGHT[1]>() {
    return {iss<+1, 0>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::BOTTOM[0]>() {
    return {iss<-1, -1>(), iss<0, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::BOTTOM[1]>() {
    return {iss<0, -1>(), iss<+1, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::TOP[0]>() {
    return {iss<-1, +1>(), iss<0, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side3D::TOP[1]>() {
    return {iss<0, +1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::LEFT[0]>() {
    return {iss<-1, -1, -1>(), iss<-1, 0, -1>(), iss<-1, -1, 0>(), iss<-1, 0, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::LEFT[1]>() {
    return {iss<-1, 0, -1>(), iss<-1, +1, -1>(), iss<-1, 0, 0>(), iss<-1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::LEFT[2]>() {
    return {iss<-1, -1, 0>(), iss<-1, 0, 0>(), iss<-1, -1, +1>(), iss<-1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::LEFT[3]>() {
    return {iss<-1, 0, 0>(), iss<-1, +1, 0>(), iss<-1, 0, +1>(), iss<-1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::RIGHT[0]>() {
    return {iss<+1, -1, -1>(), iss<+1, 0, -1>(), iss<+1, -1, 0>(), iss<+1, 0, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::RIGHT[1]>() {
    return {iss<+1, 0, -1>(), iss<+1, +1, -1>(), iss<+1, 0, 0>(), iss<+1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::RIGHT[2]>() {
    return {iss<+1, -1, 0>(), iss<+1, 0, 0>(), iss<+1, -1, +1>(), iss<+1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::RIGHT[3]>() {
    return {iss<+1, 0, 0>(), iss<+1, +1, 0>(), iss<+1, 0, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BOTTOM[0]>() {
    return {iss<-1, -1, -1>(), iss<0, -1, -1>(), iss<-1, -1, 0>(), iss<0, -1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BOTTOM[1]>() {
    return {iss<0, -1, -1>(), iss<+1, -1, -1>(), iss<0, -1, 0>(), iss<+1, -1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BOTTOM[2]>() {
    return {iss<-1, -1, 0>(), iss<0, -1, 0>(), iss<-1, -1, +1>(), iss<0, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BOTTOM[3]>() {
    return {iss<0, -1, 0>(), iss<+1, -1, 0>(), iss<0, -1, +1>(), iss<+1, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::TOP[0]>() {
    return {iss<-1, +1, -1>(), iss<0, +1, -1>(), iss<-1, +1, 0>(), iss<0, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::TOP[1]>() {
    return {iss<0, +1, -1>(), iss<+1, +1, -1>(), iss<0, +1, 0>(), iss<+1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::TOP[2]>() {
    return {iss<-1, +1, 0>(), iss<0, +1, 0>(), iss<-1, +1, +1>(), iss<0, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::TOP[3]>() {
    return {iss<0, +1, 0>(), iss<+1, +1, 0>(), iss<0, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BACK[0]>() {
    return {iss<-1, -1, -1>(), iss<0, -1, -1>(), iss<-1, 0, -1>(), iss<0, 0, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BACK[1]>() {
    return {iss<0, -1, -1>(), iss<+1, -1, -1>(), iss<0, 0, -1>(), iss<+1, 0, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BACK[2]>() {
    return {iss<-1, 0, -1>(), iss<0, 0, -1>(), iss<-1, +1, -1>(), iss<0, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::BACK[3]>() {
    return {iss<0, 0, -1>(), iss<+1, 0, -1>(), iss<0, +1, -1>(), iss<+1, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::FRONT[0]>() {
    return {iss<-1, -1, +1>(), iss<0, -1, +1>(), iss<-1, 0, +1>(), iss<0, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::FRONT[1]>() {
    return {iss<0, -1, +1>(), iss<+1, -1, +1>(), iss<0, 0, +1>(), iss<+1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::FRONT[2]>() {
    return {iss<-1, 0, +1>(), iss<0, 0, +1>(), iss<-1, +1, +1>(), iss<0, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side3D::FRONT[3]>() {
    return {iss<0, 0, +1>(), iss<+1, 0, +1>(), iss<0, +1, +1>(), iss<+1, +1, +1>()};
}

} // namespace face_indices

} // namespace zephyr::mesh