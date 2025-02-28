#pragma once

#include <array>
#include <string>

namespace zephyr::mesh {

/// @brief Именованные индексы граней AMR-ячейки.
/// Константы формируются следующим образом:
/// LEFT0 = LEFT = L,
/// LEFT1 = LEFT0 +  6,
/// LEFT2 = LEFT0 + 12
/// LEFT3 = LEFT0 + 18
/// Для остальных аналогично.
/// Сокращенное название для грани BACK -- X.
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

/// @brief Грань в название, часто используется в дебаге
inline std::string side_to_string(Side s) {
    std::string num = " face (" + std::to_string(s / 6) + ")";
    switch (s % 6) {
        case Side::LEFT:   return "left" + num;
        case Side::RIGHT:  return "right" + num;
        case Side::BOTTOM: return "bottom" + num;
        case Side::TOP:    return "top" + num;
        case Side::BACK:   return "back" + num;
        case Side::FRONT:  return "front" + num;
        default: return "other" + num;
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

/// @brief Передаем нормальный вектор, функция выбирает сторону
template <int i, int j, int k = 0>
inline constexpr Side side_by_dir();

template <>
inline constexpr Side side_by_dir<-1, 0, 0>() { return Side::LEFT; }

template <>
inline constexpr Side side_by_dir<+1, 0, 0>() { return Side::RIGHT; }

template <>
inline constexpr Side side_by_dir<0, -1, 0>() { return Side::BOTTOM; }

template <>
inline constexpr Side side_by_dir<0, +1, 0>() { return Side::TOP; }

template <>
inline constexpr Side side_by_dir<0, 0, -1>() { return Side::BACK; }

template <>
inline constexpr Side side_by_dir<0, 0, +1>() { return Side::FRONT; }


/// @brief Для AMR-ячейки вершины и грани упорядочены определеным образом,
/// поэтому индексы вершин на гранях заведомо известны.
namespace face_indices {

/// @brief Неопределенный индекс
inline constexpr int undef() {
    return -1;
}

/// @brief Индексы вершин простой грани в AMR-ячейке
/// sf -- Simple Face
template<int dim, Side side>
inline constexpr std::array<int, 4> sf() {
    return {undef(), undef(), undef(), undef()};
}

/// @brief Индексы вершин части составной грани AMR-ячейки
/// cf -- Complex Face
template<int dim, Side side>
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
inline constexpr std::array<int, 4> sf<2, Side::LEFT>() {
    return {iss<-1, -1>(), iss<-1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side::RIGHT>() {
    return {iss<+1, -1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side::BOTTOM>() {
    return {iss<-1, -1>(), iss<+1, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<2, Side::TOP>() {
    return {iss<-1, +1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::LEFT>() {
    return {iss<-1, -1, -1>(), iss<-1, +1, -1>(), iss<-1, -1, +1>(), iss<-1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::RIGHT>() {
    return {iss<+1, -1, -1>(), iss<+1, +1, -1>(), iss<+1, -1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::BOTTOM>() {
    return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, -1, +1>(), iss<+1, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::TOP>() {
    return {iss<-1, +1, -1>(), iss<+1, +1, -1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::BACK>() {
    return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, +1, -1>(), iss<+1, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> sf<3, Side::FRONT>() {
    return {iss<-1, -1, +1>(), iss<+1, -1, +1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::LEFT>() {
    return {iss<-1, -1>(), iss<-1, 0>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::LEFT1>() {
    return {iss<-1, 0>(), iss<-1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::RIGHT0>() {
    return {iss<+1, -1>(), iss<+1, 0>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::RIGHT1>() {
    return {iss<+1, 0>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::BOTTOM0>() {
    return {iss<-1, -1>(), iss<0, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::BOTTOM1>() {
    return {iss<0, -1>(), iss<+1, -1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::TOP0>() {
    return {iss<-1, +1>(), iss<0, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<2, Side::TOP1>() {
    return {iss<0, +1>(), iss<+1, +1>(), undef(), undef()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::LEFT0>() {
    return {iss<-1, -1, -1>(), iss<-1, 0, -1>(), iss<-1, -1, 0>(), iss<-1, 0, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::LEFT1>() {
    return {iss<-1, 0, -1>(), iss<-1, +1, -1>(), iss<-1, 0, 0>(), iss<-1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::LEFT2>() {
    return {iss<-1, -1, 0>(), iss<-1, 0, 0>(), iss<-1, -1, +1>(), iss<-1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::LEFT3>() {
    return {iss<-1, 0, 0>(), iss<-1, +1, 0>(), iss<-1, 0, +1>(), iss<-1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::RIGHT0>() {
    return {iss<+1, -1, -1>(), iss<+1, 0, -1>(), iss<+1, -1, 0>(), iss<+1, 0, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::RIGHT1>() {
    return {iss<+1, 0, -1>(), iss<+1, +1, -1>(), iss<+1, 0, 0>(), iss<+1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::RIGHT2>() {
    return {iss<+1, -1, 0>(), iss<+1, 0, 0>(), iss<+1, -1, +1>(), iss<+1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::RIGHT3>() {
    return {iss<+1, 0, 0>(), iss<+1, +1, 0>(), iss<+1, 0, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BOTTOM0>() {
    return {iss<-1, -1, -1>(), iss<0, -1, -1>(), iss<-1, -1, 0>(), iss<0, -1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BOTTOM1>() {
    return {iss<0, -1, -1>(), iss<+1, -1, -1>(), iss<0, -1, 0>(), iss<+1, -1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BOTTOM2>() {
    return {iss<-1, -1, 0>(), iss<0, -1, 0>(), iss<-1, -1, +1>(), iss<0, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BOTTOM3>() {
    return {iss<0, -1, 0>(), iss<+1, -1, 0>(), iss<0, -1, +1>(), iss<+1, -1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::TOP0>() {
    return {iss<-1, +1, -1>(), iss<0, +1, -1>(), iss<-1, +1, 0>(), iss<0, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::TOP1>() {
    return {iss<0, +1, -1>(), iss<+1, +1, -1>(), iss<0, +1, 0>(), iss<+1, +1, 0>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::TOP2>() {
    return {iss<-1, +1, 0>(), iss<0, +1, 0>(), iss<-1, +1, +1>(), iss<0, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::TOP3>() {
    return {iss<0, +1, 0>(), iss<+1, +1, 0>(), iss<0, +1, +1>(), iss<+1, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BACK0>() {
    return {iss<-1, -1, -1>(), iss<0, -1, -1>(), iss<-1, 0, -1>(), iss<0, 0, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BACK1>() {
    return {iss<0, -1, -1>(), iss<+1, -1, -1>(), iss<0, 0, -1>(), iss<+1, 0, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BACK2>() {
    return {iss<-1, 0, -1>(), iss<0, 0, -1>(), iss<-1, +1, -1>(), iss<0, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::BACK3>() {
    return {iss<0, 0, -1>(), iss<+1, 0, -1>(), iss<0, +1, -1>(), iss<+1, +1, -1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::FRONT0>() {
    return {iss<-1, -1, +1>(), iss<0, -1, +1>(), iss<-1, 0, +1>(), iss<0, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::FRONT1>() {
    return {iss<0, -1, +1>(), iss<+1, -1, +1>(), iss<0, 0, +1>(), iss<+1, 0, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::FRONT2>() {
    return {iss<-1, 0, +1>(), iss<0, 0, +1>(), iss<-1, +1, +1>(), iss<0, +1, +1>()};
}

template<>
inline constexpr std::array<int, 4> cf<3, Side::FRONT3>() {
    return {iss<0, 0, +1>(), iss<+1, 0, +1>(), iss<0, +1, +1>(), iss<+1, +1, +1>()};
}

} // namespace face_indices

} // namespace zephyr::mesh