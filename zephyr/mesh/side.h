#pragma once

#include <array>
#include <string>

namespace zephyr::mesh {

/// @brief Число вершин простой грани.
constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Количество подграней.
constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число граней у простой ячейки.
constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }


/// @brief Структура позволяет хранить индексы граней AMR-ячейки.
/// В большинстве контекстов работает как простой `enum Side : int`,
/// то есть позволяет получить `Side2D::LEFT`, `Side3D::BACK`, ...
/// Также позволяет получить индексы подграней, к примеру,
/// для трехмерной ячейки `Side3D::BACK[0]`, `Side3D::BACK[1]`, ...
/// @tparam dim Размерность ячейки
template <int dim>
struct Side final {
    // Только двумерные и трехмерные
    static_assert(dim == 2 || dim == 3);

    /// @brief Именованные индексы простых граней
    static const Side LEFT;    //=0
    static const Side RIGHT;   //=1
    static const Side BOTTOM;  //=2
    static const Side TOP;     //=3
    static const Side BACK;    //=4; Только для dim == 3.
    static const Side FRONT;   //=5; Только для dim == 3.

    /// @brief Сокращенные имена граней (Z = BACK)
    static const Side L, R, B, T, Z, F;

    /// @brief Индекс грани
    int val;

    /// @brief Неявное создание из типа int
    constexpr Side(int v) noexcept : val(v) { }

    /// @brief Неявное преобразование в int
    constexpr operator int() const { return val; }

    /// @brief Число простых граней ячейки
    static constexpr int count() { return 2 * dim; }

    /// @brief Число подграней сложной грани
    static constexpr int n_subfaces() { return dim == 2 ? 2 : 4; }

    /// @brief Индексы подграней на стороне в списке faces.
    std::array<Side, n_subfaces()> subfaces() const {
        if constexpr (dim == 2) {
            return {Side{val}, Side{val}[1]};
        } else {
            return {Side{val}, Side{val}[1], Side{val}[2], Side{val}[3]};
        }
    }

    /// @brief Получить подгрань сложной грани
    constexpr Side operator[](int idx) const {
        return Side(val + count() * idx);
    }

    /// @brief Перечисление простых граней двумерной ячейки
    static constexpr std::array<Side, count()> items() {
        if constexpr (dim == 2) {
            return {LEFT, RIGHT, BOTTOM, TOP};
        } else {
            return {LEFT, RIGHT, BOTTOM, TOP, BACK, FRONT};
        }
    }

    /// @brief Получить сторону по нормальному вектору направления
    /// Side3D::LEFT == Side3D::by_dir<-1, 0, 0>().
    template <int i, int j, int k = 0>
    static constexpr Side by_dir() {
        if constexpr (i < 0 && j == 0 && k == 0) { return LEFT; }
        if constexpr (i > 0 && j == 0 && k == 0) { return RIGHT; }
        if constexpr (i == 0 && j < 0 && k == 0) { return BOTTOM; }
        if constexpr (i == 0 && j > 0 && k == 0) { return TOP; }
        if constexpr (i == 0 && j == 0 && k < 0) { return BACK; }
        if constexpr (i == 0 && j == 0 && k > 0) { return FRONT; }
        return -1;
    }

    /// @brief Преобразование в строку
    std::string to_string() const {
        std::string num = "[" + std::to_string(val / count()) + "]";
        switch (val % count()) {
            case 0: return "LEFT" + num;
            case 1: return "RIGHT" + num;
            case 2: return "BOTTOM" + num;
            case 3: return "TOP" + num;
            case 4: return "BACK" + num;
            case 5: return "FRONT" + num;
            default: return "STRANGE FACE";
        }
    }

    // Для AMR-ячейки вершины и грани упорядочены определенным образом,
    // поэтому индексы вершин на гранях заведомо известны.

    /// @brief Индексы вершин простой грани в AMR-ячейке
    /// sf -- Simple Face
    constexpr std::array<short, 8> sf() const {
        if constexpr (dim == 2) {
            switch (val) {
                case L: return {iss<-1, -1>(), iss<-1, +1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case R: return {iss<+1, -1>(), iss<+1, +1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case B: return {iss<-1, -1>(), iss<+1, -1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case T: return {iss<-1, +1>(), iss<+1, +1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                default: return {undef(), undef(), undef(), undef()};
            }
        } else {
            switch (val) {
                case L: return {iss<-1, -1, -1>(), iss<-1, +1, -1>(), iss<-1, -1, +1>(), iss<-1, +1, +1>(), undef(), undef(), undef(), undef()};
                case R: return {iss<+1, -1, -1>(), iss<+1, +1, -1>(), iss<+1, -1, +1>(), iss<+1, +1, +1>(), undef(), undef(), undef(), undef()};
                case B: return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, -1, +1>(), iss<+1, -1, +1>(), undef(), undef(), undef(), undef()};
                case T: return {iss<-1, +1, -1>(), iss<+1, +1, -1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>(), undef(), undef(), undef(), undef()};
                case Z: return {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, +1, -1>(), iss<+1, +1, -1>(), undef(), undef(), undef(), undef()};
                case F: return {iss<-1, -1, +1>(), iss<+1, -1, +1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>(), undef(), undef(), undef(), undef()};
                default: return {undef(), undef(), undef(), undef(), undef(), undef(), undef(), undef()};
            }
        }
    }

    /// @brief Индексы вершин части составной грани AMR-ячейки
    /// cf -- Complex Face
    constexpr std::array<short, 8> cf() const {
        if constexpr (dim == 2) {
            switch (val) {
                case L[0]: return {iss<-1,-1>(), iss<-1, 0>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case L[1]: return {iss<-1, 0>(), iss<-1,+1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case R[0]: return {iss<+1,-1>(), iss<+1, 0>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case R[1]: return {iss<+1, 0>(), iss<+1,+1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case B[0]: return {iss<-1,-1>(), iss< 0,-1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case B[1]: return {iss< 0,-1>(), iss<+1,-1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case T[0]: return {iss<-1,+1>(), iss< 0,+1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                case T[1]: return {iss< 0,+1>(), iss<+1,+1>(), undef(), undef(), undef(), undef(), undef(), undef()};
                default:   return {undef(), undef(), undef(), undef(), undef(), undef(), undef(), undef()};
            }
        }
        else {
            switch (val) {
                case L[0]: return {iss<-1,-1,-1>(), iss<-1, 0,-1>(), iss<-1,-1, 0>(), iss<-1, 0, 0>(), undef(), undef(), undef(), undef()};
                case L[1]: return {iss<-1, 0,-1>(), iss<-1,+1,-1>(), iss<-1, 0, 0>(), iss<-1,+1, 0>(), undef(), undef(), undef(), undef()};
                case L[2]: return {iss<-1,-1, 0>(), iss<-1, 0, 0>(), iss<-1,-1,+1>(), iss<-1, 0,+1>(), undef(), undef(), undef(), undef()};
                case L[3]: return {iss<-1, 0, 0>(), iss<-1,+1, 0>(), iss<-1, 0,+1>(), iss<-1,+1,+1>(), undef(), undef(), undef(), undef()};
                case R[0]: return {iss<+1,-1,-1>(), iss<+1, 0,-1>(), iss<+1,-1, 0>(), iss<+1, 0, 0>(), undef(), undef(), undef(), undef()};
                case R[1]: return {iss<+1, 0,-1>(), iss<+1,+1,-1>(), iss<+1, 0, 0>(), iss<+1,+1, 0>(), undef(), undef(), undef(), undef()};
                case R[2]: return {iss<+1,-1, 0>(), iss<+1, 0, 0>(), iss<+1,-1,+1>(), iss<+1, 0,+1>(), undef(), undef(), undef(), undef()};
                case R[3]: return {iss<+1, 0, 0>(), iss<+1,+1, 0>(), iss<+1, 0,+1>(), iss<+1,+1,+1>(), undef(), undef(), undef(), undef()};
                case B[0]: return {iss<-1,-1,-1>(), iss< 0,-1,-1>(), iss<-1,-1, 0>(), iss< 0,-1, 0>(), undef(), undef(), undef(), undef()};
                case B[1]: return {iss< 0,-1,-1>(), iss<+1,-1,-1>(), iss< 0,-1, 0>(), iss<+1,-1, 0>(), undef(), undef(), undef(), undef()};
                case B[2]: return {iss<-1,-1, 0>(), iss< 0,-1, 0>(), iss<-1,-1,+1>(), iss< 0,-1,+1>(), undef(), undef(), undef(), undef()};
                case B[3]: return {iss< 0,-1, 0>(), iss<+1,-1, 0>(), iss< 0,-1,+1>(), iss<+1,-1,+1>(), undef(), undef(), undef(), undef()};
                case T[0]: return {iss<-1,+1,-1>(), iss< 0,+1,-1>(), iss<-1,+1, 0>(), iss< 0,+1, 0>(), undef(), undef(), undef(), undef()};
                case T[1]: return {iss< 0,+1,-1>(), iss<+1,+1,-1>(), iss< 0,+1, 0>(), iss<+1,+1, 0>(), undef(), undef(), undef(), undef()};
                case T[2]: return {iss<-1,+1, 0>(), iss< 0,+1, 0>(), iss<-1,+1,+1>(), iss< 0,+1,+1>(), undef(), undef(), undef(), undef()};
                case T[3]: return {iss< 0,+1, 0>(), iss<+1,+1, 0>(), iss< 0,+1,+1>(), iss<+1,+1,+1>(), undef(), undef(), undef(), undef()};
                case Z[0]: return {iss<-1,-1,-1>(), iss<0, -1,-1>(), iss<-1, 0,-1>(), iss< 0, 0,-1>(), undef(), undef(), undef(), undef()};
                case Z[1]: return {iss< 0,-1,-1>(), iss<+1,-1,-1>(), iss< 0, 0,-1>(), iss<+1, 0,-1>(), undef(), undef(), undef(), undef()};
                case Z[2]: return {iss<-1, 0,-1>(), iss< 0, 0,-1>(), iss<-1,+1,-1>(), iss< 0,+1,-1>(), undef(), undef(), undef(), undef()};
                case Z[3]: return {iss< 0, 0,-1>(), iss<+1, 0,-1>(), iss< 0,+1,-1>(), iss<+1,+1,-1>(), undef(), undef(), undef(), undef()};
                case F[0]: return {iss<-1,-1,+1>(), iss< 0,-1,+1>(), iss<-1, 0,+1>(), iss< 0, 0,+1>(), undef(), undef(), undef(), undef()};
                case F[1]: return {iss< 0,-1,+1>(), iss<+1,-1,+1>(), iss< 0, 0,+1>(), iss<+1, 0,+1>(), undef(), undef(), undef(), undef()};
                case F[2]: return {iss<-1, 0,+1>(), iss< 0, 0,+1>(), iss<-1,+1,+1>(), iss< 0,+1,+1>(), undef(), undef(), undef(), undef()};
                case F[3]: return {iss< 0, 0,+1>(), iss<+1, 0,+1>(), iss< 0,+1,+1>(), iss<+1,+1,+1>(), undef(), undef(), undef(), undef()};
                default:   return {undef(), undef(), undef(), undef(), undef(), undef(), undef(), undef()};
            }
        }
    }

private:
    /// @brief Неопределенный индекс вершины
    static constexpr short undef() { return -1; }

    /// @brief Индексация вершин в кубической AMR-ячейке. Повторяет индексацию функции SqCube::iss
    template<short i, short j, short k = -1>
    static constexpr short iss() {
        static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1, "Available indices: {-1, 0, 1}");
        return 9 * (k + 1) + 3 * (j + 1) + i + 1;
    }
};

// Инициализация констант для Side<2>, ошибка компиляции
// при попытке использовать Side2D::BACK
using Side2D = Side<2>;
template <> constexpr Side2D Side2D::LEFT   = 0;
template <> constexpr Side2D Side2D::RIGHT  = 1;
template <> constexpr Side2D Side2D::BOTTOM = 2;
template <> constexpr Side2D Side2D::TOP    = 3;
template <> constexpr Side2D Side2D::L = LEFT;
template <> constexpr Side2D Side2D::R = RIGHT;
template <> constexpr Side2D Side2D::B = BOTTOM;
template <> constexpr Side2D Side2D::T = TOP;

// Инициализация констант для Side<3>
using Side3D = Side<3>;
template <> constexpr Side3D Side3D::LEFT   = 0;
template <> constexpr Side3D Side3D::RIGHT  = 1;
template <> constexpr Side3D Side3D::BOTTOM = 2;
template <> constexpr Side3D Side3D::TOP    = 3;
template <> constexpr Side3D Side3D::BACK   = 4;
template <> constexpr Side3D Side3D::FRONT  = 5;
template <> constexpr Side3D Side3D::L = LEFT;
template <> constexpr Side3D Side3D::R = RIGHT;
template <> constexpr Side3D Side3D::B = BOTTOM;
template <> constexpr Side3D Side3D::T = TOP;
template <> constexpr Side3D Side3D::Z = BACK;
template <> constexpr Side3D Side3D::F = FRONT;

/// @brief Вывод названия грани в поток
template <int dim>
std::ostream& operator<<(std::ostream& os, Side<dim> side) {
    os << side.to_string();
    return os;
}

/// @brief Индекс грани в строку, dim -- размерность
inline std::string side_to_string(int s, int dim) {
    return dim == 2 ? Side2D(s).to_string() : Side3D(s).to_string();
}

/// @brief Требуемый порядок граней
static_assert(Side2D::LEFT   == 0, "Left   != 0");
static_assert(Side2D::RIGHT  == 1, "Right  != 1");
static_assert(Side2D::BOTTOM == 2, "Bottom != 2");
static_assert(Side2D::TOP    == 3, "Top    != 3");

static_assert(Side3D::LEFT   == 0, "Left   != 0");
static_assert(Side3D::RIGHT  == 1, "Right  != 1");
static_assert(Side3D::BOTTOM == 2, "Bottom != 2");
static_assert(Side3D::TOP    == 3, "Top    != 3");
static_assert(Side3D::BACK   == 4, "Back   != 4");
static_assert(Side3D::FRONT  == 5, "Front  != 5");

} // namespace zephyr::mesh