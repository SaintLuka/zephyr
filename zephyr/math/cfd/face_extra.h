#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/limiter.h>

#define sgn(x) (x > 0.0 ? 1.0 : (x == 0.0 ? 0.0 : -1.0))

namespace zephyr::math {

template<class T>
class FaceExtraSimple;

template<class T>
class FaceExtraATvL;

template<class T>
class FaceExtraTriad;

namespace FaceExtra {

using zephyr::geom::Vector3d;

/// @brief Размер структуры или класса в double
template<class T>
static constexpr int SizeOf() {
    return sizeof(T) / sizeof(double);
}

/// @brief Приведение произвольного типа к массиву double
template<class T>
static inline std::array<double, SizeOf<T>()> &arr(T &vec) {
    return reinterpret_cast<std::array<double, SizeOf<T>()> &>(vec);
}

/// @brief Приведение произвольного типа к массиву double
template<class T>
static inline const std::array<double, SizeOf<T>()> &arr(const T &vec) {
    return reinterpret_cast<const std::array<double, SizeOf<T>()> &>(vec);
}

/// @brief Создание простого экстраполятора
/// @param uc Значение в данной ячейке
/// @param ucx, ucy, ucz Производные в данной ячейке
/// @param un Значение в смежной ячейке
/// @param unx, uny, unz Производные в смежной ячейке
/// @param cell_c Центр данной ячейки
/// @param neib_c Центр смежной ячейки
/// @param face_c Центр грани
template<class T>
inline FaceExtraSimple<T> Simple(Limiter &limiter,
                                 const T &uc, const T &ucx, const T &ucy, const T &ucz,
                                 const T &un, const T &unx, const T &uny, const T &unz,
                                 const Vector3d &cell_c, const Vector3d &neib_c,
                                 const Vector3d &face_c) {
    return FaceExtraSimple<T>(limiter,
                              arr(uc), arr(ucx), arr(ucy), arr(ucz),
                              arr(un), arr(unx), arr(uny), arr(unz),
                              cell_c, neib_c, face_c);
}

/// @brief Создание экстраполятора типа ATvL
/// @param uc Значение в данной ячейке
/// @param ucx, ucy, ucz Производные в данной ячейке
/// @param un Значение в смежной ячейке
/// @param unx, uny, unz Производные в смежной ячейке
/// @param cell_c Центр данной ячейки
/// @param neib_c Центр смежной ячейки
/// @param face_c Центр грани
template<class T>
inline FaceExtraATvL<T> ATvL(const T &uc, const T &ucx, const T &ucy, const T &ucz,
                             const T &un, const T &unx, const T &uny, const T &unz,
                             const Vector3d &cell_c, const Vector3d &neib_c,
                             const Vector3d &face_c) {
    return FaceExtraATvL<T>(arr(uc), arr(ucx), arr(ucy), arr(ucz),
                            arr(un), arr(unx), arr(uny), arr(unz),
                            cell_c, neib_c, face_c);
}

/// @brief Создание экстраполятора типа ATvL
/// @param uc Значение в данной ячейке
/// @param ucx, ucy, ucz Производные в данной ячейке
/// @param un Значение в смежной ячейке
/// @param unx, uny, unz Производные в смежной ячейке
/// @param cell_c Центр данной ячейки
/// @param neib_c Центр смежной ячейки
/// @param face_c Центр грани
template<class T>
inline FaceExtraTriad<T> Triad(const T &uc, const T &ucx, const T &ucy, const T &ucz,
                             const T &un, const T &unx, const T &uny, const T &unz,
                             const Vector3d &cell_c, const Vector3d &neib_c,
                             const Vector3d &face_c) {
    return FaceExtraTriad<T>(arr(uc), arr(ucx), arr(ucy), arr(ucz),
                            arr(un), arr(unx), arr(uny), arr(unz),
                            cell_c, neib_c, face_c);
}

} // namespace FaceExtra


/// @brief Простая экстраполяция для схем второго порядка
template<class T>
class FaceExtraSimple {
public:
    static const int N = FaceExtra::SizeOf<T>();
    using Array = std::array<double, N>;

    /// @brief Создание простого экстраполятора
    /// @param uc Значение в данной ячейке
    /// @param ucx, ucy, ucz Производные в данной ячейке
    /// @param un Значение в смежной ячейке
    /// @param unx, uny, unz Производные в смежной ячейке
    /// @param cell_c Центр данной ячейки
    /// @param neib_c Центр смежной ячейки
    /// @param face_c Центр грани
    FaceExtraSimple(Limiter &_limiter,
                    const Array &uc, const Array &ucx, const Array &ucy, const Array &ucz,
                    const Array &un, const Array &unx, const Array &uny, const Array &unz,
                    const Vector3d &cell_c, const Vector3d &neib_c,
                    const Vector3d &face_c) : limiter(_limiter) {

        Vector3d dr = neib_c - cell_c;

        for (int i = 0; i < N; ++i) {
            u_c[i] = uc[i];
            u_n[i] = un[i];
            du_c[i] = ucx[i] * dr.x() + ucy[i] * dr.y() + unz[i] * dr.z();
            du_n[i] = unx[i] * dr.x() + uny[i] * dr.y() + unz[i] * dr.z();
        }
    }

    /// @brief Экстраполяция величины ul на внутреннюю сторону грани
    T m(const T &_ul) const {
        T _um;

        auto &ul = FaceExtra::arr(_ul);
        auto &um = FaceExtra::arr(_um);

        for (int i = 0; i < N; ++i) {
            double dzp = u_n[i] - u_c[i];
            double dzm = 2.0 * du_c[i] - dzp;

            um[i] = ul[i] + 0.5 * limiter(dzm, dzp) * dzp;
        }

        return _um;
    }

    /// @brief Экстраполяция величины ur на внешнюю сторону грани
    T p(const T &_ur) const {
        T _up;

        auto &ur = FaceExtra::arr(_ur);
        auto &up = FaceExtra::arr(_up);

        for (int i = 0; i < N; ++i) {
            double dzm = u_n[i] - u_c[i];
            double dzp = 2.0 * du_n[i] - dzm;

            up[i] = ur[i] - 0.5 * limiter(dzp, dzm) * dzm;
        }

        return _up;
    }

private:

    /// Параметры экстраполяции
    Limiter &limiter;
    Array u_c, u_n;
    Array du_c, du_n;
};

/// @brief Класс для TVD экстраполяции на грань.
/// Anderson, Thomas, van Leer
template<class T>
class FaceExtraATvL {
public:
    static const int N = FaceExtra::SizeOf<T>();
    using Array = std::array<double, N>;

    /// @brief Создание экстраполятора типа ATvL
    /// @param uc Значение в данной ячейке
    /// @param ucx, ucy, ucz Производные в данной ячейке
    /// @param un Значение в смежной ячейке
    /// @param unx, uny, unz Производные в смежной ячейке
    /// @param cell_c Центр данной ячейки
    /// @param neib_c Центр смежной ячейки
    /// @param face_c Центр грани
    FaceExtraATvL(const Array &uc, const Array &ucx, const Array &ucy, const Array &ucz,
                  const Array &un, const Array &unx, const Array &uny, const Array &unz,
                  const Vector3d &cell_c, const Vector3d &neib_c,
                  const Vector3d &face_c) {

        Vector3d rc = face_c - cell_c;
        Vector3d rn = neib_c - face_c;

        dx_c = rc.norm();
        dx_n = rn.norm();

        for (int i = 0; i < N; ++i) {
            u_c[i] = uc[i];
            u_n[i] = un[i];
            du_c[i] = (ucx[i] * rc.x() + ucy[i] * rc.y() + ucz[i] * rc.z()) / dx_c;
            du_n[i] = (unx[i] * rn.x() + uny[i] * rn.y() + unz[i] * rn.z()) / dx_n;
        }
    }

    /// @brief Экстраполяция величины ul на внутреннюю сторону грани
    T m(const T &_ul) const {
        T _um;

        auto &ul = FaceExtra::arr(_ul);
        auto &um = FaceExtra::arr(_um);

        double dx = dx_c + dx_n;
        double delta = dx_c / dx;
        double eta = (12.0 * delta * delta - 1.0) / 24.0;

        for (int i = 0; i < N; ++i) {
            double dzp = u_n[i] - u_c[i];
            double dzm = 2.0 * du_c[i] * dx - dzp;
            double s = limiters::vanAlbada2(dzm, dzp);

            um[i] = ul[i] + s * ((0.5 * delta - s * eta) * dzm + (0.5 * delta + s * eta) * dzp);
        }

        return _um;
    }

    /// @brief Экстраполяция величины ur на внешнюю сторону грани
    T p(const T &_ur) const {
        T _up;

        auto &ur = FaceExtra::arr(_ur);
        auto &up = FaceExtra::arr(_up);

        double dx = dx_c + dx_n;
        double delta = dx_n / dx;
        double eta = (12.0 * delta * delta - 1.0) / 24.0;

        for (int i = 0; i < N; ++i) {
            double dzm = u_n[i] - u_c[i];
            double dzp = 2.0 * du_n[i] * dx - dzm;
            double s = limiters::vanAlbada2(dzm, dzp);

            up[i] = ur[i] - s * ((0.5 * delta + s * eta) * dzm + (0.5 * delta - s * eta) * dzp);
        }

        return _up;
    }

private:
    /// Параметры экстраполяции
    double dx_c, dx_n;
    Array u_c, u_n;
    Array du_c, du_n;
};


template<class T>
class FaceExtraTriad {
public:
    static const int N = FaceExtra::SizeOf<T>();
    using Array = std::array<double, N>;

    /// @brief Создание экстраполятора типа ATvL
    /// @param uc Значение в данной ячейке
    /// @param ucx, ucy, ucz Производные в данной ячейке
    /// @param un Значение в смежной ячейке
    /// @param unx, uny, unz Производные в смежной ячейке
    /// @param cell_c Центр данной ячейки
    /// @param neib_c Центр смежной ячейки
    /// @param face_c Центр грани
    FaceExtraTriad(const Array &uc, const Array &ucx, const Array &ucy, const Array &ucz,
                  const Array &un, const Array &unx, const Array &uny, const Array &unz,
                  const Vector3d &cell_c, const Vector3d &neib_c,
                  const Vector3d &face_c) {

        Vector3d rc = face_c - cell_c;
        Vector3d rn = neib_c - face_c;

        dx_c = rc.norm();
        dx_n = rn.norm();

        for (int i = 0; i < N; ++i) {
            u_c[i] = uc[i];
            u_n[i] = un[i];
            du_c[i] = (ucx[i] * rc.x() + ucy[i] * rc.y() + ucz[i] * rc.z()) / dx_c;
            du_n[i] = (unx[i] * rn.x() + uny[i] * rn.y() + unz[i] * rn.z()) / dx_n;
        }

    }

    /// @brief Экстраполяция величины ul на внутреннюю сторону грани
    T m(const T &_ul) const {
        T _um;

        auto &ul = FaceExtra::arr(_ul);
        auto &um = FaceExtra::arr(_um);
        double dx = dx_c + dx_n;
        double eps = 10e-6;

        for (int i = 0; i < N; ++i) {
            double dzp = u_n[i] - u_c[i]; // q_i+1/2
            double dzm = 2.0 * du_c[i] * dx - dzp; // q_i-1/2

            double i2 = 0.5 * std::abs(sgn(dzm) + sgn(dzp));
            double i4 = 0.25 * std::abs(sgn(dzm) + sgn(dzp) + sgn(du_c[i]) + sgn(du_n[i]));

            double kj = i4 / 3 + (1 - i4) * i2;
            double wj = (3 - kj) / (1 - kj + eps) * i2 + (1 - i2); 

            double A = dzp * sgn(dzm);
            double B = dzm * sgn(dzp);

            double dltp = std::max(0.0, std::min(A, wj * B)) * sgn(dzp);
            double dltm = std::max(0.0, std::min(B, wj * A)) * sgn(dzm);

            um[i] = ul[i] + 0.25 * ((1 - kj) * dltm + (1 + kj) * dltp);
        }

        return _um;
    }

    /// @brief Экстраполяция величины ur на внешнюю сторону грани
    T p(const T &_ur) const {
        T _up;

        auto &ur = FaceExtra::arr(_ur);
        auto &up = FaceExtra::arr(_up);
        double dx = dx_c + dx_n;
        double eps = 10e-6;

        for (int i = 0; i < N; ++i) {
            double dzm = u_n[i] - u_c[i];
            double dzp = 2.0 * du_n[i] * dx - dzm;

            double i2 = 0.5 * std::abs(sgn(dzm) + sgn(dzp));
            double i4 = 0.25 * std::abs(sgn(dzm) + sgn(dzp) + sgn(du_c[i]) + sgn(du_n[i]));

            double kj = i4 / 3 + (1 - i4) * i2;
            double wj = (3 - kj) / (1 - kj + eps) * i2 + (1 - i2); 

            double A = dzp * sgn(dzm);
            double B = dzm * sgn(dzp);

            double dltp = std::max(0.0, std::min(A, wj * B)) * sgn(dzm);
            double dltm = std::max(0.0, std::min(B, wj * A)) * sgn(dzp);

            up[i] = ur[i] - 0.25 * ((1 - kj) * dltp + (1 + kj) * dltm);
        }

        return _up;
    }

private:
    /// Параметры экстраполяции
    double dx_c, dx_n;
    Array u_c, u_n;
    Array du_c, du_n;
};

} // math
// zephyr