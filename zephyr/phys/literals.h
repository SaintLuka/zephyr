#pragma once

namespace zephyr::phys {

/// ============================ Единицы длины ================================

constexpr double operator"" _nm(long double x) { return double(x*1.0e-9); }
constexpr double operator"" _um(long double x) { return double(x*1.0e-6); }
constexpr double operator"" _mm(long double x) { return double(x*1.0e-3); }
constexpr double operator"" _cm(long double x) { return double(x*1.0e-2); }
constexpr double operator"" _m (long double x) { return double(x); }
constexpr double operator"" _km(long double x) { return double(x*1.0e+3); }

/// =========================== Единицы объема ================================

constexpr double operator"" _um3(long double x) { return double(x*1.0e-18); }
constexpr double operator"" _mm3(long double x) { return double(x*1.0e-9); }
constexpr double operator"" _cm3(long double x) { return double(x*1.0e-6); }
constexpr double operator"" _m3 (long double x) { return double(x); }

/// ========================== Единицы времени ================================

constexpr double operator"" _ns(long double x) { return double(x*1.0e-9); }
constexpr double operator"" _us(long double x) { return double(x*1.0e-6); }
constexpr double operator"" _ms(long double x) { return double(x*1.0e-3); }
constexpr double operator"" _s (long double x) { return double(x); }

/// =========================== Единицы скорости ==============================

constexpr double operator"" _cm_s(long double x) { return double(x*1.0e-2); }
constexpr double operator"" _m_s (long double x) { return double(x); }
constexpr double operator"" _km_s(long double x) { return double(x*1.0e+3); }

/// ======================== Единицы массы и плотности ========================

constexpr double operator"" _g    (long double x) { return double(x*1.0e-3); }
constexpr double operator"" _kg   (long double x) { return double(x); }
constexpr double operator"" _kg_m3(long double x) { return double(x); }
constexpr double operator"" _g_cm3(long double x) { return double(x*1.0e3); }

/// ============================ Единицы давления =============================

constexpr double operator"" _Pa (long double x) { return double(x); }
constexpr double operator"" _kPa(long double x) { return double(x*1.0e3); }
constexpr double operator"" _bar(long double x) { return double(x*1.0e5); }
constexpr double operator"" _MPa(long double x) { return double(x*1.0e6); }
constexpr double operator"" _GPa(long double x) { return double(x*1.0e9); }
constexpr double operator"" _TPa(long double x) { return static_cast<double>(x*1.0e12); }

/// =========================== Единицы энергии ===============================

constexpr double operator"" _eV (long double x) { return double(x*1.602176634e-19); }
constexpr double operator"" _keV(long double x) { return double(x*1.602176634e-16); }
constexpr double operator"" _MeV(long double x) { return double(x*1.602176634e-13); }
constexpr double operator"" _GeV(long double x) { return double(x*1.602176634e-10); }

constexpr double operator"" _uJ (long double x) { return double(x*1.0e-6); }
constexpr double operator"" _mJ (long double x) { return double(x*1.0e-3); }
constexpr double operator"" _J  (long double x) { return double(x); }
constexpr double operator"" _kJ (long double x) { return double(x*1.0e+3); }
constexpr double operator"" _MJ (long double x) { return double(x*1.0e+6); }
constexpr double operator"" _GJ (long double x) { return double(x*1.0e+9); }

/// ============================= Сложные единицы =============================

/// @brief Удельная энергия
constexpr double operator"" _J_kg (long double x) { return double(x); }
constexpr double operator"" _kJ_kg(long double x) { return double(x*1.0e3); }
constexpr double operator"" _MJ_kg(long double x) { return double(x*1.0e6); }

/// @brief Удельная теплоемкость
constexpr double operator"" _J_kgK(long double x) { return double(x); }

// @brief Теплопроводность
constexpr double operator"" _W_mK(long double x) { return double(x); }

// @brief Градусы Цельсия
constexpr double operator"" _C(long double x) { return double(x + 273.15); }

}