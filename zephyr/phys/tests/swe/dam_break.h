#pragma once
#include <zephyr/phys/tests/swe/sw_test.h>

namespace zephyr::phys::swe {

/// @details Computes the solutions for a dam break without friction
/// [1] A. Ritter. Die Fortpflanzung der Wasserwellen. Zeitschrift des Vereines
/// Deuscher Ingenieure, 36(33):947–954, 1892. 23, 24
/// [2] J. J. Stoker. Water Waves: The Mathematical Theory with Applications.
/// Pure and Applied Mathematics. Interscience Publishers, New York, USA, 1957. 23, 24.
class DamBreak : public SwTest {
public:
	/// @brief Тесты по умолчанию: 1 - обычный распад разрыва,
	/// 2 - распад разрыва с сухим дном.
	explicit DamBreak(int sol);

	/// @brief Constructor
	DamBreak(double depth_L, double depth_R, double x_dam = 0.0);

	/// @brief Destructor
	~DamBreak() override = default;

	/// @brief Название теста
	std::string name() const override { return "DamBreak"; }

	/// @brief Левая граница области
	double x_min() const override { return -1.0; }

	/// @brief Правая граница области
	double x_max() const override { return +1.0; }

	/// @brief Конечное время
	double max_time() const override { return 0.2; }

	/// @brief Глубина (толщина слоя воды)
	double depth(double x, double t) const override;

	/// @brief Усредненная скорость
	double speed(double x, double t) const override;

private:
	/// @brief Computes the solution
	void compute();

	double h_L, h_R; // hl (resp. hr) the water heights on the left resp. right) of the dam
	double x0; //the dam location

	double h_mid, u_mid, c_mid; //water height, velocity and wave velocity in the intermediate state (only for the dam break on wet soil: Stoker's solution)
	double v; //the shock velocity (only for the Stoker's solution)
	double c_L, c_R; // cl (resp. cr) left (resp. right) wave velocity
};

} // namespace zephyr::phys
