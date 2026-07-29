#pragma once
#include <zephyr/phys/tests/swe/sw_test.h>

namespace zephyr::phys::swe {

/// @class Step
/// @brief Computes dam break with a step solutions
/// @details Class that computes the solutions for a dam break with a step without friction, see \cite Han14.
class Step : public SwTest {
public:
	/// @brief Constructor
	Step();

	/// @brief Destructor
	~Step() override = default;

	/// @brief Название теста
	std::string name() const override { return "Step"; }

	/// @brief Левая граница области
	double x_min() const override { return -1.0; }

	/// @brief Правая граница области
	double x_max() const override { return +1.0; }

	/// @brief Конечное время
	double max_time() const override { return 0.1; }

	/// @brief Поверхность дна
	IBed::Ptr topography() const override;

	/// @brief Уровень дна
	double bed(double x) const override;

	/// @brief Глубина (толщина слоя воды)
	double depth(double x, double t) const override;

	/// @brief Усредненная скорость
	double speed(double x, double t) const override;

private:
	double xdam, step_size; //the dam location and size of the step
	double h_left, h_right; // hl (resp. hr) the water heights on the left (resp. right) of the dam
	double u_left, u_right;
	double h_1, h_2; //water heights before, at and after the step respectively
	double u_1, u_2; //speed before, at and after the step respectively
};

} // zephyr::phys::swe