#pragma once

namespace zephyr::phys {

/// @brief Одномерный тест на мелкую воду
class SwTest1D {
public:
	/// @brief Destructor
	virtual ~SwTest1D() = default;

	/// @brief Глубина (толщина слоя воды)
	virtual double depth(double x, double t) const = 0;

	/// @brief Усредненная скорость
	virtual double speed(double x, double t) const = 0;

	/// @brief Уровень поверхности
	double level(double x, double t) const {
		return bed(x, t) + depth(x, t);
	}

	/// @brief Уровень дна
	virtual double bed(double x) const { return 0.0; }

	/// @brief Уровень дна
	virtual double bed(double x, double t) const { return bed(x); }
};

} // namespace zephyr::phys
