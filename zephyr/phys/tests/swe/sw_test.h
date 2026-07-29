#pragma once
#include <zephyr/geom/bed.h>

namespace zephyr::phys::swe {

using geom::IBed;
using geom::Vector2d;
using geom::Vector3d;

inline constexpr double GRAV = 9.81;

/// @brief Тест на мелкую воду
class SwTest {
public:
	/// @brief Destructor
	virtual ~SwTest() = default;

	/// @brief Название теста
	virtual std::string name() const { return "SwTest"; }

	/// @brief Левая граница области
	virtual double x_min() const { return -1.0; }

	/// @brief Правая граница области
	virtual double x_max() const { return +1.0; }

	/// @brief Нижняя граница области
	virtual double y_min() const { return -1.0; }

	/// @brief Верхняя граница области
	virtual double y_max() const { return +1.0; }

	/// @brief Конечный момент времени
	virtual double max_time() const { return std::numeric_limits<double>::max(); }

	/// @brief Поверхность дна
	virtual IBed::Ptr topography() const {
		return geom::ConstBed::create(bed(0.0));
	}

	// ------------------ Одномерные тесты ------------------

	/// @brief Уровень дна
	virtual double bed(double x) const { return 0.0; }

	/// @brief Уровень дна
	virtual double bed(double x, double t) const { return bed(x); }

	/// @brief Глубина (толщина слоя воды)
	virtual double depth(double x, double t) const = 0;

	/// @brief Усредненная скорость
	virtual double speed(double x, double t) const = 0;

	/// @brief Уровень поверхности
	double level(double x, double t) const {
		return bed(x, t) + depth(x, t);
	}

	// ------------------ Двумерные тесты ------------------

	/// @brief Уровень дна
	virtual double bed(const Vector3d& v) const {
		return bed(v.x());
	}

	/// @brief Уровень дна
	virtual double bed(const Vector3d& v, double t) const {
		return bed(v);
	}

	/// @brief Глубина (толщина слоя воды)
	virtual double depth(const Vector3d& v, double t) const {
		return depth(v.x(), t);
	};

	/// @brief Усредненная скорость
	virtual Vector2d speed(const Vector3d& v, double t) const {
		return {speed(v.x(), t), 0.0};
	};

	/// @brief Уровень поверхности
	double level(const Vector3d& v, double t) const {
		return bed(v, t) + depth(v, t);
	}
};

/// @brief Двумерный тест на мелкую воду. Запретили одномерные функции.
class SwTest2D : public SwTest {
public:
	/// @brief Destructor
	~SwTest2D() override = default;

	// ----------- Одномерные тесты ------------------

	/// @brief Уровень дна
	double bed(double x) const override {
		throw std::runtime_error("call bed() 1D version");
	}

	/// @brief Уровень дна
	double bed(double x, double t) const override {
		throw std::runtime_error("call bed() 1D version");
	}

	/// @brief Глубина (толщина слоя воды)
	double depth(double x, double t) const override {
		throw std::runtime_error("call depth() 1D version");
	}

	/// @brief Усредненная скорость
	double speed(double x, double t) const override {
		throw std::runtime_error("call speed() 1D version");
	}
};

} // namespace zephyr::phys::swe
