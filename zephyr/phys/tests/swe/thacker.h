#pragma once
#include <zephyr/phys/tests/swe/sw_test.h>

namespace zephyr::phys::swe {

/// @class Thacker1D
/// @brief Computes %Thacker1D solution
/// @details Class that computes the solution for %Thacker1D parabola, see \cite Thacker81.
class Thacker1D : public SwTest {
public:
	/// @brief Тест по умолчанию
	Thacker1D();

	/// @brief Constructor
	Thacker1D(double coeff, double bottom);

	/// @brief Destructor
	~Thacker1D() override = default;

	/// @brief Название теста
	std::string name() const override { return "Thacker1D"; }

	/// @brief Левая граница области
	double x_min() const override { return -lx; }

	/// @brief Правая граница области
	double x_max() const override { return +lx; }

	/// @brief Конечное время
	double max_time() const override;

	/// @brief Поверхность дна
	IBed::Ptr topography() const override;

	/// @brief Уровень дна
	double bed(double x) const override;

	/// @brief Глубина (толщина слоя воды)
	double depth(double x, double t) const override;

	/// @brief Усредненная скорость
	double speed(double x, double t) const override;

private:
	double coeff, bottom;
	double omega;
	double a, h0;
	double lx;
};

/// @class Thacker2D
/// @brief Computes %Thacker1D solutions in 2D
/// @details Class that computes the solutions for %Thacker1D paraboloid, see \cite Thacker81.
class Thacker2D : public SwTest2D {
public:
	/// @brief Constructor
	/// @param sol 1 - вращение плоскости, 2 - колебание параболы
	explicit Thacker2D(int sol);

	/// @brief Destructor
	~Thacker2D() override = default;

	/// @brief Название теста
	std::string name() const override { return "Thacker2D"; }

	/// @brief Левая граница области
	double x_min() const override { return -1.0; }

	/// @brief Правая граница области
	double x_max() const override { return +1.0; }

	/// @brief Нижняя граница области
	double y_min() const override { return -1.0; }

	/// @brief Верхняя граница области
	double y_max() const override { return +1.0; }

	/// @brief Конечное время
	double max_time() const override;

	/// @brief Поверхность дна
	IBed::Ptr topography() const override;

	/// @brief Уровень дна
	double bed(const Vector3d& v) const override;

	/// @brief Глубина (толщина слоя воды)
	double depth(const Vector3d& v, double t) const override;

	/// @brief Усредненная скорость
	Vector2d speed(const Vector3d& v, double t) const override;

private:
	int sol;
	double coeff, bottom; ///< Характеристики дна
	double a, h0, omega;  ///< Производные параметры
	double eta, r0;		  ///< Характеристики начальных данных

};

} //  zephyr::phys::swe
