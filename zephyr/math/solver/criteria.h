#pragma once
#include <cmath>
#include <variant>

namespace zephyr::math {

/// @brief Критерий адаптации по перепаду величин
struct SlopeCriterion {
	static constexpr double std_ratio = 0.4;
	static constexpr double inf = std::numeric_limits<double>::max();

	double dens_split = 0.05;  ///< Адаптируется, если хотя бы один больше
	double dens_merge = 0.02;  ///< Огрубляется, если все перепады меньше

	double pres_split = 0.05;
	double pres_merge = 0.02;

	/// @brief Стандартная адаптация по плотности и давлению
	static constexpr SlopeCriterion Default() {
		return {};
	}

	/// @brief Адаптация только по плотности
	static SlopeCriterion Density(double split_threshold) {
		return {
			.dens_split = split_threshold,
			.dens_merge = std_ratio * split_threshold,
			.pres_split = inf,
			.pres_merge = inf
		};
	}

	/// @brief Адаптация только по плотности
	static SlopeCriterion Density(double split_threshold, double merge_threshold) {
		return {
			.dens_split = split_threshold,
			.dens_merge = std::min(merge_threshold, 0.5 * split_threshold),
			.pres_split = inf,
			.pres_merge = inf
		};
	}

	/// @brief Адаптация только по давлению
	static SlopeCriterion Pressure(double split_threshold) {
		return {
			.dens_split = inf,
			.dens_merge = inf,
			.pres_split = split_threshold,
			.pres_merge = std_ratio * split_threshold
		};
	}

	/// @brief Адаптация только по давлению
	static SlopeCriterion Pressure(double split_threshold, double merge_threshold) {
		return {
			.dens_split = inf,
			.dens_merge = inf,
			.pres_split = split_threshold,
			.pres_merge = std::min(merge_threshold, 0.5 * split_threshold),
		};
	}

	/// @brief Адаптация по плотности и давлению, стандартный порог огрубления
	static SlopeCriterion DensPres(double dens_split, double pres_split) {
		return {
			.dens_split = dens_split,
			.dens_merge = std_ratio * dens_split,
			.pres_split = pres_split,
			.pres_merge = std_ratio * pres_split
		};
	}
};

/// @brief Критерий адаптации по перепаду величин
struct ChiCriterion {
	static constexpr double std_ratio = 0.4; // ??
	static constexpr double inf = std::numeric_limits<double>::max();

	double epsilon = 0.01;    ///< Порог сглаживания осцилляций

	double dens_split = 0.15;  ///< Адаптируется, если хотя chi больше
	double dens_merge = 0.15;  ///< Огрубляется, если chi меньше

	double pres_split = 0.1;  // ??
	double pres_merge = 0.1;  // ??

	/// @brief Стандартная адаптация по плотности и давлению
	static ChiCriterion Default() {
		return {};
	}

	/// @brief Адаптация только по плотности
	static ChiCriterion Density(double split_threshold) {
		return {
			.dens_split = split_threshold,
			.dens_merge = std_ratio * split_threshold,
			.pres_split = inf,
			.pres_merge = inf
		};
	}

	/// @brief Адаптация только по плотности
	static ChiCriterion Density(double split_threshold, double merge_threshold) {
		return {
			.dens_split = split_threshold,
			.dens_merge = std::min(merge_threshold, 0.5 * split_threshold),
			.pres_split = inf,
			.pres_merge = inf
		};
	}

	/// @brief Адаптация только по давлению
	static ChiCriterion Pressure(double split_threshold) {
		return {
			.dens_split = inf,
			.dens_merge = inf,
			.pres_split = split_threshold,
			.pres_merge = std_ratio * split_threshold
		};
	}

	/// @brief Адаптация только по давлению
	static ChiCriterion Pressure(double split_threshold, double merge_threshold) {
		return {
			.dens_split = inf,
			.dens_merge = inf,
			.pres_split = split_threshold,
			.pres_merge = std::min(merge_threshold, 0.5 * split_threshold),
		};
	}

	/// @brief Адаптация по плотности и давлению, стандартный порог огрубления
	static ChiCriterion DensPres(double dens_split, double pres_split) {
		return {
			.dens_split = dens_split,
			.dens_merge = std_ratio * dens_split,
			.pres_split = pres_split,
			.pres_merge = std_ratio * pres_split
		};
	}
};

/// @brief Один из критериев адаптации
using Criterion = std::variant<SlopeCriterion, ChiCriterion>;

} // namespace zephyr::math