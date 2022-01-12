#pragma once

/**
 * @file
 * @brief Definition of mathematical and physical constants.
 */

/**
 * @brief Namespace containing the definition of useful mathematical and physical constants.
 */
namespace Constants {
	inline constexpr double Pi = 3.14159265358979311599796346854; ///< Pi in double precision
	inline constexpr double HalfPi = 0.5 * Pi; ///< Half of pi
	inline constexpr double PiPlusOne = Pi + 1.0; ///< Pi plus 1
	inline constexpr double AlphaEm = 0.0072973525693; ///< Fine-structure constant
	inline constexpr double ElectronMass = 0.51099895; ///< Electron mass in MeV
}
