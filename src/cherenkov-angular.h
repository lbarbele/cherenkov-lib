/**
 * @file
 * @brief Defines in the namespace Cherenkov::Angular functions used to compute the angular distribution of 
 * Cherenkov light in air showers.
 */

#pragma once

#include <cmath>

/**
 * @enum mapper::ParticleType
 * @brief Enumeration of primary particle types.
 */
enum class ParticleType {
	Proton,
	Gamma
};

/**
 * @brief Namespace Cherenkov
 */
namespace Cherenkov {
	/** @brief Functions to calculate the angular distribution of %Cherenkov photons in air showers.
	 *
	 * This namespace contains functions implementing the parametrization of the 
	 * angular distribution of Cherenkov photons published in 
	 * [https://doi.org/10.1140/epjc/s10052-021-08971-7](Eur. Phys. J. C 2021 81: 195).
	 *
	 */
	namespace Angular {
	
		/** @brief Evaluates the angular distribution of Cherenkov photons emitted by air showers at a given angle.
		  * 
		  * @param theta Angle measured with respect to the shower axis
		  * @param showerAge Age \fs\f of the shower at the emission altitude
		  * @param refractiveIndex Refractive index at the emission altitude
		  * @param showerEnergyTeV Energy of the primary particle in TeV
		  * @param primaryParticle Type of particle that initiated the shower
		  */
		double
		PDF(
			const double theta,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);

		/** @brief Computes the cumulative distribution of Cherenkov photons with angles below a given value.
		 * 
		 * @param theta Angle measured with respect to the shower axis
		 * @param showerAge Age \fs\f of the shower at the emission altitude
		 * @param refractiveIndex Refractive index at the emission altitude
		 * @param showerEnergyTeV Energy of the primary particle in TeV
		 * @param primaryParticle Type of particle that initiated the shower
		 */  
		double
		CDF(
			const double theta,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);


		/** @brief Computes the angular distribution of Cherenkov photons emitted by air showers.
		 * 
		 * @param lowAngle Lower value of the integration interval
		 * @param highAngle Upper value of the integration interval
		 * @param showerAge Age \fs\f of the shower at the emission altitude
		 * @param refractiveIndex Refractive index at the emission altitude
		 * @param showerEnergyTeV Energy of the primary particle in TeV
		 * @param primaryParticle Type of particle that initiated the shower
		 */
		double
		Integral(
			const double lowAngle,
			const double highAngle,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);
	}
}
