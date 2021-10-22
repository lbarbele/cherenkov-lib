/**
 * @file
 * @brief Defines in the namespace Cherenkov::Angular methods used to evaluate the longitudinal profiles of Cherenkov light emission in air showers.
 */

#pragma once

namespace Cherenkov {
	/** @brief Functions to calculate the longitudinal profiles of Cherenkov-light emission.
	 *
	 */
	namespace Longi {
		/** @brief Evaluates the average number of Cherenkov photons per electron per cm.
		 * 
		 * @param showerAge Age \fs\f of the shower at the emission altitude
		 * @param cutoffEnergy Cutoff energy (minimum value) considered for the shower energy spectrum of electrons [MeV]
		 * @param refractiveIndex Refractive index at the emission altitude
		 * @param lambdaMin Lower value of the wavelenght interval of Cherenkov light [cm]
		 * @param lambdaMax Upper value of the wavelenght interval of Cherenkov light [cm]
		 *
		 */
		double
		YieldFactor(
			const double showerAge, 
			const double cutoffEnergy, 
			const double refractiveIndex,
			const double lambdaMin,
			const double lambdaMax);
	}
}
