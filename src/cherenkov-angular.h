/**
 * @file
 * @brief Defines in the namespace Cherenkov::Angular functions used to compute the angular distribution of 
 * Cherenkov light in air showers.
 */

#pragma once

#include <constants.h>
#include <particle-type.h>

/**
 * @brief Namespace Cherenkov
 */
namespace Cherenkov {
	
	class AngularDistribution {
	public:
		double fEnergyTeV, fMaxTheta;
		ParticleType fParticleType;
		
		double fShowerAge, fRefractiveIndex;
		
		double fNu, fOmega1, fOmega2, fEpsilon, fNormalization;
		double fThetaEm, fCosThetaEm, fSinThetaEm;
		double fTwoThetaEm, fCosTwoThetaEm, fSinTwoThetaEm;
	
		void UpdateParameters(const double showerAge, const double refractiveIndex);
		double FunctionK(const double theta) const;
		
		double UnnormalizedIntegral(const double lowAngle, const double highAngle) const;
		double UnnormalizedIntegralLeft(const double lowAngle, const double highAngle) const;
		double UnnormalizedIntegralRight(const double lowAngle, const double highAngle) const;
		
		double UnnormalizedIntegralSineLeft(const double lowAngle, const double highAngle) const;
		double UnnormalizedIntegralSineRight(const double lowAngle, const double highAngle) const;
	
	public:
		AngularDistribution(const double showerEnergyTeV, const ParticleType particleType, const double maxTheta = Constants::Pi);
		double PDF(const double theta, const double showerAge, const double refractiveIndex);
		double CDF(const double theta, const double showerAge, const double refractiveIndex);
		double Integral(const double lowAngle, const double highAngle, const double showerAge, const double refractiveIndex);
		double IntegralSine(const double lowAngle, const double highAngle, const double showerAge, const double refractiveIndex);
	};
}
