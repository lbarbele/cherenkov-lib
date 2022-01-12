/**
 * @file
 * @brief Defines in the namespace Cherenkov::Angular functions used to compute the angular distribution of 
 * Cherenkov light in air showers.
 */

#pragma once

#include <particle-type.h>

/**
 * @brief Namespace Cherenkov
 */
namespace Cherenkov {
	
	class AngularDistribution {
	private:
		static const double fPeakWidth;
		static const double fMinRefractiveIndex;
		static const double fMinShowerAge;
		static const double fMaxShowerAge;
		static const double fMaxTheta;
		static const int fSeriesTruncationOrder;
		
		double fEnergyTeV;
		ParticleType fParticleType;
		
		double fNu, fOmega1, fOmega2, fEpsilon, fNormalization;
		double fThetaEm, fCosThetaEm, fSinThetaEm;
		double fTwoThetaEm, fCosTwoThetaEm, fSinTwoThetaEm;
	
		void UpdateParameters(const double showerAge, const double refractiveIndex);
		double FunctionK(const double theta) const;
		double UnnormalizedIntegral(const double lowAngle, const double highAngle) const;
	
	public:
		AngularDistribution(const double showerEnergyTeV, const ParticleType particleType);
		double PDF(const double theta, const double showerAge, const double refractiveIndex);
		double CDF(const double theta, const double showerAge, const double refractiveIndex);
		double Integral(const double lowAngle, const double highAngle, const double showerAge, const double refractiveIndex);
		double IntegralWithSine(const double lowAngle, const double highAngle, const double showerAge, const double refractiveIndex);
	};
}
