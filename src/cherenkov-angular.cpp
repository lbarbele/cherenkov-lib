#include <cherenkov-angular.h>
#include <constants.h>
#include <math.h>

#include <cmath>
#include <array>
#include <iostream>


const double Cherenkov::AngularDistribution::fPeakWidth = 1e-4;
const double Cherenkov::AngularDistribution::fMinRefractiveIndex = 1.0 + 3e-5;
const double Cherenkov::AngularDistribution::fMinShowerAge = 0.5;
const double Cherenkov::AngularDistribution::fMaxShowerAge = 1.5;
const double Cherenkov::AngularDistribution::fMaxTheta = Constants::Pi/3.0;
const int Cherenkov::AngularDistribution::fSeriesTruncationOrder = 8;



Cherenkov::AngularDistribution::AngularDistribution(
  const double showerEnergyTeV,
  const ParticleType particleType
) : 
  fEnergyTeV(showerEnergyTeV),
  fParticleType(particleType)
{
}



void
Cherenkov::AngularDistribution::UpdateParameters(
	const double showerAge,
	const double refractiveIndex
)
{
	const double nMinusOne = std::max(refractiveIndex, fMinRefractiveIndex) - 1.0;
	const double truncatedAge = std::max(fMinShowerAge, std::min(showerAge, fMaxShowerAge));
	const double logShowerAge = std::log(truncatedAge);
	
	// Compute the parametrized parameters
	switch(fParticleType) {
	case ParticleType::Proton:
		fNu = 0.21155 * std::pow(nMinusOne, -0.16639) + 1.21803 * logShowerAge;
		fOmega1 = 1.0 / (4.513 * std::pow(nMinusOne, 0.45092) * std::pow(fEnergyTeV, -0.008843) -
			0.058687 * logShowerAge);
		fOmega2 = fOmega1 * (1.0 -  1.0 / (0.90725 + 0.41722 * truncatedAge));
		fEpsilon = 0.009528 + 0.022552 * std::pow(fEnergyTeV, -0.4207);
		break;
	case ParticleType::Gamma:
		fNu = 0.34329 * std::pow(nMinusOne, -0.10683) + 1.46852 * logShowerAge;
		fOmega1 = 1.0 / (1.4053 * std::pow(nMinusOne, 0.32382) - 0.048841 * logShowerAge);
		fOmega2 = fOmega1 * (1.0 - 1.0 / (0.95734 + 0.26472 * truncatedAge));
		fEpsilon = 0.0031206;
		break;
	}
	
	// The maximum Cherenkov emission angle (depends on the refractive index only)
	fThetaEm = std::acos(1.0/refractiveIndex);
	fCosThetaEm = 1.0/refractiveIndex;
	fSinThetaEm = std::sqrt((refractiveIndex-1.0)*(refractiveIndex+1.0))/refractiveIndex;
	fTwoThetaEm = 2.0*fThetaEm;
	fCosTwoThetaEm = (2.0*fCosThetaEm - refractiveIndex) * fCosThetaEm;
	fSinTwoThetaEm = 2.0*fSinThetaEm/refractiveIndex;
	
	// The normalization factor (depends on parameters and thetaEm)
	fNormalization = UnnormalizedIntegral(0.0, fMaxTheta);
}



double
Cherenkov::AngularDistribution::FunctionK(const double theta)
const
{
	return std::pow(theta,fNu-1.0) *
		(std::exp(-theta*fOmega1) + fEpsilon*std::exp(-theta*fOmega2));
}



double
Cherenkov::AngularDistribution::UnnormalizedIntegral(
	const double lowAngle,
	const double highAngle
) const
{
	// This variable holds the value that will be returned
	double integral = 0;

	// if lowAngle is to the left of the peak, integrate between lowAngle and min(thetaEm, highAngle)
	if (lowAngle < fThetaEm) {
	
		// compute delta = thetaEm - theta, ensuring highAngle <= theta
		const double deltaThetaLow = fThetaEm - lowAngle;
		const double deltaThetaHigh = std::max(fThetaEm - highAngle, 0.0);
		
		double integralLeft = 0;
		
		// Avoid the computation of x*log(x) at x = 0, which approaches 0 as x -> 0
		if (deltaThetaLow > 1e-20) {
			integralLeft += (std::cos(lowAngle) - fCosThetaEm) *
				(Constants::Pi - std::log(deltaThetaLow/fThetaEm));
		} 

		if (deltaThetaHigh > 1e-20) {
			integralLeft -= (std::cos(fThetaEm-deltaThetaHigh) - fCosThetaEm) * 
				(Constants::Pi - std::log(deltaThetaHigh/fThetaEm));
		}
		
		// The integral is computed as a truncated Taylor series, below are the coefficients
		const double seriesFactor[4] = {fCosThetaEm, fSinThetaEm, -fCosThetaEm, -fSinThetaEm};
		
		double seriesLow = 0; /* related to lowAngle */
		double seriesHigh = 0; /* related to highAngle */
		
		// Apply Horner's rule to compute the series
		for (int k = fSeriesTruncationOrder; k >= 1; k--) {
			seriesLow = seriesFactor[k%4]/double(k) + seriesLow*deltaThetaLow/double(k+1);
			seriesHigh = seriesFactor[k%4]/double(k) + seriesHigh*deltaThetaHigh/double(k+1);
		}
		
		// Add the series to the integral and apply the multiplication factor
		integralLeft = (integralLeft + seriesLow*deltaThetaLow - seriesHigh*deltaThetaHigh) *
			FunctionK(fThetaEm)/fSinThetaEm;
			
		// Add to the overall integral value
		integral += integralLeft;
	}
	
	// If highAngle is to the right of the peak, integrate between max(thetaEm,lowAngle) and highAngle
	if (highAngle > fThetaEm) {
	
		// Compute the integration limits after the change of variables
		const double tMin = std::sqrt(std::max(1.0 - fThetaEm/lowAngle, 0.0));
		const double tMax = std::sqrt(1.0 - fThetaEm/highAngle);
		
		// Actually compute the integral
		double integralRight = (4.0 / fThetaEm) * Math::RombergIntegral(tMin, tMax, 1e-10, [=](double t) {
			const double theta = fThetaEm / ((1.0-t)*(1.0+t));
			return std::pow(theta, fNu + 1) *
				(t <= 1e-10 ? t*Constants::HalfPi + 1.0 - std::pow(t, t) : t * (Constants::HalfPi - std::log(t))) *
				(std::exp(-fOmega1*theta) + fEpsilon*std::exp(-fOmega2*theta));
		});
		
		// Add to the overall integral value
		integral += integralRight;
	}
	
	return integral;
}



double
Cherenkov::AngularDistribution::PDF(
	const double theta, 
	const double showerAge,
	const double refractiveIndex
)
{
	// Check value of theta
	if (theta < 0) {
		return 0;
	} else if (theta > fMaxTheta) {
		return 0;
	}

	// Update the parameters depending on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	if (theta <= fThetaEm - fPeakWidth) {
		// Region left of the peak: compute the function as is
		return 	FunctionK(fThetaEm) * (Constants::Pi - std::log(1.0-theta/fThetaEm)) *
			(std::sin(theta)/fSinThetaEm) / fNormalization;
	}
	else if (theta < fThetaEm + fPeakWidth) {
		// Very close to the peak: compute the average in an interval containing the peak
		return UnnormalizedIntegral(fThetaEm - fPeakWidth, fThetaEm + fPeakWidth) /
			(2.0*fPeakWidth*fNormalization);
	}
	else {
		// Right of the peak: compute the function as is
		return FunctionK(theta) * (Constants::Pi - std::log(1.0 - fThetaEm/theta)) / fNormalization;
	}
}



double
Cherenkov::AngularDistribution::CDF(
	const double theta, 
	const double showerAge,
	const double refractiveIndex
)
{
	// Check value of theta
	if (theta < 0) {
		return 0;
	} else if (theta > fMaxTheta) {
		return 1;
	}
	
	// Update the parameters depending on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	// Return the (normalized) integral of the distribution between 0 and theta (given)
	return UnnormalizedIntegral(0, theta) / fNormalization;
}



double
Cherenkov::AngularDistribution::Integral(
	const double lowAngle,
	const double highAngle, 
	const double showerAge,
	const double refractiveIndex
)
{
	// Update the parameters depending on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	// Return the (normalized) integral of the distribution between the given angles
	return UnnormalizedIntegral(std::max(lowAngle,0.0), std::min(highAngle, fMaxTheta)) /
		fNormalization;
}



double
Cherenkov::AngularDistribution::IntegralWithSine(
	const double lowAngle, 
	const double highAngle, 
	const double showerAge,
	const double refractiveIndex
)
{
	// Reorder angles, if necessary
	if (lowAngle > highAngle) {
		return -IntegralWithSine(highAngle, lowAngle, showerAge, refractiveIndex);
	}
	
	// Check value of theta
	if (highAngle > fMaxTheta) {
		return IntegralWithSine(lowAngle, fMaxTheta, showerAge, refractiveIndex);;
	} else if (lowAngle < 0) {
		return IntegralWithSine(0, highAngle, showerAge, refractiveIndex);
	}
	
	// Update the parameters depending on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	// This variable holds the value that will be returned
	double integral = 0;

	// if lowAngle is to the left of the peak, integrate between lowAngle and min(thetaEm, highAngle)
	if (lowAngle < fThetaEm) {
	
		// compute delta = 2*(thetaEm - theta), ensuring highAngle <= theta
		const double deltaThetaLow = 2.0*(fThetaEm - lowAngle);
		const double deltaThetaHigh = 2.0*std::max(fThetaEm - highAngle, 0.0);
		
		// Initialize the integral
		double integralLeft = deltaThetaLow - deltaThetaHigh;

		// Avoid the computation of x*log(x) at x = 0 (note that x*log(x) -> 0 as x -> 0)
		if (deltaThetaHigh > 1e-20) {
			integralLeft -= (deltaThetaHigh + std::sin(fTwoThetaEm-deltaThetaHigh) - fSinTwoThetaEm) * 
				(Constants::Pi - std::log(deltaThetaHigh/fTwoThetaEm));
		}
		
		if (deltaThetaLow > 1e-20) {
			integralLeft += (deltaThetaLow + std::sin(2.0*lowAngle) - fSinTwoThetaEm) * 
				(Constants::Pi - std::log(deltaThetaLow/fTwoThetaEm));
		}
		
		// The integral is computed as a truncated Taylor series, below are the coefficients
		const double seriesFactor[4] = {fSinTwoThetaEm, -fCosTwoThetaEm, 
			-fSinTwoThetaEm, fCosTwoThetaEm};
		
		double seriesLow = 0; /* related to lowAngle */
		double seriesHigh = 0; /* related to highAngle */
		
		// Compute the series using Horner's rule
		for (int k = fSeriesTruncationOrder; k >= 1; k--) {
			seriesLow = seriesFactor[k%4]/double(k) + seriesLow * deltaThetaLow / double(k+1);
			seriesHigh = seriesFactor[k%4]/double(k) + seriesHigh * deltaThetaHigh / double(k+1);
		}
		
		// Add the series to the integral and apply any necessary multiplication factors
		integralLeft = (integralLeft + seriesLow*deltaThetaLow - seriesHigh*deltaThetaHigh) *
			0.25 * FunctionK(fThetaEm) / fSinThetaEm;
		
		// Add the result to the overall integral
		integral += integralLeft;
	}
	
	// If highAngle is to the right of the peak, integrate between max(thetaEm,lowAngle) and highAngle
	if (highAngle > fThetaEm) {
	
		// Compute the integration limits after the change of variables
		const double tMin = std::sqrt(std::max(1.0 - fThetaEm/lowAngle, 0.0));
		const double tMax = std::sqrt(1.0 - fThetaEm/highAngle);
		
		// Actually compute the integral
		double integralRight = (4.0 / fThetaEm) * Math::RombergIntegral(tMin, tMax, 1e-10, [=](double t) {
			const double theta = fThetaEm / ((1.0-t)*(1.0+t));
			return std::sin(theta) * std::pow(theta, fNu + 1) *
				(t <= 1e-10 ? t*Constants::HalfPi + 1.0 - std::pow(t, t) : t * (Constants::HalfPi - std::log(t))) *
				(std::exp(-fOmega1*theta) + fEpsilon*std::exp(-fOmega2*theta));
		});
	
		// Add the result to the overall integral value
		integral += integralRight;
	}
	
	// Return the integral taking into account the normalization of the distribution function
	return integral / fNormalization;
}


