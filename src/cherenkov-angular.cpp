#include <cherenkov-angular.h>
#include <constants.h>

#include <cmath>

// Helper functions (prototypes) - only available within this file

namespace {
	inline
	void
	UpdateParameters(
		const double showerAge,
		const double refractiveIndex,
		const double showerEnergyTeV,
		const ParticleType part);
	
	inline
	double
	FunctionK(
		const double theta);

	inline
	double
	SmallAngleIntegral(
		const double lowAngle,
		const double highAngle);

	inline
	double
	LargeAngleIntegral(
		const double lowAngle,
		const double highAngle);

	inline
	double
	UnnormalizedIntegral(
		const double lowAngle,
		const double highAngle);
		
	inline
	double
	UnnormalizedDensity(
		const double theta);
	
	// Parametrized parameters
	double _nu(0), _w1(0), _w2(0), _ep(0);
	
	// Input variables
	double _showerAge(-1), _refractiveIndex(-1), _showerEnergyTeV(-1);
	ParticleType _primaryParticle(ParticleType::Proton);
	
	// Compound quantities
	double _thetaEm(0), _normalization(0);
	double _cosThetaEm(0), _sinThetaEm(0);
	double _r1(0), _r2(0), _powThetaEmNu(0);
	double _seriesCoefficient[4] = {0,0,0,0};
}



//
// Functions declared in cherenkov-angular.h
//
double
Cherenkov::Angular::PDF(
	const double theta,
	const double showerAge,
	const double refractiveIndex,
	const double showerEnergyTeV,
	const ParticleType primaryParticle)
{
  // Check if theta is outside the range within which the function is defined
  if (theta <= 0 || theta >= Constants::MaxTheta) {
    return 0;
  }
  
  // Compute the parameters of the parametrized function
  UpdateParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  return UnnormalizedDensity(theta) / _normalization;
}

double
Cherenkov::Angular::CDF(
	const double theta,
	const double showerAge,
	const double refractiveIndex,
	const double showerEnergyTeV,
	const ParticleType primaryParticle)
{
  // Check limits of theta
  if (theta <= 0) {
    return 0;
  } else if (theta >= Constants::MaxTheta) {
    return 1;
  }

  // Compute the parameters of the parametrized function
  UpdateParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  return UnnormalizedIntegral(0.0, theta) / _normalization;
}

double
Cherenkov::Angular::Integral(
	const double lowAngle,
  const double highAngle,
  const double showerAge,
  const double refractiveIndex,
  const double showerEnergyTeV,
  const ParticleType primaryParticle)
{
	// Check input angles
	if (lowAngle > highAngle) {
		return Cherenkov::Angular::Integral(highAngle, lowAngle, showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);
	} else if (highAngle <= 0 || lowAngle >= Constants::MaxTheta) {
		return 0;
	} 

  // Compute the parameters of our parametrized function
  UpdateParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  return UnnormalizedIntegral(std::max(lowAngle,0.0), std::min(highAngle,Constants::MaxTheta)) / _normalization;
}



//
// Definition of the helper functions
//
namespace {
	void
	UpdateParameters(
		const double showerAge,
    const double refractiveIndex,
    const double showerEnergyTeV,
    const ParticleType primaryParticle)
	{
		// This fuction computes a modified version of the original parameter set:
		//
		// nu (same as in the paper)
		// w1 = 1/t1 (inverse of theta_1)
		// w2 = (t2-t1)/(t1*t2) (also modified)
		// ep (same as in the paper)
		//
		// Values of showerAge and refractive index are truncated to the following intervals:
		//
		// 0.5 <= showerAge <= 1.5
		// 3.0 x 10^(-5) <= refractiveIndex - 1
		//
		// With this, it is ensured that w1 and w2 are not too large (so exponentials do not result
		// in underflow) and also that w1,w2 > 0.
		//
		// In the computation of thetaEm and its sines, no restriction is put on refractiveIndex
		
		// Check if the input varibles have changed. If not, there is no need to update the parameters
		if (showerAge == _showerAge
		&& refractiveIndex == _refractiveIndex
		&& showerEnergyTeV == _showerEnergyTeV
		&& primaryParticle == _primaryParticle
		) {
			return;
		} 
		
		// Store new value of the input variables
		_showerAge = showerAge;
		_refractiveIndex = refractiveIndex;
		_showerEnergyTeV = showerEnergyTeV;
		_primaryParticle = primaryParticle;

		// Compute auxiliar quantities
		double delta = std::max(_refractiveIndex - 1.0, 3e-5); // bounded to 1e-5!
		double s = std::max(0.5, std::min(_showerAge, 1.5));
		double logs = std::log(s);
		
		// Compute new values for the parametrized parameters
		switch(primaryParticle) {
			// Parametrization for a primary proton
			case ParticleType::Proton:
				_nu = 0.21155 * std::pow(delta, -0.16639) + 1.21803 * logs;
				_w1 = 1.0 / (4.513 * std::pow(delta, 0.45092) * std::pow(showerEnergyTeV, -0.008843) - 0.058687* logs);
				_w2 = _w1 * (1.0 -  1.0 / (0.90725 + 0.41722 * s));
				_ep = 0.009528 + 0.022552 * std::pow(showerEnergyTeV, -0.4207);
				break;
			// Parametrization for a primary gamma		
			case ParticleType::Gamma:
				_nu = 0.34329 * std::pow(delta, -0.10683) + 1.46852 * logs;
				_w1 = 1.0 / (1.4053 * std::pow(delta, 0.32382) - 0.048841 * logs);
				_w2 = _w1 * (1.0 - 1.0 / (0.95734 + 0.26472 * s));
				_ep = 0.0031206;
				break;
		}

		// (global) Quantities depending only on the refractive index
		_thetaEm = std::acos(1.0/_refractiveIndex);
		_cosThetaEm = 1.0/_refractiveIndex;
		_sinThetaEm = std::sqrt((_refractiveIndex-1.0)*(1.0+_refractiveIndex))/_refractiveIndex;
		_seriesCoefficient[0] = _cosThetaEm;
		_seriesCoefficient[1] = _sinThetaEm;
		_seriesCoefficient[2] = - _cosThetaEm;
		_seriesCoefficient[3] = - _sinThetaEm;
		
		// (global) Compound quantities depending on the parametrized parameters
		_powThetaEmNu = std::pow(_thetaEm, _nu);
		_r1 = _thetaEm * _w1;
		_r2 = _thetaEm * _w2;
		_normalization = UnnormalizedIntegral(0.0, Constants::MaxTheta);
	}

	double
	FunctionK(
		const double theta)
	{
		// Note this function requires theta > 0 for nu < 1, but there is not checked here!
		return std::pow(theta, _nu-1.0) * (std::exp(-theta * _w1) + _ep * std::exp(-theta * _w2));
	}

	double
	SmallAngleIntegral(
		const double lowAngle,
		const double highAngle)
	{
		// We do not check values of the input parameters here! It is required that:
		// 0 <= lowAngle,highAngle <= thetaEm <= pi
		
		// TODO: make this value adjustable
		static const int numberOfTerms = 10;

		// Auxiliar quantities
		double deltaLow = _thetaEm - lowAngle;
		double deltaHigh = _thetaEm - highAngle;
		
		double integral = 0;

		// Avoid the computation of 0*log(0)
		if (deltaLow > 1e-20) {
		  integral += (std::cos(lowAngle) - _cosThetaEm) 
		    * (Constants::Pi - std::log(deltaLow/_thetaEm));
		} 

		// Avoid the computation of 0*log(0)
		if (deltaHigh > 1e-20) {
		  integral -= (std::cos(highAngle) - _cosThetaEm)
		    * (Constants::Pi - std::log(deltaHigh/_thetaEm));
		}
		
		// Compute Taylor series using Horner's method
		double seriesLow = _seriesCoefficient[numberOfTerms%4]/numberOfTerms;
		double seriesHigh = _seriesCoefficient[numberOfTerms%4]/numberOfTerms;
		for (int k = numberOfTerms-1; k >= 1; k--) {
			seriesLow = _seriesCoefficient[k%4]/k + seriesLow * deltaLow / double(k+1);
			seriesHigh = _seriesCoefficient[k%4]/k + seriesHigh * deltaHigh / double(k+1);
		}
		integral += seriesLow*deltaLow - seriesHigh*deltaHigh;
	
		// Finish the computation
		integral *= FunctionK(_thetaEm) / _sinThetaEm;

		return integral;
	}



	double
	LargeAngleIntegral(
		const double lowAngle,
		const double highAngle)
	{
		// This is a default step size we use during integration. In practice the value will be close to
		// this one but somewhat smaller. Note that errors from the composite Simpson's rule are 
		// proportional to stepSize^4
		// TODO: make this value adjustable
		static const double defaultStepSize = 5e-4;
		
		const double tLeft = std::sqrt(1 - _thetaEm/lowAngle);
		const double tRight = std::sqrt(1 - _thetaEm/highAngle);
		
		int numberOfSteps = 1 + (tRight - tLeft) / defaultStepSize;
		double stepSize = (tRight - tLeft) / double(numberOfSteps);
		
		auto integrand = [](const double & t)->double{
			double v = (1.0-t)*(1.0+t);
			double value = std::exp(-_r1/v) + _ep * std::exp(-_r2/v);
			value *= t <= 1e-10 ?
				Constants::HalfPi * t + 1.0 - std::pow(t, t)
				: t * (Constants::HalfPi - std::log(t));
			value *= std::pow(v, -(_nu + 1));
			return value;
		};
		
		double integral = 0;
		for (int iStep = 0; iStep < numberOfSteps; iStep++) {
			integral +=	integrand(tLeft + iStep * stepSize)	+ 2.0 * integrand(tLeft + (iStep + 0.5)*stepSize);
		}
		integral = ( 2.0 * integral + integrand(tRight) - integrand(tLeft) ) * stepSize / 6.0;
		integral *= 4.0 * _powThetaEmNu;
		
		return integral;
	}

	double
	UnnormalizedIntegral(
		const double lowAngle,
		const double highAngle)
	{
		// This function requires lowAngle and highAngle to be consistent with 0 <= lowAngle < highAngle
		if (highAngle <= _thetaEm) {
		  // integration interval is contained in [0,thetaEm]
		  return SmallAngleIntegral(lowAngle, highAngle);
		} else if (lowAngle >= _thetaEm) {
		  // integration interval is contained in [thetaEm, infinity)
		  return LargeAngleIntegral(lowAngle, highAngle);
		} else {
		  // integration interval is contained in [0, "infinity") with thetaEm in (lowAngle,highAngle)
		  // so we split the computation in two parts, each correspoding to one of the cases above
		  return SmallAngleIntegral(lowAngle, _thetaEm) + LargeAngleIntegral(_thetaEm, highAngle);
		}
	}

	double
	UnnormalizedDensity(
		const double theta)
	{
		static const double deltaTheta = 1e-5;
		
		if (theta < _thetaEm - deltaTheta) {
	 		// when theta is smaller than thetaEm and no too close from it, compute the density as is
		  return ( std::sin(theta) / _sinThetaEm )
		    * (Constants::Pi - std::log(1.0 - theta/_thetaEm) )
		    * FunctionK(_thetaEm);
		} else if (theta < _thetaEm + deltaTheta) {
			// when theta is too close from thetaEm, return the average value so to avoid a divergence
		  return UnnormalizedIntegral(_thetaEm - deltaTheta, _thetaEm + deltaTheta) * 0.5 / deltaTheta;
		} else {
			// for theta larger than thetaEm and not too close from it, compute the density as defined
		  return (Constants::Pi - std::log(1.0 - _thetaEm/theta)) * FunctionK(theta);
		}
	}

}
