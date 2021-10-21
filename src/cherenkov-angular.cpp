#include "cerang.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>

std::array<double, 4> 
CherenkovParameters(const double showerAge,
                    const double refractiveIndex,
                    const double showerEnergyTeV,
                    const ParticleType part);
                    
double
FunctionK(const double theta,
          const std::array<double, 4> parameters);
          
double
SmallAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters);
                   
double
LargeAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters);
                   
double
CherenkovUnnormalizedIntegral(const double lowAngle,
                              const double highAngle,
                              const double thetaEm, 
                              const std::array<double,4> parameters);
                              
double
CherenkovUnnormalizedDensity(const double theta,
                             const double thetaEm, 
                             const std::array<double,4> parameters);










std::array<double, 4> 
CherenkovParameters(const double showerAge,
                    const double refractiveIndex,
                    const double showerEnergyTeV,
                    const ParticleType primaryParticle)
{
	// we do not check input parameters here! that is, it is assumed that:
  //
  // 0 < showerAge < 3 (a value of 0 would make log(s) diverge, 3 is unphysical)
  // 1.0 < refractiveIndex (pow() will return NaN if refractiveIndex = 1)
  // showerEnergy > 0 (pow() of showerEnergy = 0 is also NaN)
  //
  // during tests, we ensure that output is correct if the above bounds are respected
  
  // the output of this function is a modified version of the original parameters,
  // changes are made for numerical convenience
  //
  // nu (same as in the paper)
  // w1 = 1/t1 (inverse of theta_1)
  // w2 = (t2-t1)/(t1*t2) (also modified)
  // ep (same as in the paper)

  // Creat the object that will be returned
  std::array<double, 4> arrayOfParameters = {0, 0, 0, 0};

  // Give friendly names to the elements of the array, for readability
  auto & nu = arrayOfParameters[0];
  auto & w1 = arrayOfParameters[1];
  auto & w2 = arrayOfParameters[2];
  auto & ep = arrayOfParameters[3];

  // Compute auxiliar quantities
  double delta = refractiveIndex - 1.0;
  double logAge = std::log(showerAge);
  
  switch(primaryParticle) {
  
  	// Parametrization for a primary proton
  	case ParticleType::Proton:
			nu = 0.21155 * std::pow(delta, -0.16639) + 1.21803 * logAge;
			w1 = 1.0 / (4.513 * std::pow(delta, 0.45092) * std::pow(showerEnergyTeV, -0.008843) - 0.058687* logAge);
			w2 = w1 * (1.0 -  1.0 / (0.90725 + 0.41722 * showerAge));
			ep = 0.009528 + 0.022552 * std::pow(showerEnergyTeV, -0.4207);
			break;
			
		// Parametrization for a primary gamma		
  	case ParticleType::Gamma:
		  nu = 0.34329 * std::pow(delta, -0.10683) + 1.46852 * logAge;
			w1 = 1.0 / (1.4053 * std::pow(delta, 0.32382) - 0.048841 * logAge);
			w2 = w1 * (1.0 - 1.0 / (0.95734 + 0.26472 * showerAge));
			ep = 0.0031206;
			break;
			
  }

	// Avoid negative values for w1 and w2
	if (w1 < 0) {
		w1 = 0;
	}
	
	if (w2 < 0) {
		w2 = 0;
	}

  // Done!
  return arrayOfParameters;
}

double
FunctionK(const double theta,
          const std::array<double, 4> parameters)
{
	// Note this function requires theta > 0 for nu < 1, but there is not checked here!
  return std::pow(theta, parameters[0]-1) * (
  	std::exp(-theta * parameters[1])
  	+ parameters[3] * std::exp(-theta*parameters[2])
  );
}

double
SmallAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters)
{
	// We do not check values of the input parameters here! It is required that:
	// 0 <= lowAngle,highAngle <= thetaEm <= pi
	// numberOfTerms > 0 can be small
	
	// Let this number be fixed by now
	static const int numberOfTerms = 10;
	
  double integral = 0;

  double cosThetaEm = std::cos(thetaEm);
  double sinThetaEm = std::sin(thetaEm);

	// Avoid the computation of 0*log(0)
  if (thetaEm - lowAngle > 1e-20) {
    integral += (std::cos(lowAngle) - cosThetaEm) 
      * (Constants::Pi - std::log(1.0 - lowAngle/thetaEm));
  }

	// Avoid the computation of 0*log(0)
  if (thetaEm - highAngle > 1e-20) {
    integral -= (std::cos(highAngle) - cosThetaEm)
      * (Constants::Pi - std::log(1.0 - highAngle/thetaEm));
  }

  // Compute the power series
  double deltaLow = thetaEm - lowAngle;
  double deltaHigh = thetaEm - highAngle;
  double powerLow = 1;
  double powerHigh = 1;
  double factorial = 1;
  double multiplier[4] = {cosThetaEm, sinThetaEm, -cosThetaEm, -sinThetaEm};

  for (int k = 1; k <= numberOfTerms; k++) {
    powerLow  *= deltaLow;
    powerHigh *= deltaHigh;
    factorial *= k;

    integral -= multiplier[k%4] * (powerHigh - powerLow) / (k * factorial);
  }

  integral *= FunctionK(thetaEm, parameters) / sinThetaEm;

  return integral;
}



double
LargeAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters)
{
	// This is a default step size we use during integration. In practice the value will be close to
	// this one but somewhat smaller. Note that errors from the composite Simpson's rule are 
	// proportional to stepSize^4
	static const double defaultStepSize = 1e-3;

	// Parametrized quantities
  const double & nu = parameters[0];
  const double & w1 = parameters[1];
  const double & w2 = parameters[2];
  const double & ep = parameters[3];

	// Auxiliary parameters
  const double r1 = thetaEm * w1;
  const double r2 = thetaEm * w2;
  const double rRight = thetaEm / lowAngle;
  const double rLeft = thetaEm / highAngle;

	// Compute number of steps by ensuring that stepSize <= defaultStepSize, but close to it
  int numberOfSteps = 1 + (rRight - rLeft) / defaultStepSize;
  double stepSize = (rRight - rLeft) / double(numberOfSteps);
  double halfStep = 0.5*stepSize;

	// Auxiliary function: computes x*(pi + 1 - log(x)), taking care of values close to 0 (assumes x >= 0)
	auto xPiPlusOneLogx = [](const double x)->double{
		// If x is too close to 0, we expand x*log(x) in a Taylor series and take the first order term only
    // That is: x*log(x) ~ x^x - 1
    // In this way we ensure that the function is always finite, because the C++ standard defines x^0 = 1,
    // whatever the value of 0 (even NaN).
		if (x <= 1e-10) {
			return Constants::PiPlusOne * x + 1.0 - std::pow(x, x);
		} else {
			return x * (Constants::PiPlusOne - std::log(x));
		}
	};

	// Auxiliary function that computes the integrand
  auto integrand = [&](const double & r)->double{

    double value = 0;
    value += (nu + 1 - r1/r) * std::exp(-r1/r);
    value += ep * (nu + 1 - r2/r) * std::exp(-r2/r);
    value *= xPiPlusOneLogx(1-r);
    value /= std::pow(r, nu + 2);

    return value;
  };

	// Initial value for the integral
  double integral = 0;
  
  // This is composite Simpson's rule for the integral that doesn't have a closed form expression
  for (int iStep = 0; iStep < numberOfSteps; iStep++) {
  	integral +=
  	  integrand(rLeft + iStep * stepSize)
  		+ 2 * integrand(rLeft + (2*iStep + 1) * halfStep);
  }
  integral = - ( 2 * integral + integrand(rRight) - integrand(rLeft) ) * stepSize / 6.0;
  
  // Below terms are from the integral by parts we did to reach a finite integrate for all r
  integral += xPiPlusOneLogx(1-rLeft)
  	* (std::exp(-r1/rLeft) + ep*std::exp(-r2/rLeft))
  	/ std::pow(rLeft, nu+1);
  	
  integral -= xPiPlusOneLogx(1-rRight)
  	* (std::exp(-r1/rRight) + ep*std::exp(-r2/rRight))
  	/ std::pow(rRight, nu+1);

	// This factor comes after the change of variables theta -> r = thetaEm / theta
  integral *= std::pow(thetaEm, nu);

  return integral;
}

double
CherenkovUnnormalizedIntegral(const double lowAngle,
                              const double highAngle,
                              const double thetaEm, 
                              const std::array<double,4> parameters)
{
	// This function requires lowAngle and highAngle to be consistent with 0 <= lowAngle < highAngle
  if (highAngle <= thetaEm) {
    // integration interval is contained in [0,thetaEm]
    return SmallAngleIntegral(lowAngle, highAngle, thetaEm, parameters);
  } else if (lowAngle >= thetaEm) {
    // integration interval is contained in [thetaEm, infinity)
    return LargeAngleIntegral(lowAngle, highAngle, thetaEm, parameters);
  } else {
    // integration interval is contained in [0, "infinity") with thetaEm in (lowAngle,highAngle)
    // so we split the computation in two parts, each correspoding to one of the cases above
    return SmallAngleIntegral(lowAngle, thetaEm, thetaEm, parameters)
      + LargeAngleIntegral(thetaEm, highAngle, thetaEm, parameters);
  }
}

double
CherenkovUnnormalizedDensity(const double theta,
                             const double thetaEm, 
                             const std::array<double,4> parameters)
{
  static const double delta = 1e-5;
  
  if (theta < thetaEm - delta) {
 		// when theta is smaller than thetaEm and no too close from it, compute the density as is
    return ( std::sin(theta) / std::sin(thetaEm) )
      * (Constants::Pi - std::log(1.0 - theta/thetaEm) )
      * FunctionK(thetaEm, parameters);
  } else if (theta < thetaEm + delta) {
  	// when theta is too close from thetaEm, return the average value so to avoid a divergence
    return CherenkovUnnormalizedIntegral(thetaEm - delta, thetaEm + delta, thetaEm, parameters) / (2*delta);
  } else {
  	// for theta larger than thetaEm and not too close from it, compute the density as defined
    return (Constants::Pi - std::log(1.0 - thetaEm/theta))
      * FunctionK(theta, parameters);
  }
}



double
CherenkovPDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV,
             const ParticleType primaryParticle)
{
  // Check if theta is outside the range within which the function is defined
  if (0 <= theta || theta >= Constants::MaxTheta) {
    return 0;
  }
  
  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  // Compute the maximum Cherenkov emission angle
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedDensity(theta, thetaEm, parameters)
    / CherenkovUnnormalizedIntegral(0.0, Constants::MaxTheta, thetaEm, parameters);
}



double
CherenkovCDF(const double theta,
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

  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  // Compute the maximum Cherenkov emission angle, used in the function below
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedIntegral(0.0, theta, thetaEm, parameters)
    / CherenkovUnnormalizedIntegral(0.0, Constants::MaxTheta, thetaEm, parameters);
}



double
CherenkovIntegral(const double lowAngle,
                  const double highAngle,
                  const double showerAge,
                  const double refractiveIndex,
                  const double showerEnergyTeV,
                  const ParticleType primaryParticle)
{
	// Check input angles
	if (lowAngle > highAngle) {
		return -CherenkovIntegral(highAngle, lowAngle, showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);
	} else if (highAngle <= 0) {
		return 0;
	}

  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV, primaryParticle);

  // Compute the maximum Cherenkov emission angle, used in the function below
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedIntegral(
  	std::max(lowAngle,0.0), 
  	std::min(highAngle,Constants::MaxTheta),
  	thetaEm,
  	parameters
  ) / CherenkovUnnormalizedIntegral(0.0, Constants::MaxTheta, thetaEm, parameters);
}
