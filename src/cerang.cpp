#include "cerang.h"

#include <cmath>
#include <iostream>

static const double kPiOnTwo = std::asin(1.0);
static const double kPi = std::acos(-1.0);
static const double kPiPlusOne = kPi + 1;

//
// Helper functions
//
double
FunctionK(const double theta,
          const std::array<double, 4> parameters)
{
  return std::pow(theta, parameters[0]-1)
    * std::exp(-theta/parameters[1])
    * (1 + parameters[3] * std::exp(theta/parameters[2]));
}

double
SmallAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters,
                   const int numberOfTerms)
{

  double integral = 0;

  double cosThetaEm = std::cos(thetaEm);
  double sinThetaEm = std::sin(thetaEm);

  if (lowAngle != thetaEm) {
    integral += (std::cos(lowAngle) - cosThetaEm) 
      * (kPi - std::log(1 - lowAngle/thetaEm));
  }

  if (highAngle != thetaEm) {
    integral -= (std::cos(highAngle) - cosThetaEm)
      * (kPi - std::log(1 - highAngle/thetaEm));
  }

  // Compute the power series
  double deltaLow = thetaEm - lowAngle;
  double deltaHigh = thetaEm - highAngle;
  double powerLow = 1;
  double powerHigh = 1;
  double factorial = 1;
  double multiplier[4] = {cosThetaEm, sinThetaEm, -cosThetaEm, -sinThetaEm};

  for (int k = 1; k <= 2; k++) {
    powerLow  *= deltaLow;
    powerHigh *= deltaHigh;
    factorial *= k;

    integral -= multiplier[k%4] * (powerHigh - powerLow) / (k * factorial);
  }

  integral *= FunctionK(thetaEm, parameters)/sinThetaEm;

  return integral;
}

double
LargeAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters,
                   const double defaultStepSize)
{

  const double & nu = parameters[0];
  const double & t1 = parameters[1];
  const double & t2 = parameters[2];
  const double & ep = parameters[3];

  double ra = thetaEm / t1;
  double rb = t2 < t1 ? 0 : ra * (1.0 - t1/t2);
  double rl = thetaEm / lowAngle;
  double rr = thetaEm / highAngle;

  int numberOfSteps = 1 + (rl - rr) / defaultStepSize;
  double stepSize = (rl - rr) / double(numberOfSteps);
  double halfStep = 0.5*stepSize;

  auto integrand = [&](const double & r)->double{

    if (1-r <= 0) {
      return 0;
    }

    double value = 0;
    value += (nu + 1 - ra/r) * std::exp(-ra/r);
    value += ep * (nu + 1 - rb/r) * std::exp(-rb/r);
    value *= kPiPlusOne - std::log(1-r);
    value *= 1-r;
    value /= std::pow(r, nu + 2);

    return value;
  };

  double integral = 0;

  double functionRight = integrand(rr);
  double functionLeft = 0;

  for (double r = rr + halfStep; r < rl; r += stepSize) {
    functionLeft = functionRight;
    functionRight = integrand(r+halfStep);

    integral -= (functionLeft + 4*integrand(r) + functionRight);
  }
  integral *= stepSize / 6.0;

  integral += (1-rr)
    * (kPiPlusOne - std::log(1-rr))
    * (std::exp(-ra/rr) + ep*std::exp(-rb/rr))
    / std::pow(rr, nu+1);

  if (rl != 1) {
    integral -= (1-rl)
      * (kPiPlusOne - std::log(1-rl))
      * (std::exp(-ra/rl) + ep*std::exp(-rb/rl))
      / std::pow(rl, nu+1);
  }

  integral *= std::pow(thetaEm, nu);

  return integral;
}

double
CherenkovUnnormalizedDensity(const double theta,
                             const double thetaEm, 
                             const std::array<double,4> parameters)
{
  static const double delta = 1e-5;

  if (theta < 0) {
    return 0;
  } else if (theta < thetaEm - delta) {
    return ( std::sin(theta) / std::sin(thetaEm) )
      * (kPi - std::log(1 - theta/thetaEm) )
      * FunctionK(thetaEm, parameters);
  } else if(theta < thetaEm + delta) {
    return CherenkovUnnormalizedIntegral(thetaEm-delta, thetaEm+delta, thetaEm, parameters) / (2 * delta);
  } else {
    return (kPi - std::log(1 - thetaEm/theta))
      * FunctionK(theta, parameters);
  }
}

double
CherenkovUnnormalizedIntegral(const double lowAngle,
                              const double highAngle,
                              const double thetaEm, 
                              const std::array<double,4> parameters)
{

  if (lowAngle > highAngle) {
    return -CherenkovUnnormalizedIntegral(highAngle, lowAngle, thetaEm, parameters);
  } else if (highAngle < 0) {
    return 0.;
  } else if (lowAngle < 0) {
    return CherenkovUnnormalizedIntegral(0, highAngle, thetaEm, parameters);
  }

  if (highAngle <= thetaEm) {
    // integration interval is contained in [0,thetaEm]
    return SmallAngleIntegral(lowAngle, highAngle, thetaEm, parameters);
  } else if (lowAngle >= thetaEm) {
    // integration interval is contained in [thetaEm, "infinity")
    return LargeAngleIntegral(lowAngle, highAngle, thetaEm, parameters);
  } else {
    // integration interval is contained in [0, "infinity")
    // so we split the compmutation to fall into both cases above separately
    return SmallAngleIntegral(lowAngle, thetaEm, thetaEm, parameters)
      + LargeAngleIntegral(thetaEm, highAngle, thetaEm, parameters);
  }
}

//
// Functions defined in the API
//
std::array<double, 4> 
CherenkovParameters(const double showerAge,
                    const double refractiveIndex,
                    const double showerEnergyTeV)
{
  // No checks on the values of the arguments are made here, aiming for
  // performance. Those are expected to be done externally. It is expected that:
  //
  // showerAge > 0 (a value of 0 would make log(s) diverge)
  // refractiveIndex > 1 (the pow() will return NaN if refractiveIndex = 1)
  // showerEnergy > 0 (pow() of showerEnergy = 0 is also NaN)

  // Creat the object that will be returned
  std::array<double, 4> arrayOfParameters = {0, 0, 0, 0};

  // Check values of age and refractive index
  if (showerAge <= 0 || showerAge > 3 || refractiveIndex <= 1) {
    return arrayOfParameters;
  }

  // Give friendly names to the elements of the array, for readability
  auto & nu = arrayOfParameters[0];
  auto & t1 = arrayOfParameters[1];
  auto & t2 = arrayOfParameters[2];
  auto & ep = arrayOfParameters[3];

  // Compute auxiliar quantities
  double delta = refractiveIndex - 1.0;
  double logAge = std::log(showerAge);

  // Parametrization for primary gamma:
  nu = 0.34329 * std::pow(delta, -0.10683) + 1.46852 * logAge;
  t1 =  1.4053 * std::pow(delta, 0.32382) - 0.048841 * logAge;
  t2 = t1 * (0.95734 + 0.26472 * showerAge);
  ep = 0.0031206;

  // Parametrization for primary proton:
  // nu = 0.21155 * std::pow(delta, -0.16639) + 1.21803 * logAge;
  // t1 = 4.513 * std::pow(delta, 0.45092) * std::pow(showerEnergyTeV, -0.008843) - 0.058687* logAge;
  // t2 = t1 * (0.90725 + 0.41722 * showerAge);
  // ep = 0.009528 + 0.022552 * std::pow(showerEnergyTeV, -0.4207);

  // Check parameter values

  // From the form of our parametrization, nu, in principle can assume and value
  // and epsilon (here "ep") is always positive. Therefore we don't need to check
  // these parameters (unless some unexpected error occurs!).

  // As t1 cannot be <= 0, we set a minimum value to it. Numerically, this limit is
  // "almost" equivalent to t1 = 0
  if (t1 <= 0) {
    t1 = 1e-10;
  }

  // We cannot have t2 < t1 because, if we did, our density function would grow
  // as propto exp(theta) for large values of theta. Thus, we establish t2 >= t1
  if (t2 < t1) {
    t2 = t1;
  }

  // Done!
  return arrayOfParameters;
}



double
CherenkovPDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV)
{
  // Check if theta is outside range where the density function is defined
  if (theta <= 0) {
    return 0;
  } else if (theta >= kPiOnTwo) {
    return 1;
  }

  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);

  // Compute the maximum Cherenkov emission angle, used in the function below
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedDensity(theta, thetaEm, parameters)
    / CherenkovUnnormalizedIntegral(0.0, kPiOnTwo, thetaEm, parameters);
}



double
CherenkovCDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV)
{
  // Check if theta is outside range where the density function is defined
  if (theta <= 0) {
    return 0;
  } else if (theta >= kPiOnTwo) {
    return 1;
  }

  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);

  // Compute the maximum Cherenkov emission angle, used in the function below
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedIntegral(0.0, theta, thetaEm, parameters)
    / CherenkovUnnormalizedIntegral(0.0, kPiOnTwo, thetaEm, parameters);
}



double
CherenkovIntegral(const double lowAngle,
                  const double highAngle,
                  const double showerAge,
                  const double refractiveIndex,
                  const double showerEnergyTeV)
{
  // Compute the parameters of our parametrized function
  auto parameters = CherenkovParameters(showerAge, refractiveIndex, showerEnergyTeV);

  // Compute the maximum Cherenkov emission angle, used in the function below
  double thetaEm = std::acos(1.0/refractiveIndex);

  return CherenkovUnnormalizedIntegral(lowAngle, std::max(highAngle, kPiOnTwo), thetaEm, parameters)
    / CherenkovUnnormalizedIntegral(0.0, kPiOnTwo, thetaEm, parameters);
}