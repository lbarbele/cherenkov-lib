#pragma once

#include <array>

// Functions provided by this library
std::array<double, 4> 
CherenkovParameters(const double showerAge,
                    const double refractiveIndex,
                    const double showerEnergyTeV);

double
CherenkovPDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV);

double
CherenkovCDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV);

double
CherenkovIntegral(const double lowAngle,
                  const double highAngle,
                  const double showerAge,
                  const double refractiveIndex,
                  const double showerEnergyTeV);

// Internal functions: remove them from here and make them static in cerang.cpp
double
FunctionK(const double theta,
          const std::array<double, 4> parameters);

double
SmallAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters,
                   const int numberOfTerms = 10);

double
LargeAngleIntegral(const double lowAngle,
                   const double highAngle,
                   const double thetaEm,
                   const std::array<double, 4> parameters,
                   const double defaultStepSize = 1e-3);

double
CherenkovUnnormalizedDensity(const double theta,
                             const double thetaEm, 
                             const std::array<double,4> parameters);

double
CherenkovUnnormalizedIntegral(const double lowAngle,
                              const double highAngle,
                              const double thetaEm, 
                              const std::array<double,4> parameters);
