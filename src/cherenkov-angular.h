#pragma once

#include <cmath>

enum class ParticleType {
	Proton,
	Gamma
};

double
CherenkovPDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV,
             const ParticleType primaryParticle);
             
double
CherenkovCDF(const double theta,
             const double showerAge,
             const double refractiveIndex,
             const double showerEnergyTeV,
             const ParticleType primaryParticle);
             
double
CherenkovIntegral(const double lowAngle,
                  const double highAngle,
                  const double showerAge,
                  const double refractiveIndex,
                  const double showerEnergyTeV,
                  const ParticleType primaryParticle);
