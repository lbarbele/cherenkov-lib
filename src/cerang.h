#pragma once

#include <cmath>

enum class ParticleType {
	Proton,
	Gamma
};

namespace Constants {
	inline constexpr double Pi = 3.14159265358979311599796346854;
	inline constexpr double PiPlusOne = Pi + 1.0;
	inline constexpr double MaxTheta = 0.5*Pi;
}

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
