#pragma once

#include <cmath>

enum class ParticleType {
	Proton,
	Gamma
};

namespace Cherenkov {
	namespace Angular {
		double
		PDF(
			const double theta,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);
				         
		double
		CDF(
			const double theta,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);
				         
		double
		Integral(
			const double lowAngle,
			const double highAngle,
			const double showerAge,
			const double refractiveIndex,
			const double showerEnergyTeV,
			const ParticleType primaryParticle);
	}
}
