#pragma once

#include <array>

namespace Cherenkov {
	namespace Longi {
		double
		YieldFactor(
			const double showerAge, 
			const double cutoffEnergy, 
			const double refractiveIndex,
			const double lambdaMin,
			const double lambdaMax);
	}
}
