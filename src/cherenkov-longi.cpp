#include <array>
#include <cmath>

#include <gsl/gsl_sf_hyperg.h>

#include <cherenkov-longi.h>
#include <constants.h>

//
//
//
double
Cherenkov::Longi::YieldFactor(
	const double showerAge,
	const double cutoffEnergy,
	const double refractiveIndex,
	const double lambdaMin,
	const double lambdaMax)
{
	const double a1 = 6.42522 - 1.53183 * showerAge;
	const double a2 = 168.168 - 42.1368 * showerAge;
	
	const double normalization = gsl_sf_hyperg_2F1(1, showerAge, 1 + showerAge, (a2-a1)/(a2+cutoffEnergy))
			/ ( showerAge	* std::pow(a2 + cutoffEnergy, showerAge) );
	
	const double thresholdEnergy = Constants::ElectronMass
		* refractiveIndex 
		/ std::sqrt( (refractiveIndex + 1) * (refractiveIndex - 1) );
		
	double factor = 0;
	
	factor += (1 - thresholdEnergy/a1) * (1 + thresholdEnergy/a1)
		* gsl_sf_hyperg_2F1(1, showerAge, 1+showerAge, (a2-a1)/(a2+thresholdEnergy));
		
	factor += (thresholdEnergy/a1) * (thresholdEnergy/a1)
		* (1 + showerAge * a1 / a2)
		* gsl_sf_hyperg_2F1(1, showerAge, 1+showerAge, a2/(a2+thresholdEnergy));
		
	factor -= showerAge * (thresholdEnergy / a1) * (1 + thresholdEnergy / a2);
	
	factor *= 4 * Constants::Pi * Constants::AlphaEm;
	factor *= (1.0/lambdaMin - 1.0/lambdaMax);
	
	factor /= normalization
		* showerAge
		* std::pow(a2 + thresholdEnergy, showerAge);
	
	return factor;
}
