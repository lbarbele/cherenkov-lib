#include <cherenkov-angular.h>
#include <constants.h>
#include <math.h>

#include <cmath>

Cherenkov::AngularDistribution::AngularDistribution(
  const double showerEnergyTeV,
  const ParticleType particleType,
  const double maxTheta
) : 
  fEnergyTeV(showerEnergyTeV),
  fParticleType(particleType),
  fMaxTheta(maxTheta),
  fShowerAge(-1),
  fRefractiveIndex(-1)
{
}



void
Cherenkov::AngularDistribution::UpdateParameters(
	const double showerAge,
	const double refractiveIndex
)
{
  // Parameters are recomputed only if shower age or refractive index change
  if (showerAge == fShowerAge && refractiveIndex == fRefractiveIndex) {
    return;
  }
  
  // Copy the values of shower age and refractive index
  fShowerAge = showerAge;
  fRefractiveIndex = refractiveIndex;
  
  // Truncate the values of shower age and refractive index within the following bounds:
  // 0.5 <= showerAge <= 1.5
  // refractiveIndex >= 1.0 + 3e-5
	const double nMinusOne = std::max(refractiveIndex, 1.0 + 3e-5) - 1.0;
	const double truncatedAge = std::max(0.5, std::min(showerAge, 1.5));
	const double logShowerAge = std::log(truncatedAge);
	
	// Compute the parametrized parameters according to the particle type
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
	
	// Compute thetaEm (maximum Cherenkov emission angle) and its sines and cosines
	const double v = 1. - 1./refractiveIndex;
	fThetaEm = v < 0.001 ? std::sqrt(2.*v)*(1.+(v/12.)*(1.+9.*v/40.)) : std::acos(1./refractiveIndex);
  fTwoThetaEm = 2.*fThetaEm;
	fCosThetaEm = 1./refractiveIndex;
	fSinThetaEm = std::sqrt((1.-fCosThetaEm)*(1.+fCosThetaEm));
	fCosTwoThetaEm = (2.*fCosThetaEm - refractiveIndex) * fCosThetaEm;
	fSinTwoThetaEm = 2.*fSinThetaEm*fCosThetaEm;
	
	// Compute the normalization factor
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
Cherenkov::AngularDistribution::UnnormalizedIntegralLeft(
	const double lowAngle,
	const double highAngle
) const
{
  const double deltaLow = fThetaEm-lowAngle;
  const double deltaHigh = fThetaEm-highAngle;
  const double deltaCosLow = std::cos(lowAngle) - fCosThetaEm;
  const double deltaCosHigh = std::cos(highAngle) - fCosThetaEm;
  const double coeffs[4] = {fCosThetaEm, fSinThetaEm, -fCosThetaEm, -fSinThetaEm};

  double integral = 0;
  double powLow = 1;
  double powHigh = 1;
  double term;
  for (int k = 1; k < 1000; k++) {
    powLow *= deltaLow/k;
    powHigh *= deltaHigh/k;
    term = (powLow - powHigh) / k;
    integral += coeffs[k%4] * term;
    if (std::fabs(term) <= 1e-10*std::fabs(integral))
      break;
  }

  integral += Constants::Pi * (deltaCosLow - deltaCosHigh);
  integral += deltaHigh < 1e-10 ? 
    std::log(std::pow(deltaHigh/fThetaEm+1e-200,deltaCosHigh)) :
    deltaCosHigh * std::log(deltaHigh/fThetaEm); 
  integral -= deltaLow < 1e-10 ? 
    std::log(std::pow(deltaLow/fThetaEm+1e-200,deltaCosLow)) :
    deltaCosLow * std::log(deltaLow/fThetaEm);
    
  integral *= FunctionK(fThetaEm) / fSinThetaEm;
    
  return integral;
}



double
Cherenkov::AngularDistribution::UnnormalizedIntegralRight(
	const double lowAngle,
	const double highAngle
) const
{
  // Integrand (after the change t = sqrt(1-thetaEm/theta)) used for numerical integration
  auto integrand = [=](double t) {
    const double theta = fThetaEm / ((1.0-t)*(1.0+t));
    return std::pow(theta, fNu + 1) *
	    (t <= 1e-10 ? t*Constants::HalfPi + 1.0 - std::pow(t, t) : t * (Constants::HalfPi - std::log(t))) *
	    (std::exp(-fOmega1*theta) + fEpsilon*std::exp(-fOmega2*theta));
  };
  
  // Integral given as an expansion in powers of (1 - thetaEm/theta)
  auto uSeries = [=](const double theta, const double omega)->double{
    const double u = 1.0 - fThetaEm/theta;
    if (u < 1e-18)
      return 0;
  
    const double pimlog = u > 1e-5 ? u*(Constants::Pi - std::log(u)) :
      u*Constants::Pi - std::log(std::pow(u+1e-200,u));
    
    double bm2 = 0;
    double bm1 = std::exp(-omega*fThetaEm);
    double sum = (pimlog+u)*bm1;
    double a1 = fNu-1-omega*fThetaEm;
    double a2 = -1;
    double b, term, upn;
    for (int n = 2; n < 100; n++) {
      a1 += 2;
      a2++;
      b = a1*bm1 - a2*(fNu+a2)*bm2*upn;
      upn = u/n;
      b *= upn;
      term = b * (pimlog + upn);
      sum += term;
      if (std::fabs(term) <= 1e-10*std::fabs(sum))
        break;
      bm2 = bm1;
      bm1 = b;
    }

    return sum * std::pow(fThetaEm, fNu);
  };
  
  // Integral given as an expansion in powers of thetaEm/theta
  auto rSeries = [=](const double theta, const double omega)->double{
    const double wt = omega*theta;
    const double r = fThetaEm / theta;
    const double wtem = fThetaEm * omega;
    
    double h = 0;
    
    if (wt > fNu + 1) {
      const double fpMin = 1e-30;
      double b = wt + 1 - fNu;
      double c = 1.0 / fpMin;
      double d = 1.0 / b;
      h = -d;
      for (int n = 1; n < 1000; n++) {
        double an = -n * (n - fNu);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < fpMin)
          d = fpMin;
        c = b + an/c;
        if (std::fabs(c) < fpMin)
          c = fpMin;
        d = 1.0/d;
        double del = d*c;
        h *= del;
        if (std::fabs(del - 1) < 1e-10)
          break;
      }
    } else {
      double term = 1.0/fNu;
      h = term;
      for (int n = 1; n < 1000; n++) {
        term *= wt / (n + fNu);
        h += term;
        if (std::fabs(term) < 1e-10*std::fabs(h))
          break;
      }
      h -= std::exp(wt - fNu*std::log(wt) + std::lgamma(fNu));
    }
    
    double rPower = 1;
    double sum = Constants::Pi * h;
    for (int k = 1; k < 1000; k++) {
      rPower *= r;
      h = (rPower + wtem*h) / (fNu - k);
      sum += h / double(k);
      if (std::fabs(h / double(k)) <= 1e-10*std::fabs(sum))
        break;
    }
    
    return sum * std::exp(-wt + fNu * std::log(theta));
  };
  
  const double thresholdAngle = 2*fThetaEm;

  double integral = 0;
  
  if (lowAngle < thresholdAngle) {
    const double actualHighAngle = std::min(highAngle, thresholdAngle);
    integral += uSeries(actualHighAngle, fOmega1) - uSeries(lowAngle, fOmega1) +
      fEpsilon * (uSeries(actualHighAngle, fOmega2) - uSeries(lowAngle, fOmega2));
  }

  if (highAngle > thresholdAngle) {
    const double actualLowAngle = std::max(lowAngle, thresholdAngle);
    if (std::fabs(fNu - std::round(fNu)) < 0.001) {
      const double tMin = std::sqrt(1.0 - fThetaEm/actualLowAngle);
      const double tMax = std::sqrt(1.0 - fThetaEm/highAngle);
      integral += (4.0 / fThetaEm) * Math::RombergIntegral(tMin, tMax, 1e-10, integrand);
    } else {
      integral += rSeries(highAngle, fOmega1) - rSeries(actualLowAngle, fOmega1) +
        fEpsilon * (rSeries(highAngle, fOmega2) - rSeries(actualLowAngle, fOmega2));
    }
  }
  
	return integral;
}



double
Cherenkov::AngularDistribution::UnnormalizedIntegral(
	const double lowAngle,
	const double highAngle
) const
{
	double integral = 0;

	// if lowAngle is to the left of the peak, integrate between lowAngle and min(thetaEm, highAngle)
	if (lowAngle < fThetaEm) {
    integral += UnnormalizedIntegralLeft(lowAngle, std::min(highAngle, fThetaEm));
	}
	
	// If highAngle is to the right of the peak, integrate between max(thetaEm,lowAngle) and highAngle
	if (highAngle > fThetaEm) {
	  integral += UnnormalizedIntegralRight(std::max(lowAngle, fThetaEm), highAngle);
	}
	
	return integral;
}



double
Cherenkov::AngularDistribution::UnnormalizedIntegralSineLeft(
	const double lowAngle,
	const double highAngle
) const
{
  const double deltaLow = 2.0*(fThetaEm - lowAngle);
  const double deltaHigh = 2.0*(fThetaEm - highAngle);
  const double factorLow = fSinTwoThetaEm - std::sin(2*lowAngle) - deltaLow;
  const double factorHigh = fSinTwoThetaEm - std::sin(2*highAngle) - deltaHigh;
  const double coeffs[4] = {fSinTwoThetaEm, -fCosTwoThetaEm, -fSinTwoThetaEm, fCosTwoThetaEm};
  
  double integral = 0;
  double powLow = 1;
  double powHigh = 1;
  double term, onePerk;
  for (int k = 1; k < 1000; k++) {
    powLow *= deltaLow / k;
    powHigh *= deltaHigh / k;
    term = (powLow - powHigh) / k;
    integral += coeffs[k%4] * term;
    if (std::fabs(term) <= 1e-10*std::fabs(integral))
    break;
  }
  
  integral += deltaLow - deltaHigh;
  integral += Constants::Pi * (factorHigh-factorLow);
  integral -= deltaHigh < 1e-10 ?
    std::log(std::pow(deltaHigh/fTwoThetaEm+1e-200,factorHigh)) :
    factorHigh * std::log(deltaHigh/fTwoThetaEm);
  integral += deltaLow < 1e-10 ?
    std::log(std::pow(deltaLow/fTwoThetaEm+1e-200,factorLow)) :
    factorLow * std::log(deltaLow/fTwoThetaEm);
  integral *= 0.25 * FunctionK(fThetaEm) / fSinThetaEm;
  
  return integral;
}



double
Cherenkov::AngularDistribution::UnnormalizedIntegralSineRight(
	const double lowAngle,
	const double highAngle
) const
{
  // for numerical integration
  auto integrand = [=](double t) {
	  const double theta = fThetaEm / ((1.0-t)*(1.0+t));
	  return std::sin(theta) * std::pow(theta, fNu + 1) *
		  (t <= 1e-10 ? t*Constants::HalfPi + 1.0 - std::pow(t, t) : t * (Constants::HalfPi - std::log(t))) *
		  (std::exp(-fOmega1*theta) + fEpsilon*std::exp(-fOmega2*theta));
  };
  
  // Integral given as an expansion in powers of (1 - thetaEm/theta)
  auto uSeries = [=](const double theta, const double omega)->double{
    const double u = 1.0 - fThetaEm/theta;
    if (u < 1e-18)
      return 0;
  
    const double wtem = omega*fThetaEm;
    const double pimlog = u > 1e-5 ? u*(Constants::Pi - std::log(u)) :
      u*Constants::Pi - std::log(std::pow(u+1e-200,u));
    
    double bm2 = 0;
    double am2 = 0;
    double bm1 = std::exp(-omega*fThetaEm);
    double am1 = bm1 * fCosThetaEm;
    bm1 *= fSinThetaEm;
    double sum = (pimlog + u) * bm1;
    double c1 = fNu-1-omega*fThetaEm;
    double c2 = -1;
    double upn = u;
    double a, b, term;
    for (int n = 2; n < 100; n++) {
      c1 += 2;
      c2++;
      a = c1*am1 - c2*(fNu+c2)*am2*upn - fThetaEm*bm1;
      b = c1*bm1 - c2*(fNu+c2)*bm2*upn + fThetaEm*am1;
      upn = u/n;
      a *= upn;
      b *= upn;
      term = (pimlog + upn) * b;
      sum += term;
      if (std::fabs(term) <= 1e-10*std::fabs(sum))
        break;
      bm2 = bm1;
      bm1 = b;
      am2 = am1;
      am1 = a;
    }

    return sum * std::pow(fThetaEm, fNu);
  };
  
  // Integral given as an expansion in powers of thetaEm/theta
  auto rSeries = [=](const double theta, const double omega) {
    const double r = fThetaEm/theta;
    const double sint = std::sin(theta);
    const double cost = std::cos(theta);
    const double invt = 1.0/theta;
    
    double kap = 1.0/fNu;
    double sig = 0;
    double old;
    double c = kap;
    double s = sig;
    double den = fNu * invt;
    for (int n = 1; n < 1000; n++) {
      den += invt;
      old = kap;
      kap = (omega*kap - sig) / den;
      sig = (omega*sig + old) / den;
      c += kap;
      s += sig;
      if (std::fabs(kap) <= 1e-10*std::fabs(c) && std::fabs(sig) <= 1e-10*std::fabs(s))
        break;
    }
    
    double b = sint*c - cost*s;
    double rp = 1;
    double series = Constants::Pi * b;
    den = fNu;
    for (int k = 1; k < 1000; k++) {
      den--;
      rp *= r;
      old = s;
      s = (omega*s + c) * fThetaEm / den;
      c = (rp + fThetaEm*(omega*c - old)) / den;
      b = (sint*c - cost*s) / k;
      series += b;
      if (std::fabs(b) < 1e-10*std::fabs(series))
        break;
    }
    
    return series * std::exp(-omega*theta + fNu*std::log(theta));
  };
  
  const double thresholdAngle = 2*fThetaEm;

  double integral = 0;
  
  if (lowAngle < thresholdAngle) {
    const double actualHighAngle = std::min(highAngle, thresholdAngle);
    integral += uSeries(actualHighAngle, fOmega1) - uSeries(lowAngle, fOmega1) +
      fEpsilon * (uSeries(actualHighAngle, fOmega2) - uSeries(lowAngle, fOmega2));
  }

  if (highAngle > thresholdAngle) {
    const double actualLowAngle = std::max(lowAngle, thresholdAngle);
    if (std::fabs(fNu - std::round(fNu)) < 0.001) {
      const double tMin = std::sqrt(1.0 - fThetaEm/actualLowAngle);
      const double tMax = std::sqrt(1.0 - fThetaEm/highAngle);
      integral += (4.0 / fThetaEm) * Math::RombergIntegral(tMin, tMax, 1e-10, integrand);
    } else {
      integral += rSeries(highAngle, fOmega1) - rSeries(actualLowAngle, fOmega1) +
        fEpsilon * (rSeries(highAngle, fOmega2) - rSeries(actualLowAngle, fOmega2));
    }
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
  static const double peakWidth = 1e-5;

	// Check the value of theta
	if (theta <= 0) {
		return 0;
	} else if (theta > fMaxTheta) {
		return 0;
  }

	// Update the parameters that depend on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	if (theta <= fThetaEm - peakWidth) {
		// Region to the left of the peak: compute the distribution function as defined
		return 	FunctionK(fThetaEm) * (Constants::Pi - std::log(1.0-theta/fThetaEm)) *
			(std::sin(theta)/fSinThetaEm) / fNormalization;
	}
	else if (theta < fThetaEm + peakWidth) {
		// Very close to the peak: compute the average in an interval containing the peak
		return UnnormalizedIntegral(fThetaEm - peakWidth, fThetaEm + peakWidth) /
			(2.0*peakWidth*fNormalization);
	}
	else {
		// Right of the peak: compute the distribution function as defined
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
	if (theta <= 0) {
		return 0;
	} else if (theta >= fMaxTheta) {
		return 1;
	}
	
	// Update the parameters that depend on refractive index and shower age
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
  // Reorder angles, if necessary
	if (lowAngle > highAngle) {
		return -Integral(highAngle, lowAngle, showerAge, refractiveIndex);
	}
	
	// Update the parameters that depend on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	// Return the (normalized) integral of the distribution between the given angles
	return UnnormalizedIntegral(std::max(lowAngle,0.), std::min(highAngle,fMaxTheta)) /
		fNormalization;
}



double
Cherenkov::AngularDistribution::IntegralSine(
	const double lowAngle, 
	const double highAngle, 
	const double showerAge,
	const double refractiveIndex
)
{
	// Reorder angles, if necessary
	if (lowAngle > highAngle) {
		return -IntegralSine(highAngle, lowAngle, showerAge, refractiveIndex);
	}
	
	// Update the parameters that depend on refractive index and shower age
	UpdateParameters(showerAge, refractiveIndex);
	
	double integral = 0;
	
	if (lowAngle < fThetaEm) {
	  integral += UnnormalizedIntegralSineLeft(std::max(lowAngle,0.), std::min(highAngle,fThetaEm));
	}
	
	if (highAngle > fThetaEm) {
	  integral += UnnormalizedIntegralSineRight(std::max(lowAngle,fThetaEm), std::min(highAngle,fMaxTheta));
	}
	
	// return the normalized integral
	return integral / fNormalization;
}


