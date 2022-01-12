#pragma once

#include <cmath>
#include <array>

namespace Math {

	template<class Functor>
	double
	RombergIntegral(const double a, const double b, const double eps, Functor fcn)
	{
		constexpr int maxIterations = 30;
		
		std::array<double, maxIterations> v = {0.0};
		
		double stepSize = b-a;
		v[0] = 0.5*stepSize*(fcn(a) + fcn(b));
		
		double valueBefore = v[0];
		int numberOfSteps = 1;
		
		for (int k = 1; k < maxIterations; k++) {
			numberOfSteps *= 2;
			stepSize /= 2.0;
			
			for (int i = 1; i < numberOfSteps; i+=2)
				v[k] += fcn(a + i*stepSize);
			v[k] = v[k]*stepSize + 0.5*v[k-1];
			
			for (int j = k-1; j >= 0; j--)
				v[j] = v[j+1] + (v[j+1] - v[j]) / (std::pow(4,k-j) - 1.0);
			
			if (std::fabs(valueBefore - v[0]) < eps)
				return v[0];
			else
				valueBefore = v[0];
		}
		
		return -1;
	}
	
}
