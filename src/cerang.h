#pragma once

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
