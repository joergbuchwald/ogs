/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 10:56 AM
 */

#pragma once

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
/// \f$\rho_{vS}\f$
double saturatedVaporDensity(double const T);

/// \f$\frac{\partial \rho_{vS}}{\partial T}\f$
double saturatedVaporDensityDerivative(double const T);

double humidity(double const T, double const p, double const water_density);

double waterVaporDensity(double const T, double const p,
                         double const water_density);

double dwaterVaporDensitydT(double const T, double const p,
                            double const water_density);

/// \f$D_v \f$
double waterVaporDiffusion(double const T, double const Dvr);

/// \f$ D_{vr} \f$. For the FEBEX bentonite
double relativeDiffusionCoefficient(double const S, double const porosity,
                                    double const tortuosity);

double Dpv(double const T, double const p, double const water_density,
           double const S, double const porosity, double const tortuosity);

double DTv(double const T, double const p, double const water_density,
           double const S, double const porosity, double const tortuosity);
}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib

