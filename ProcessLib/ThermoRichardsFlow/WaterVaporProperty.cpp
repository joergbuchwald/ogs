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

#include "WaterVaporProperty.h"

#include <cmath>

#include "MaterialLib/PhysicalConstant.h"
namespace ProcessLib
{
namespace ThermoRichardsFlow
{
/// \f$\rho_{vS}\f$
double saturatedVaporDensity(double const T)
{
    return 1.0e-3 * exp(19.81 - 4975.9 / T);
}

/// \f$\frac{\partial \rho_{vS}}{\partial T}\f$
double saturatedVaporDensityDerivative(double const T)
{
    return 4.9759 * exp(19.81 - 4975.9 / T) / (T * T);
}

double humidity(double const T, double const p, double const water_density)
{
    return std::exp(
        p / (MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
             T * water_density));
}

double waterVaporDensity(double const T, double const p,
                         double const water_density)
{
    return humidity(T, p, water_density) * saturatedVaporDensity(T);
}

double dwaterVaporDensitydT(double const T, double const p,
                            double const water_density)
{
    double const h = humidity(T, p, water_density);
    double const rho_v = h * saturatedVaporDensity(T);
    double const drho_vS_dT = saturatedVaporDensityDerivative(T);

    return h * drho_vS_dT - rho_v * p /
                                (water_density * T * T *
                                 MaterialLib::PhysicalConstant::
                                     SpecificGasConstant::WaterVapour);
}

/// \f$D_v \f$
double waterVaporDiffusion(double const T, double const Dvr)
{
    double base_heat_diffusion_coefficient = 2.16e-5;

    return base_heat_diffusion_coefficient *
           std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                    1.8) *
           Dvr;
}

/// \f$ D_{vr} \f$. For the FEBEX bentonite
double relativeDiffusionCoefficient(double const S, double const porosity,
                                    double const tortuosity)
{
    return tortuosity * porosity * (1 - S);
}

// The Penman-Millington-Quirk (PMQ) model
#ifdef PMQ
double relativeDiffusionCoefficient(double const S, double const porosity,
                                    double const tortuosity)
{
    return tortuosity * porosity * (1 - S) * (1 - S);
}
#endif

double Dpv(double const T, double const p, double const water_density,
           double const S, double const porosity, double const tortuosity)
{
    double const Dvr = relativeDiffusionCoefficient(S, porosity, tortuosity);
    double const Dv = waterVaporDiffusion(T, Dvr);
    double const rho_v = waterVaporDensity(T, p, water_density);
    return Dv * rho_v /
           (water_density *
            MaterialLib::PhysicalConstant::SpecificGasConstant::WaterVapour *
            T);
}

double DTv(double const T, double const p, double const water_density,
           double const S, double const porosity, double const tortuosity)
{
    double const Dvr = relativeDiffusionCoefficient(S, porosity, tortuosity);
    double const Dv = waterVaporDiffusion(T, Dvr);
    return Dv * dwaterVaporDensitydT(T, p, water_density);
}
}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
