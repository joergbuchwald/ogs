/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
struct MechanicsModel
{
    MechanicsModel() { DBUG("Prototype constructor..."); }
    virtual double storageContribution(Eigen::Vector3d /*E*/,
                                       Eigen::Vector3d /*nu*/)
    {
        return 0.0;
    }
    virtual double thermalExpansivityContribution(
        Eigen::
            Matrix<double, 3, 3> /*solid_linear_thermal_expansion_coefficient*/,
        Eigen::Vector3d /*E*/,
        Eigen::Vector3d /*nu*/)
    {
        return 0.0;
    }

    virtual double bulkCompressibilityFromYoungsModulus(Eigen::Vector3d E,
                                                        Eigen::Vector3d nu)
    {
        auto const nu12 = nu[0];
        auto const nu23 = nu[1];
        auto const nu13 = nu[2];
        auto const beta_ortho =
            (E[0] * E[1] * E[2] == 0.0)
                ? 0.0
                : (E[0] * E[1] + E[0] * E[2] * (1 - 2 * nu23) +
                   E[1] * E[2] * (1 - 2 * nu12 - 2 * nu13)) /
                      (E[0] * E[1] * E[2]);
        return beta_ortho;
    }
    virtual bool custommodel()
    {
        return false;
    }

        virtual ~MechanicsModel() = default;
};

struct Uniaxial : MechanicsModel
{
    Uniaxial() { DBUG("using uniaxial simplified mechanics model"); }

    double storageContribution(Eigen::Vector3d E, Eigen::Vector3d nu) override
    {
        auto const nu12 = nu[0];
        auto const nu23 = nu[1];
        auto const nu13 = nu[2];
        auto const nu21 = nu12 * E[1] / E[0];
        auto const nu32 = nu23 * E[2] / E[1];
        auto const nu31 = nu13 * E[2] / E[0];
        auto const D = 1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                       2 * nu12 * nu23 * nu31;
        return D / (E[2] * (1 - nu12 * nu21));
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> solid_linear_thermal_expansion_coefficient,
        Eigen::Vector3d E,
        Eigen::Vector3d nu) override
    {
        auto const nu12 = nu[0];
        auto const nu23 = nu[1];
        auto const nu13 = nu[2];
        auto const nu21 = nu12 * E[1] / E[0];
        auto const D = (1 - nu12 * nu21);
        return -(solid_linear_thermal_expansion_coefficient(2, 2) +
                 solid_linear_thermal_expansion_coefficient(0, 0) *
                     (nu13 + nu12 * nu23) / D +
                 solid_linear_thermal_expansion_coefficient(1, 1) *
                     (nu23 + nu13 * nu21) / D);
    }
    bool custommodel() override
    {
        return false;
    }
};

struct Hydrostatic : MechanicsModel
{
    Hydrostatic() { DBUG("using hydrostatic simplified mechanics model"); }

    double storageContribution(Eigen::Vector3d E, Eigen::Vector3d nu) override
    {
        DBUG("storage contribution uhydrostatic class");
        const double beta_ortho =
            this->bulkCompressibilityFromYoungsModulus(E, nu);
        return beta_ortho;
    }

    double thermalExpansivityContribution(
        Eigen::Matrix<double, 3, 3> solid_linear_thermal_expansion_coefficient,
        Eigen::Vector3d /*E*/,
        Eigen::Vector3d /*nu*/) override
    {
        return -solid_linear_thermal_expansion_coefficient.trace();
    }
    bool custommodel() override
    {
        return false;
    }
};

struct Rigid : MechanicsModel
{
    Rigid() { DBUG("using rigid mechanics model"); }

    double storageContribution(Eigen::Vector3d /*E*/,
                               Eigen::Vector3d /*nu*/) override
    {
        return 0.0;
    }

    double thermalExpansivityContribution(
        Eigen::
            Matrix<double, 3, 3> /*solid_linear_thermal_expansion_coefficient*/,
        Eigen::Vector3d /*E*/,
        Eigen::Vector3d /*nu*/) override
    {
        return 0.0;
    }
    bool custommodel() override
    {
        return false;
    }
};

struct Custom : MechanicsModel
{
    Custom() { DBUG("using custom mechanics model"); }

    double storageContribution(Eigen::Vector3d /*E*/,
                               Eigen::Vector3d /*nu*/) override
    {
        OGS_FATAL("bitte hier nicht");
        return 0.0;
    }

    double thermalExpansivityContribution(
        Eigen::
            Matrix<double, 3, 3> /*solid_linear_thermal_expansion_coefficient*/,
        Eigen::Vector3d /*E*/,
        Eigen::Vector3d /*nu*/) override
    {
        OGS_FATAL("bitte hier nicht");
        return 0.0;
    }
    bool custommodel() override
    {
        return true;
    }
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
