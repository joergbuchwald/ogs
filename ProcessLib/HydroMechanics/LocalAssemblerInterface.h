/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    virtual std::size_t setIPDataInitialConditions(
        std::string const& name, double const* values,
        int const integration_order) = 0;

    virtual std::vector<double> getSigma() const = 0;

    virtual std::vector<double> getEpsilon() const = 0;

    virtual std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    // TODO move to NumLib::ExtrapolatableElement
    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned /*integration_point*/) const = 0;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
