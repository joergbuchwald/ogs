/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2020, 9:01 AM
 */

#pragma once
#include "NonisothermalRichardsFlowMechanicsProcess.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
/// Monolithic scheme for T, H, and M processes.
template <int DisplacementDim>
class NonisothermalRichardsFlowMechanicsMonolithicTHM final
    : public NonisothermalRichardsFlowMechanicsProcess<
          NonisothermalRichardsFlowMechanicsMonolithicTHM<DisplacementDim>,
          DisplacementDim>
{
public:
    NonisothermalRichardsFlowMechanicsMonolithicTHM(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>&&
            process_data,
        SecondaryVariableCollection&& secondary_variables)
        : NonisothermalRichardsFlowMechanicsProcess<
              NonisothermalRichardsFlowMechanicsMonolithicTHM<DisplacementDim>,
              DisplacementDim>(
              std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(process_data), std::move(secondary_variables), true)
    {
    }
};

extern template class NonisothermalRichardsFlowMechanicsMonolithicTHM<2>;
extern template class NonisothermalRichardsFlowMechanicsMonolithicTHM<3>;

}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
