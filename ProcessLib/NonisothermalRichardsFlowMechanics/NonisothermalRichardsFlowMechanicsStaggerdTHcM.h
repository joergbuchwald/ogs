/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2020, 8:57 AM
 */

#pragma once

#include "DataOfStaggeredTHcM.h"
#include "NonisothermalRichardsFlowMechanicsProcess.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
/// Staggered scheme among TH and M processes, where HM is solved with the
/// monolithic scheme.
template <int DisplacementDim>
class NonisothermalRichardsFlowMechanicsStaggerdTHcM final
    : public NonisothermalRichardsFlowMechanicsProcess<
          NonisothermalRichardsFlowMechanicsStaggerdTHcM<DisplacementDim>,
          DisplacementDim>
{
public:
    NonisothermalRichardsFlowMechanicsStaggerdTHcM(
        std::string name, MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>&&
            process_data,
        SecondaryVariableCollection&& secondary_variables,
        DataOfStaggeredTHcM&& data_of_staggeredTHcM)
        : NonisothermalRichardsFlowMechanicsProcess<
              NonisothermalRichardsFlowMechanicsStaggerdTHcM<DisplacementDim>,
              DisplacementDim>(
              std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(process_data), std::move(secondary_variables), false),
          data_of_staggeredTHcM_(std::move(data_of_staggeredTHcM))
    {
    }

private:
    DataOfStaggeredTHcM data_of_staggeredTHcM_;
};

extern template class NonisothermalRichardsFlowMechanicsStaggerdTHcM<2>;
extern template class NonisothermalRichardsFlowMechanicsStaggerdTHcM<3>;

}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
