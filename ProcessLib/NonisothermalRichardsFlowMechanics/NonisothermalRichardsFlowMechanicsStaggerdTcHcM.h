/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2020, 8:34 AM
 */

#pragma once

#include "DataOfStaggeredTcHcM.h"
#include "NonisothermalRichardsFlowMechanicsProcess.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
/// Staggered scheme among T, H, and M processes.
template <int DisplacementDim>
class NonisothermalRichardsFlowMechanicsStaggerdTcHcM final
    : public NonisothermalRichardsFlowMechanicsProcess<
          NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>,
          DisplacementDim>
{
public:
    NonisothermalRichardsFlowMechanicsStaggerdTcHcM(
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
        SecondaryVariableCollection&& secondary_variables,
        DataOfStaggeredTcHcM&& data_of_staggeredTcHcM)
        : NonisothermalRichardsFlowMechanicsProcess<
              NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>,
              DisplacementDim>(
              std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(process_data), std::move(secondary_variables), false),
          data_of_staggeredTcHcM_(std::move(data_of_staggeredTcHcM))
    {
    }

private:
    DataOfStaggeredTcHcM data_of_staggeredTcHcM_;
};

extern template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<2>;
extern template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<3>;

}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
