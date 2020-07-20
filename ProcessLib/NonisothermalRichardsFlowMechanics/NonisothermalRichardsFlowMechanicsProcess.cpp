/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 20, 2020, 9:13 AM
 */

#include "NonisothermalRichardsFlowMechanicsProcess.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
// TEST #include "NonisothermalRichardsFlowMechanicsFEM.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/RichardsMechanics/CreateLocalAssemblers.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
template <typename Derived, int DisplacementDim>
NonisothermalRichardsFlowMechanicsProcess<Derived, DisplacementDim>::
    NonisothermalRichardsFlowMechanicsProcess(
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
        bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      process_data_(std::move(process_data))
{
}

template <typename Derived, int DisplacementDim>
bool NonisothermalRichardsFlowMechanicsProcess<
    Derived, DisplacementDim>::isLinear() const
{
    return false;
}

template <typename Derived, int DisplacementDim>
void NonisothermalRichardsFlowMechanicsProcess<Derived, DisplacementDim>::
    initializeConcreteProcessData(MeshLib::Mesh const& mesh)
{
    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function) {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             local_assemblers_,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtSigma);

    add_secondary_variable("swelling_stress",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtSwellingStress);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerIF::getIntPtEpsilon);

    add_secondary_variable("velocity", DisplacementDim,
                           &LocalAssemblerIF::getIntPtDarcyVelocity);

    add_secondary_variable("saturation", 1,
                           &LocalAssemblerIF::getIntPtSaturation);

    add_secondary_variable("porosity", 1, &LocalAssemblerIF::getIntPtPorosity);

    //
    // enable output of internal variables defined by material models
    //
    ProcessLib::Deformation::solidMaterialInternalToSecondaryVariables<
        LocalAssemblerIF>(process_data_.solid_materials,
                          add_secondary_variable);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& mesh_property =
            *mesh.getProperties().template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data = getIntegrationPointMetaData(mesh, name);

        // Check the number of components.
        if (ip_meta_data.n_components !=
            mesh_property.getNumberOfGlobalComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data ({:d}) than in "
                "the integration point field data for '{:s}': {:d}.",
                ip_meta_data.n_components, name,
                mesh_property.getNumberOfGlobalComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : local_assemblers_)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &mesh_property[position],
                    ip_meta_data.integration_order);
            if (integration_points_read == 0)
            {
                OGS_FATAL(
                    "No integration points read in the integration point "
                    "initial conditions set function.");
            }
            position += integration_points_read * ip_meta_data.n_components;
        }
    }
}

template <typename Derived, int DisplacementDim>
void NonisothermalRichardsFlowMechanicsProcess<Derived, DisplacementDim>::
    postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                double const t, double const dt,
                                const int process_id)
{
    DBUG("PostTimestep  NonisothermalRichardsFlowMechanicsProcess.");

    if (hasMechanicalProcess(process_id))
    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerIF::postTimestep, local_assemblers_,
            pv.getActiveElementIDs(), *_local_to_global_index_map,
            *x[process_id], t, dt);
    }
}

template <typename Derived, int DisplacementDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
NonisothermalRichardsFlowMechanicsProcess<
    Derived, DisplacementDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(local_to_global_index_map_single_component_.get(),
                           manage_storage);
}
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
