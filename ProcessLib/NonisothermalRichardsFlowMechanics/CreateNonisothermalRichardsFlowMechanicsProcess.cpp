/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 20, 2020, 12:56 PM
 */

#include "CreateNonisothermalRichardsFlowMechanicsProcess.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "DataOfStaggeredTHcM.h"
#include "DataOfStaggeredTcHcM.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NonisothermalRichardsFlowMechanicsMonolithicTHM.h"
#include "NonisothermalRichardsFlowMechanicsStaggerdTHcM.h"
#include "NonisothermalRichardsFlowMechanicsStaggerdTcHcM.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsProcessData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::porosity, MaterialPropertyLib::biot_coefficient,
        MaterialPropertyLib::bishops_effective_stress,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation, MaterialPropertyLib::tortuosity};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density,
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::thermal_conductivity,
        MaterialPropertyLib::thermal_diffusion_enhancement_factor};
    std::array const required_gas_properties = {
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::thermal_conductivity};
    std::array const required_solid_properties = {
        MaterialPropertyLib::specific_heat_capacity,
        MaterialPropertyLib::thermal_conductivity,
        MaterialPropertyLib::thermal_expansivity, MaterialPropertyLib::density};

    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Gas"),
                                required_gas_properties);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_solid_properties);
    }
}

template <int DisplacementDim>
std::unique_ptr<Process> createNonisothermalRichardsFlowMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "NONISOTHERMAL_RICHARDSFLOW_MECHANICS");
    DBUG("Create NonisothermalRichardsFlowMechanicsProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");

    auto const has_vapor_diffusion =
        //! \ogs_file_param{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__has_vapor_diffusion}
        config.getConfigParameter<bool>("has_vapor_diffusion", false);

    // Process variable.

    //! \ogs_file_param{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_T = nullptr;
    ProcessVariable* variable_p = nullptr;
    ProcessVariable* variable_u = nullptr;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    // monolithic scheme
    if ((staggered_scheme && (*staggered_scheme == "monolithic")) ||
        !staggered_scheme)
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {
             //! \ogs_file_param_special{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        variable_u = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else if ((staggered_scheme && (*staggered_scheme == "staggered_T_H_M")))
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"temperature"s, "pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_T = &process_variables[0][0].get();
        variable_p = &process_variables[1][0].get();
        variable_u = &process_variables[2][0].get();
    }
    else if ((staggered_scheme && (*staggered_scheme == "staggered_TH_M")))
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"temperature"s, "pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_T = &process_variables[0][0].get();
        variable_p = &process_variables[0][1].get();
        variable_u = &process_variables[1][0].get();
    }

    DBUG("Associate temperature with process variable '{:s}'.",
         variable_T->getName());
    if (variable_T->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '{:s}' is not a scalar variable but "
            "has "
            "{:d} components.",
            variable_T->getName(),
            variable_T->getNumberOfGlobalComponents());
    }

    DBUG("Associate pressure with process variable '{:s}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG(
        "Check the media properties of NONISOTHERMAL_RICHARDSFLOW_MECHANICS "
        "process ...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    // Initial stress conditions
    auto const initial_stress = ParameterLib::findOptionalTagParameter<double>(
        //! \ogs_file_param_special{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__initial_stress}
        config, "initial_stress", parameters,
        // Symmetric tensor size, 4 or 6, not a Kelvin vector.
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value,
        &mesh);

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__NONISOTHERMAL_RICHARDSFLOW_MECHANICS__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>
        process_data{materialIDs(mesh),
                     std::move(media_map),
                     std::move(solid_constitutive_relations),
                     initial_stress,
                     specific_body_force,
                     mass_lumping};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    if ((staggered_scheme && (*staggered_scheme == "monolithic")) ||
        !staggered_scheme)
    {
        return std::make_unique<
            NonisothermalRichardsFlowMechanicsMonolithicTHM<DisplacementDim>>(
            std::move(name), mesh, std::move(jacobian_assembler), parameters,
            integration_order, std::move(process_variables),
            std::move(process_data), std::move(secondary_variables));
    }
    else if ((staggered_scheme && (*staggered_scheme == "staggered_T_H_M")))
    {
        DataOfStaggeredTcHcM data_of_staggeredTcHcM{0, 1, 2,
                                                    has_vapor_diffusion};

        return std::make_unique<
            NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>>(
            std::move(name), mesh, std::move(jacobian_assembler), parameters,
            integration_order, std::move(process_variables),
            std::move(process_data), std::move(secondary_variables),
            std::move(data_of_staggeredTcHcM));
    }
    else if ((staggered_scheme && (*staggered_scheme == "staggered_TH_M")))
    {
        DataOfStaggeredTHcM data_of_staggeredTHcM{0, 1, has_vapor_diffusion};
        return std::make_unique<
            NonisothermalRichardsFlowMechanicsStaggerdTHcM<DisplacementDim>>(
            std::move(name), mesh, std::move(jacobian_assembler), parameters,
            integration_order, std::move(process_variables),
            std::move(process_data), std::move(secondary_variables),
            std::move(data_of_staggeredTHcM));
    }
    else
    {
        OGS_FATAL(
            "The coupling scheme of {:s} is not available for "
            "NONISOTHERMAL_RICHARDSFLOW_MECHANICS ",
            *staggered_scheme);
    }

    return {};
}

template std::unique_ptr<Process>
createNonisothermalRichardsFlowMechanicsProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process>
createNonisothermalRichardsFlowMechanicsProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
