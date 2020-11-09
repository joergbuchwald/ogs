/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoRichardsFlowProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ThermoRichardsFlowProcess.h"
#include "ThermoRichardsFlowProcessData.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::porosity, MaterialPropertyLib::biot_coefficient,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density,
        MaterialPropertyLib::bulk_modulus};
    std::array const required_solid_properties = {MaterialPropertyLib::density};

    // Thermal properties are not checked because they can be phase property or
    // meduim property (will be enabled later).
    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_solid_properties);
    }
}

void checkProcessVariableComponents(ProcessVariable const& variable,
                                    const int dim)
{
    DBUG("Associate displacement with process variable '{:s}'.",
         variable.getName());

    if (variable.getNumberOfGlobalComponents() != dim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable.getName(),
            variable.getNumberOfGlobalComponents(),
            dim);
    }
}

std::unique_ptr<Process> createThermoRichardsFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_RICHARDS_FLOW");
    DBUG("Create ThermoRichardsFlowProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_T;
    ProcessVariable* variable_p;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_RICHARDS_FLOW__process_variables__pressure}
             "pressure"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL(
            "So far, only the monolithic scheme is implemented for "
            "THERMO_RICHARDS_FLOW");
    }

    checkProcessVariableComponents(*variable_T, 1);
    checkProcessVariableComponents(*variable_p, 1);

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != mesh.getDimension())
        {
            OGS_FATAL(
                "specific body force (gravity vector) has {:d} components, "
                "but mesh dimension is {:d}",
                b.size(), mesh.getDimension());
        }
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG(
        "Check the media properties of ThermoRichardsFlow process "
        "...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    bool mass_lumping = false;
    if (auto const mass_lumping_ptr =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__mass_lumping}
            config.getConfigParameterOptional<bool>("mass_lumping"))
    {
        DBUG("Using mass lumping for the Richards flow equation.");
        mass_lumping = *mass_lumping_ptr;
    }

    bool has_water_vaporization = false;
    if (auto const has_water_vaporization_ptr =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_MECHANICS__has_water_vaporization}
        config.getConfigParameterOptional<bool>("has_water_vaporization"))
    {
        DBUG("Consider water vaporization.");
        has_water_vaporization = *has_water_vaporization_ptr;
    }

    ThermoRichardsFlowProcessData process_data{
        materialIDs(mesh), std::move(media_map), specific_body_force,
        mass_lumping, has_water_vaporization};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<ThermoRichardsFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}
}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
