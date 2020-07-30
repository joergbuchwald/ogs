/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2020, 9:32 AM
 */

#include "NonisothermalRichardsFlowMechanicsStaggerdTcHcM.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NonIsothermalRichardsFlowMechanicsTcHcMFEM.h"
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
template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<
    DisplacementDim>::constructDofTableImplementation()
{
    // So far, we use the same order element for T,H, and M process
    // Create single component dof in every of the mesh's nodes.
    this->_mesh_subset_all_nodes = std::make_unique<MeshLib::MeshSubset>(
        this->_mesh, this->_mesh.getNodes());

    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *this->_mesh_subset_all_nodes};
    // For T, H processes and also for extrapolation.
    this->local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    sparsity_pattern_T_or_H_ = NumLib::computeSparsityPattern(
        *this->local_to_global_index_map_single_component_, this->_mesh);

    // For displacement equation.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        this->getProcessVariables(data_of_staggeredTcHcM_.M_process_id)[0]
            .get()
            .getNumberOfGlobalComponents(),
        [&]() { return *this->_mesh_subset_all_nodes; });

    std::vector<int> const vec_n_components{DisplacementDim};
    this->_local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    assert(this->_local_to_global_index_map);
    assert(this->local_to_global_index_map_single_component_);
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::
    getMatrixSpecificationsImplementation(const int process_id) const
{
    if (process_id == data_of_staggeredTcHcM_.M_process_id)
    {
        auto const& l = *this->_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For T or H process.
    auto const& l = *this->local_to_global_index_map_single_component_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_T_or_H_};
}

template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<
    DisplacementDim>::initializeBoundaryConditionsImplementation()
{
    this->initializeProcessBoundaryConditionsAndSourceTerms(
        *this->local_to_global_index_map_single_component_,
        data_of_staggeredTcHcM_.T_process_id);
    this->initializeProcessBoundaryConditionsAndSourceTerms(
        *this->local_to_global_index_map_single_component_,
        data_of_staggeredTcHcM_.H_process_id);

    this->initializeProcessBoundaryConditionsAndSourceTerms(
        *this->_local_to_global_index_map,
        data_of_staggeredTcHcM_.M_process_id);
}

template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::
    setInitialConditionsConcreteProcessImplementation(GlobalVector const& /*x*/,
                                                      double const /*t*/)
{
}

template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::
    initializeConcreteProcessImplementation(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order)
{
    using nlohmann::json;

    const int deformation_variable_id = 0;
    const unsigned M_shape_function_orger =
        this->getProcessVariables(
                data_of_staggeredTcHcM_.M_process_id)[deformation_variable_id]
            .get()
            .getShapeFunctionOrder();
    const int pressure_variable_id = 0;
    const unsigned other_shape_function_orger =
        this->getProcessVariables(
                data_of_staggeredTcHcM_.H_process_id)[pressure_variable_id]
            .get()
            .getShapeFunctionOrder();

    ProcessLib::RichardsMechanics::createLocalAssemblers<
        DisplacementDim, NonIsothermalRichardsFlowMechanicsAssemblerTcHcM>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        M_shape_function_orger, other_shape_function_orger,
        this->local_assemblers_, mesh.isAxiallySymmetric(), integration_order,
        this->process_data_, data_of_staggeredTcHcM_);

    this->initializeConcreteProcessData(mesh);
    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerIF::initialize, this->local_assemblers_,
        *this->_local_to_global_index_map);
}

template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::
    assembleConcreteProcessImplementation(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    if (process_id == data_of_staggeredTcHcM_.T_process_id)
    {
        INFO(
            "Assemble the equations for the thermal process (T) of "
            "NonisothermalRichardsFlowMechanics");
    }
    if (process_id == data_of_staggeredTcHcM_.H_process_id)
    {
        INFO(
            "Assemble the equations for the hydraulic process (H) of "
            "NonisothermalRichardsFlowMechanics");
    }

    ProcessLib::ProcessVariable const& pv =
        this->getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        this->_global_assembler, &VectorMatrixAssembler::assemble,
        this->local_assemblers_, pv.getActiveElementIDs(), getDOFTables(), t,
        dt, x, xdot, process_id, M, K, b);
}

template <int DisplacementDim>
void NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::
    assembleWithJacobianConcreteProcessImplementation(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    INFO(
        "Assemble the equations for the mechanical process (M) of "
        "NonisothermalRichardsFlowMechanics");

    ProcessLib::ProcessVariable const& pv =
        this->getProcessVariables(process_id)[0];

    auto const dof_tables = getDOFTables();
    GlobalExecutor::executeSelectedMemberDereferenced(
        this->_global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        this->local_assemblers_, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        xdot, dxdot_dx, dx_dx, process_id, M, K, b, Jac);
}

template <int DisplacementDim>
std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
NonisothermalRichardsFlowMechanicsStaggerdTcHcM<DisplacementDim>::getDOFTables()
    const
{
    // declare a vector with an initialization
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables{{*this->_local_to_global_index_map,
                    *this->_local_to_global_index_map,
                    *this->_local_to_global_index_map}};

    // Assign the elements of the vector with the references to the  desired dof
    // tables, respectively.
    dof_tables[data_of_staggeredTcHcM_.T_process_id] =
        *this->local_to_global_index_map_single_component_;
    dof_tables[data_of_staggeredTcHcM_.H_process_id] =
        *this->local_to_global_index_map_single_component_;
    dof_tables[data_of_staggeredTcHcM_.M_process_id] =
        *this->_local_to_global_index_map;

    return dof_tables;
}

template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<2>;
template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<3>;
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
