/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LocalAssemblerInterface.h"
#include "ProcessLib/Process.h"
#include "ThermoRichardsFlowProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
/// Global assembler for the monolithic scheme of the non-isothermal Richards
/// flow coupled with mechanics.
///
class ThermoRichardsFlowProcess final : public Process
{
public:
    ThermoRichardsFlowProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoRichardsFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    /**
     * Get the size and the sparse pattern of the global matrix in order to
     * create the global matrices and vectors for the system equations of this
     * process.
     *
     * @param process_id Process ID. If the monolithic scheme is applied,
     *                               process_id = 0. For the staggered scheme,
     *                               process_id = 0 represents the
     *                               hydraulic (H) process, while process_id = 1
     *                               represents the mechanical (M) process.
     * @return Matrix specifications including size and sparse pattern.
     */
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

private:
    using LocalAssemblerIF = LocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void initializeBoundaryConditions() override;

    void setInitialConditionsConcreteProcess(GlobalVector const& x,
                                             double const t,
                                             int const /*process_id*/) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     double const t, double const dt,
                                     const int process_id) override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            GlobalVector const& xdot,
                                            const double t, double const dt,
                                            int const process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    ThermoRichardsFlowProcessData _process_data;

    std::vector<std::unique_ptr<LocalAssemblerIF>> _local_assemblers;

    /// Sparsity pattern for the flow equation, and it is initialized only if
    /// the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_linear_element;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          GlobalVector const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

    MeshLib::PropertyVector<double>* _heat_flux = nullptr;
    MeshLib::PropertyVector<double>* _hydraulic_flow = nullptr;
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
