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

#pragma once

#include "ProcessLib/Process.h"
#include "ProcessLib/RichardsMechanics/LocalAssemblerInterface.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
template <typename Derived, int DisplacementDim>
class NonisothermalRichardsFlowMechanicsProcess : public Process
{
public:
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
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override;
    //! @}

    /**
     * Get the size and the sparse pattern of the global matrix in order to
     * create the global matrices and vectors for the system equations of this
     * process.
     *
     * @param process_id Process ID.
     * @return Matrix specifications including size and sparse pattern.
     */

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override
    {
        return static_cast<Derived const*>(this)
            ->getMatrixSpecificationsImplementation(process_id);
    }

public:  // CRTP implementation
    bool hasMechanicalProcessImplementation(int const /*process_id*/) const
    {
        return true;
    }

    NumLib::LocalToGlobalIndexMap& getDOFTableImplementation(
        const int /*process_id*/) const
    {
        return *_local_to_global_index_map;
    }

    void constructDofTableImplementation() {}
    void initializeBoundaryConditionsImplementation() {}
    MathLib::MatrixSpecifications getMatrixSpecificationsImplementation(
        const int /*process_id*/) const
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    void setInitialConditionsConcreteProcessImplementation(
        GlobalVector const& /*x*/, double const /*t*/)
    {
    }
    void initializeConcreteProcessImplementation(
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        MeshLib::Mesh const& /*mesh*/,
        unsigned const /*integration_order*/)
    {
    }

    void assembleConcreteProcessImplementation(
        const double /*t*/, double const /*dt*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<GlobalVector*> const& /*xdot*/, int const /*process_id*/,
        GlobalMatrix& /*M*/, GlobalMatrix& /*K*/, GlobalVector& /*b*/)
    {
    }

    void assembleWithJacobianConcreteProcessImplementation(
        const double /*t*/, double const /*dt*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<GlobalVector*> const& /*xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, int const /*process_id*/, GlobalMatrix& /*M*/,
        GlobalMatrix& /*K*/, GlobalVector& /*b*/, GlobalMatrix& /*Jac*/)
    {
    }

    void computeSecondaryVariableImplementation(
        double const /*t*/, double const /*dt*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<GlobalVector*> const& /*xdot*/, int const /*process_id*/)
    {
    }

protected:
    using LocalAssemblerIF =
        RichardsMechanics::LocalAssemblerInterface<DisplacementDim>;

    RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>
        process_data_;

    std::vector<std::unique_ptr<LocalAssemblerIF>> local_assemblers_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        local_to_global_index_map_single_component_;

    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 3> xs_previous_timestep_;
    MeshLib::PropertyVector<double>* nodal_forces_ = nullptr;
    MeshLib::PropertyVector<double>* hydraulic_flow_ = nullptr;

    void constructDofTable() override
    {
        static_cast<Derived*>(this)->constructDofTableImplementation();
    }

    void initializeConcreteProcessData(MeshLib::Mesh const& mesh);

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        static_cast<Derived*>(this)->initializeConcreteProcessImplementation(
            dof_table, mesh, integration_order);
    }

    void initializeBoundaryConditions() override
    {
        static_cast<Derived*>(this)
            ->initializeBoundaryConditionsImplementation();
    }

    void setInitialConditionsConcreteProcess(GlobalVector const& x,
                                             double const t) override
    {
        static_cast<Derived*>(this)
            ->setInitialConditionsConcreteProcessImplementation(x, t);
    }

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override
    {
        static_cast<Derived*>(this)->assembleConcreteProcessImplementation(
            t, dt, x, xdot, process_id, M, K, b);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override
    {
        static_cast<Derived*>(this)
            ->assembleWithJacobianConcreteProcessImplementation(
                t, dt, x, xdot, dxdot_dx, dx_dx, process_id, M, K, b, Jac);
    }

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     double const t, double const dt,
                                     const int process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override
    {
        return static_cast<Derived const*>(this)->getDOFTableImplementation(
            process_id);
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, int const process_id) override
    {
        static_cast<Derived*>(this)->computeSecondaryVariableImplementation(
            t, dt, x, xdot, process_id);
    }

    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const override;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process.
    bool hasMechanicalProcess(int const process_id) const
    {
        return static_cast<Derived const*>(this)
            ->hasMechanicalProcessImplementation(process_id);
    }
};
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
