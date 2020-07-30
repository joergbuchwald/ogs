/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 20, 2020, 7:54 AM
 */

#pragma once

#include <vector>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
namespace MPL = MaterialPropertyLib;

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename Derived, typename ShapeFunctionDisplacement,
          typename ShapeFunctionPressure, typename IntegrationMethod,
          int DisplacementDim>
class NonIsothermalRichardsFlowMechanicsFEM
    : public RichardsMechanics::RichardsMechanicsLocalAssembler<
          ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
          DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using KelvinVectorType = typename BMatricesType::KelvinVectorType;

    using IpData = RichardsMechanics::IntegrationPointData<
        BMatricesType, ShapeMatricesTypeDisplacement, ShapeMatricesTypePressure,
        DisplacementDim, ShapeFunctionDisplacement::NPOINTS>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    NonIsothermalRichardsFlowMechanicsFEM(
        NonIsothermalRichardsFlowMechanicsFEM const&) = delete;
    NonIsothermalRichardsFlowMechanicsFEM(
        NonIsothermalRichardsFlowMechanicsFEM&&) = delete;

    NonIsothermalRichardsFlowMechanicsFEM(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>&
            process_data)
        : RichardsMechanics::RichardsMechanicsLocalAssembler<
              ShapeFunctionDisplacement, ShapeFunctionPressure,
              IntegrationMethod, DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }

    void assembleForStaggeredScheme(double const t, double const dt,
                                    Eigen::VectorXd const& local_x,
                                    Eigen::VectorXd const& local_xdot,
                                    int const process_id,
                                    std::vector<double>& local_M_data,
                                    std::vector<double>& local_K_data,
                                    std::vector<double>& local_b_data) override
    {
        static_cast<Derived*>(this)->assembleForStaggeredSchemeImplementation(
            t, dt, local_x, local_xdot, process_id, local_M_data, local_K_data,
            local_b_data);
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double dxdot_dx,
        const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override
    {
        static_cast<Derived*>(this)
            ->assembleWithJacobianForStaggeredSchemeImplementation(
                t, dt, local_x, local_xdot, dxdot_dx, dx_dx, process_id,
                local_M_data, local_K_data, local_b_data, local_Jac_data);
    }

    void assembleForStaggeredSchemeImplementation(
        double const /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_xdot*/, int const /*process_id*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/)
    {
    }

    void assembleWithJacobianForStaggeredSchemeImplementation(
        double const /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, int const /*process_id*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
    {
    }
};

}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib

#include "NonIsothermalRichardsFlowMechanicsFEM-impl.h"
