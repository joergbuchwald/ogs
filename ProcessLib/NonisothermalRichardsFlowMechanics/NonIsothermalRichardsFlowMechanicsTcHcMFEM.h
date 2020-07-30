/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 30, 2020, 1:15 PM
 */

#pragma once

#include "DataOfStaggeredTcHcM.h"
#include "NonIsothermalRichardsFlowMechanicsFEM.h"
#include "ProcessLib/RichardsMechanics/IntegrationPointData.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
class NonIsothermalRichardsFlowMechanicsAssemblerTcHcM
    : public NonIsothermalRichardsFlowMechanicsFEM<
          NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
              ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
              DisplacementDim>,
          ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
          DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using KelvinVectorType = typename BMatricesType::KelvinVectorType;

    using IpData = RichardsMechanics::IntegrationPointData<
        BMatricesType, ShapeMatricesTypeDisplacement, ShapeMatricesType,
        DisplacementDim, ShapeFunctionDisplacement::NPOINTS>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    NonIsothermalRichardsFlowMechanicsAssemblerTcHcM(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        RichardsMechanics::RichardsMechanicsProcessData<DisplacementDim>&
            process_data,
        DataOfStaggeredTcHcM& data_of_staggeredTcHcM)
        : NonIsothermalRichardsFlowMechanicsFEM<
              NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
                  ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
                  DisplacementDim>,
              ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
              DisplacementDim>(e, local_matrix_size, is_axially_symmetric,
                               integration_order, process_data),
          data_of_staggeredTcHcM_(data_of_staggeredTcHcM)
    {
    }

    void assembleForStaggeredSchemeImplementation(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data);

    void assembleWithJacobianForStaggeredSchemeImplementation(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double dxdot_dx,
        const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

private:
    DataOfStaggeredTcHcM& data_of_staggeredTcHcM_;

    void assembleForEnergyBalanceEquations(double const t, double const dt,
                                           Eigen::VectorXd const& local_x,
                                           int const process_id,
                                           std::vector<double>& local_M_data,
                                           std::vector<double>& local_K_data,
                                           std::vector<double>& local_b_data);

    void assembleForMassBalanceEquations(double const t, double const dt,
                                         Eigen::VectorXd const& local_x,
                                         Eigen::VectorXd const& local_xdot,
                                         int const process_id,
                                         std::vector<double>& local_M_data,
                                         std::vector<double>& local_K_data,
                                         std::vector<double>& local_b_data);

    void assembleWithJacobianForMomentumEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    int getVariableIndex(const int process_id) const;

    using TemperatureVectorTye =
        typename ShapeMatricesType::template VectorType<ShapeFunction::NPOINTS>;
    using PressureVectorTye =
        typename ShapeMatricesType::template VectorType<ShapeFunction::NPOINTS>;
    using DisplacementVectorTye =
        typename ShapeMatricesTypeDisplacement::template VectorType<
            ShapeFunctionDisplacement::NPOINTS * DisplacementDim>;

    using PressureMatrixTye =
        typename ShapeMatricesType::template MatrixType<ShapeFunction::NPOINTS,
                                                        ShapeFunction::NPOINTS>;
    using TemperatureMatrixTye =
        typename ShapeMatricesType::template MatrixType<ShapeFunction::NPOINTS,
                                                        ShapeFunction::NPOINTS>;
    using DisplacementMatrixTye =
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            ShapeFunctionDisplacement::NPOINTS * DisplacementDim,
            ShapeFunctionDisplacement::NPOINTS * DisplacementDim>;
};
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib

#include "NonIsothermalRichardsFlowMechanicsTcHcMFEM-impl.h"
