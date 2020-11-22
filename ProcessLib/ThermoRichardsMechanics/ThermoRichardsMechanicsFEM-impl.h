/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include <cassert>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ThermoRichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                      ShapeFunctionPressure, IntegrationMethod,
                                      DisplacementDim>::
    ThermoRichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoRichardsMechanicsProcessData<DisplacementDim>& process_data)
    : process_data_(process_data),
      integration_method_(integration_order),
      element_(e),
      is_axially_symmetric_(is_axially_symmetric)
{
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    ip_data_.reserve(n_integration_points);
    secondary_data_.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method_);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, integration_method_);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            process_data_.solid_materials, process_data_.material_ids,
            e.getID());

    auto const& medium = process_data_.media_map->getMedium(element_.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        ip_data_.emplace_back(solid_material);
        auto& ip_data = ip_data_[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data_[ip].integration_weight =
            integration_method_.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
        {
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;
        }

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        // Initial porosity. Could be read from intergration point data or mesh.
        ip_data.porosity =
            medium->property(MPL::porosity)
                .template initialValue<double>(
                    x_position,
                    std::numeric_limits<
                        double>::quiet_NaN() /* t independent */);

        ip_data.transport_porosity = ip_data.porosity;
        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            ip_data.transport_porosity =
                medium->property(MPL::transport_porosity)
                    .template initialValue<double>(
                        x_position,
                        std::numeric_limits<
                            double>::quiet_NaN() /* t independent */);
        }

        secondary_data_.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(integration_method_.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            element_.getID());
    }

    if (name == "sigma_ip")
    {
        if (process_data_.initial_stress != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                process_data_.initial_stress->name);
        }
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::sigma_eff);
    }

    if (name == "saturation_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, ip_data_,
                                                         &IpData::saturation);
    }
    if (name == "porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, ip_data_,
                                                         &IpData::porosity);
    }
    if (name == "transport_porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(
            values, ip_data_, &IpData::transport_porosity);
    }
    if (name == "swelling_stress_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::sigma_sw);
    }
    if (name == "epsilon_ip")
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, ip_data_, &IpData::eps);
    }
    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunctionPressure,
                                           IntegrationMethod, DisplacementDim>::
    setInitialConditionsConcrete(std::vector<double> const& local_x,
                                 double const t,
                                 bool const /*use_monolithic_scheme*/,
                                 int const /*process_id*/)
{
    assert(local_x.size() ==
           temperature_size + pressure_size + displacement_size);

    auto p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    MPL::VariableArray variables;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& N_p = ip_data_[ip].N_p;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        // Note: temperature dependent saturation model is not considered so
        // far.
        ip_data_[ip].saturation_prev =
            medium->property(MPL::PropertyType::saturation)
                .template value<double>(
                    variables, x_position, t,
                    std::numeric_limits<double>::quiet_NaN());
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunctionPressure,
                                           IntegrationMethod, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto const local_matrix_dim =
        displacement_size + pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto const T =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);
    auto const p_L =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    auto const u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto const T_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto const p_L_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_xdot.data() + pressure_index,
                                  pressure_size);
    auto const u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            local_matrix_dim, local_matrix_dim>>(
        local_Jac_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<local_matrix_dim>>(
            local_rhs_data, local_matrix_dim);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    typename ShapeMatricesTypePressure::NodalMatrixType M_TT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         temperature_size);
    typename ShapeMatricesTypePressure::NodalMatrixType K_TT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         temperature_size);
    typename ShapeMatricesTypePressure::NodalMatrixType K_Tp =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(temperature_size,
                                                         pressure_size);
    typename ShapeMatricesTypePressure::NodalMatrixType M_pT =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         temperature_size);
    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = ip_data_[ip].integration_weight;

        auto const& N_u_op = ip_data_[ip].N_u_op;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const& N_p = ip_data_[ip].N_p;
        auto const& dNdx_p = ip_data_[ip].dNdx_p;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, is_axially_symmetric_);

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N_p, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N_p, T_dot_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N_p, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& eps = ip_data_[ip].eps;
        auto const& sigma_eff = ip_data_[ip].sigma_eff;
        auto& S_L = ip_data_[ip].saturation;
        auto const S_L_prev = ip_data_[ip].saturation_prev;
        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            t, x_position, dt, T_ip);

        auto const beta_SR =
            (1 - alpha) /
            ip_data_[ip].solid_material.getBulkModulus(t, x_position, &C_el);
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto const K_LR =
            liquid_phase.property(MPL::PropertyType::bulk_modulus)
                .template value<double>(variables, x_position, t, dt);

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        auto const& b = process_data_.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);

        auto const chi = [medium, x_position, t, dt](double const S_L) {
            MPL::VariableArray variables;
            variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(variables, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        // Set volumetric strain rate for the general case without swelling.
        variables[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(identity2.transpose() * B * u);
        variables_prev[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(identity2.transpose() * B * (u - u_dot * dt));

        auto& phi = ip_data_[ip].porosity;
        {  // Porosity update

            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                ip_data_[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        // Swelling and possibly volumetric strain rate update.
        auto& sigma_sw = ip_data_[ip].sigma_sw;
        {
            auto const& sigma_sw_prev = ip_data_[ip].sigma_sw_prev;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;
            if (solid_phase.hasProperty(
                    MPL::PropertyType::swelling_stress_rate))
            {
                using DimMatrix = Eigen::Matrix<double, 3, 3>;
                auto const sigma_sw_dot =
                    MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                        solid_phase
                            .property(MPL::PropertyType::swelling_stress_rate)
                            .template value<DimMatrix>(
                                variables, variables_prev, x_position, t, dt));
                sigma_sw += sigma_sw_dot * dt;

                // !!! Misusing volumetric strain for mechanical volumetric
                // strain just to update the transport porosity !!!
                std::get<double>(variables[static_cast<int>(
                    MPL::Variable::volumetric_strain)]) +=
                    identity2.transpose() * C_el.inverse() * sigma_sw;
                std::get<double>(variables_prev[static_cast<int>(
                    MPL::Variable::volumetric_strain)]) +=
                    identity2.transpose() * C_el.inverse() * sigma_sw_prev;
            }

            if (solid_phase.hasProperty(MPL::PropertyType::transport_porosity))
            {
                variables_prev[static_cast<int>(
                    MPL::Variable::transport_porosity)] =
                    ip_data_[ip].transport_porosity_prev;

                ip_data_[ip].transport_porosity =
                    solid_phase.property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                    ip_data_[ip].transport_porosity;
            }
            else
            {
                variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                    phi;
            }
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        // For stress dependent permeability.
        variables[static_cast<int>(MPL::Variable::stress)]
            .emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    (ip_data_[ip].sigma_eff +
                     alpha * chi_S_L * identity2 * p_cap_ip)
                        .eval()));
        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const Ki_over_mu = K_intrinsic / mu;
        GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_linear_thermal_expansion_coefficient) *
                T_dot_ip * dt;

        auto C = ip_data_[ip].updateConstitutiveRelationThermal(
            t, x_position, dt, u, T_ip, dthermal_strain);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
            p_FR - (sigma_eff + sigma_sw).dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -= (B.transpose() * (sigma_eff + sigma_sw) -
                           N_u_op.transpose() * rho * b) *
                          w;

        //
        // displacement equation, pressure part
        //
        Kup.noalias() += B.transpose() * alpha * chi_S_L * identity2 * N_p * w;

        auto const dchi_dS_L =
            medium->property(MPL::PropertyType::bishops_effective_stress)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() -= B.transpose() * alpha *
                          (chi_S_L + dchi_dS_L * p_cap_ip * dS_L_dp_cap) *
                          identity2 * N_p * w;

        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op.transpose() * phi * rho_LR * dS_L_dp_cap * b * N_p * w;

        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const dsigma_sw_dS_L =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template dValue<DimMatrix>(
                            variables, variables_prev,
                            MPL::Variable::liquid_saturation, x_position, t,
                            dt));
            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() +=
                B.transpose() * dsigma_sw_dS_L * dS_L_dp_cap * N_p * w;
        }
        //
        // pressure equation, displacement part.
        //
        Kpu.noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                         identity2.transpose() * B * w;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() +=
            dNdx_p.transpose() * k_rel * rho_Ki_over_mu * dNdx_p * w;

        const double alphaB_minu_phi = (alpha > phi) ? alpha - phi : 0.0;
        double const a0 = alphaB_minu_phi * beta_SR;
        double const specific_storage_a_p = S_L * (phi / K_LR + S_L * a0);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi / K_LR + 2 * S_L * a0);
        double const dspecific_storage_a_S_dp_cap =
            -a0 * (S_L + p_cap_ip * dS_L_dp_cap);

        storage_p_a_p.noalias() +=
            N_p.transpose() * rho_LR * specific_storage_a_p * N_p * w;
        if (p_cap_dot_ip != 0)  // prevent division by zero.
        {
            storage_p_a_S.noalias() -= N_p.transpose() * rho_LR *
                                       specific_storage_a_S * (S_L - S_L_prev) /
                                       (dt * p_cap_dot_ip) * N_p * w;
        }

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N_p.transpose() * p_cap_dot_ip * rho_LR *
                          dspecific_storage_a_p_dp_cap * N_p * w;

        storage_p_a_S_Jpp.noalias() -=
            N_p.transpose() * rho_LR *
            ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
             specific_storage_a_S * dS_L_dp_cap) /
            dt * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() -= N_p.transpose() * rho_LR * dS_L_dp_cap * alpha *
                          identity2.transpose() * B * u_dot * N_p * w;

        double const dk_rel_dS_l =
            medium->property(MPL::PropertyType::relative_permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        GlobalDimVectorType const grad_p_cap = -dNdx_p * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_Ki_over_mu * grad_p_cap *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * rho_Ki_over_mu * b *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        //
        // pressure equation, temperature part.
        //
        // Note:  d (rho_l * beta T_)/dp * dotT
        // Add the thermal expansion term is the point is in the fully
        // saturated zone.
        if (p_cap_ip <= 0.0)  // p_l>0.0
        {
            double const fluid_volumetric_thermal_expansion_coefficient =
                MPL::getLiquidThermalExpansivity(liquid_phase, variables,
                                                 rho_LR, x_position, t, dt);
            const double eff_thermal_expansion =
                alphaB_minu_phi *
                    solid_linear_thermal_expansion_coefficient.trace() +
                phi * fluid_volumetric_thermal_expansion_coefficient;

            M_pT.noalias() -=
                N_p.transpose() * rho_LR * eff_thermal_expansion * N_p * w;
        }

        //
        // temperature equation.
        //
        {
            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(variables, x_position, t, dt);

            auto const specific_heat_capacity_solid =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(variables, x_position, t, dt);

            M_TT.noalias() +=
                w *
                (rho_SR * specific_heat_capacity_solid * (1 - phi) +
                 (S_L * rho_LR * specific_heat_capacity_fluid) * phi) *
                N_p.transpose() * N_p;

            auto const thermal_conductivity_solid =
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(variables, x_position, t, dt);

            auto const thermal_conductivity_fluid =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(variables, x_position, t, dt) *
                S_L;

            GlobalDimMatrixType const thermal_conductivity =
                MaterialPropertyLib::formEffectiveThermalConductivity<
                    DisplacementDim>(thermal_conductivity_solid,
                                     thermal_conductivity_fluid, phi);

            GlobalDimVectorType const velocity_l =
                GlobalDimVectorType(-Ki_over_mu * (dNdx_p * p_L - rho_LR * b));

            K_TT.noalias() +=
                (dNdx_p.transpose() * thermal_conductivity * dNdx_p +
                 N_p.transpose() * velocity_l.transpose() * dNdx_p * rho_LR *
                     specific_heat_capacity_fluid) *
                w;

            //
            // temperature equation, pressure part
            //
            K_Tp.noalias() -= rho_LR * specific_heat_capacity_fluid *
                              N_p.transpose() * (dNdx_p * T).transpose() *
                              Ki_over_mu * dNdx_p * w;
        }
    }

    if (process_data_.apply_mass_lumping)
    {
        storage_p_a_p = storage_p_a_p.colwise().sum().eval().asDiagonal();
        storage_p_a_S = storage_p_a_S.colwise().sum().eval().asDiagonal();
        storage_p_a_S_Jpp =
            storage_p_a_S_Jpp.colwise().sum().eval().asDiagonal();
    }

    //
    // -- Jacobian
    //
    // temperature equation.
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += M_TT / dt + K_TT;
    // temperature equation, pressure part
    local_Jac
        .template block<temperature_size, pressure_size>(temperature_index,
                                                         pressure_index)
        .noalias() += K_Tp;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p_a_p / dt + storage_p_a_S_Jpp;

    // pressure equation, temperature part (contributed by thermal expansion).
    local_Jac
        .template block<pressure_size, temperature_size>(pressure_index,
                                                         temperature_index)
        .noalias() += M_pT / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    //
    // -- Residual
    //
    // temperature equation
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        M_TT * T_dot + K_TT * T;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L + (storage_p_a_p + storage_p_a_S) * p_L_dot +
        Kpu * u_dot + M_pT * T_dot;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;
}

template <int Components, typename StoreValuesFunction>
std::vector<double> transposeInPlace(
    StoreValuesFunction const& store_values_function)
{
    std::vector<double> result;
    store_values_function(result);
    // Transpose. For time being Eigen's transposeInPlace doesn't work for
    // non-square mapped matrices.
    MathLib::toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Components, Eigen::RowMajor>>(
        result, result.size() / Components, Components) =
        MathLib::toMatrix<
            Eigen::Matrix<double, Components, Eigen::Dynamic, Eigen::RowMajor>>(
            result, Components, result.size() / Components)
            .transpose()
            .eval();

    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getSigma() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values) {
            return getIntPtSigma(0, {}, {}, values);
        });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        ip_data_, &IpData::sigma_eff, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getSwellingStress() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values) {
            return getIntPtSwellingStress(0, {}, {}, values);
        });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSwellingStress(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    auto const n_integration_points = ip_data_.size();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& sigma_sw = ip_data_[ip].sigma_sw;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_sw);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getEpsilon() const
{
    constexpr int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;

    return transposeInPlace<kelvin_vector_size>(
        [this](std::vector<double>& values) {
            return getIntPtEpsilon(0, {}, {}, values);
        });
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
        ip_data_, &IpData::eps, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = ip_data_[ip].v_darcy;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getSaturation() const
{
    std::vector<double> result;
    getIntPtSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::saturation, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getPorosity() const
{
    std::vector<double> result;
    getIntPtPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(ip_data_,
                                                     &IpData::porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getTransportPorosity() const
{
    std::vector<double> result;
    getIntPtTransportPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtTransportPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::transport_porosity, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        ip_data_, &IpData::dry_density_solid, cache);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunctionPressure,
                                           IntegrationMethod, DisplacementDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                std::vector<double> const& local_xdot,
                                double const t, double const dt,
                                bool const use_monolithic_scheme,
                                int const /*process_id*/)
{
    const int displacement_offset =
        use_monolithic_scheme ? displacement_index : 0;

    auto const T =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);
    auto const T_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

    auto const u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
                                      displacement_size);
    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    MPL::VariableArray variables;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());
    auto const& solid_phase = medium->phase("Solid");

    int const n_integration_points = integration_method_.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = ip_data_[ip].N_p;
        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;
        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N_p, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N_p, T_dot_ip);
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, is_axially_symmetric_);

        auto& eps = ip_data_[ip].eps;
        eps.noalias() = B * u;
        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_linear_thermal_expansion_coefficient) *
                T_dot_ip * dt;

        ip_data_[ip].updateConstitutiveRelationThermal(t, x_position, dt, u,
                                                       T_ip, dthermal_strain);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunctionPressure,
                                           IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_dot)
{
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);
    auto const p_L = local_x.template segment<pressure_size>(pressure_index);
    auto const u =
        local_x.template segment<displacement_size>(displacement_index);

    auto const T_dot =
        local_x_dot.template segment<temperature_size>(temperature_index);
    auto const p_L_dot =
        local_x_dot.template segment<pressure_size>(pressure_index);
    auto const u_dot =
        local_x_dot.template segment<displacement_size>(displacement_index);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    auto const& medium = process_data_.media_map->getMedium(element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(element_.getID());

    unsigned const n_integration_points =
        integration_method_.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = ip_data_[ip].N_p;
        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, is_axially_symmetric_);

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N_p, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N_p, T_dot_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_dot_ip;
        NumLib::shapeFunctionInterpolate(-p_L_dot, N_p, p_cap_dot_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        auto& eps = ip_data_[ip].eps;
        auto& S_L = ip_data_[ip].saturation;
        auto const S_L_prev = ip_data_[ip].saturation_prev;
        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
        variables_prev[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L_prev;

        auto const chi = [medium, x_position, t, dt](double const S_L) {
            MPL::VariableArray variables;
            variables.fill(std::numeric_limits<double>::quiet_NaN());
            variables[static_cast<int>(MPL::Variable::liquid_saturation)] = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(variables, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            t, x_position, dt, T_ip);

        auto const beta_SR =
            (1 - alpha) /
            ip_data_[ip].solid_material.getBulkModulus(t, x_position, &C_el);
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        // Set volumetric strain rate for the general case without swelling.
        variables[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(identity2.transpose() * B * u);
        variables_prev[static_cast<int>(MPL::Variable::volumetric_strain)]
            .emplace<double>(identity2.transpose() * B * (u - u_dot * dt));

        auto& phi = ip_data_[ip].porosity;
        {  // Porosity update
            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                ip_data_[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        // Swelling and possibly volumetric strain rate update.
        auto& sigma_sw = ip_data_[ip].sigma_sw;
        {
            auto const& sigma_sw_prev = ip_data_[ip].sigma_sw_prev;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;
            if (solid_phase.hasProperty(
                    MPL::PropertyType::swelling_stress_rate))
            {
                using DimMatrix = Eigen::Matrix<double, 3, 3>;
                auto const sigma_sw_dot =
                    MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                        solid_phase
                            .property(MPL::PropertyType::swelling_stress_rate)
                            .template value<DimMatrix>(
                                variables, variables_prev, x_position, t, dt));
                sigma_sw += sigma_sw_dot * dt;

                // !!! Misusing volumetric strain for mechanical volumetric
                // strain just to update the transport porosity !!!
                std::get<double>(variables[static_cast<int>(
                    MPL::Variable::volumetric_strain)]) +=
                    identity2.transpose() * C_el.inverse() * sigma_sw;
                std::get<double>(variables_prev[static_cast<int>(
                    MPL::Variable::volumetric_strain)]) +=
                    identity2.transpose() * C_el.inverse() * sigma_sw_prev;
            }

            if (solid_phase.hasProperty(MPL::PropertyType::transport_porosity))
            {
                variables_prev[static_cast<int>(
                    MPL::Variable::transport_porosity)] =
                    ip_data_[ip].transport_porosity_prev;

                ip_data_[ip].transport_porosity =
                    solid_phase.property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                    ip_data_[ip].transport_porosity;
            }
            else
            {
                variables[static_cast<int>(MPL::Variable::transport_porosity)] =
                    phi;
            }
        }
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);
        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        // For stress dependent permeability.
        variables[static_cast<int>(MPL::Variable::stress)]
            .emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    (ip_data_[ip].sigma_eff +
                     alpha * chi_S_L * identity2 * p_cap_ip)
                        .eval()));
        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        auto const& sigma_eff = ip_data_[ip].sigma_eff;
        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] =
            p_FR - (sigma_eff + sigma_sw).dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        ip_data_[ip].dry_density_solid = (1 - phi) * rho_SR;

        eps.noalias() = B * u;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));

        MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
            dthermal_strain =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_linear_thermal_expansion_coefficient) *
                T_dot_ip * dt;

        ip_data_[ip].updateConstitutiveRelationThermal(t, x_position, dt, u,
                                                       T_ip, dthermal_strain);

        auto const& b = process_data_.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = ip_data_[ip].dNdx_p;
        ip_data_[ip].v_darcy.noalias() =
            -K_over_mu * dNdx_p * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
        sigma_avg += sigma_eff;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*process_data_.element_saturation)[element_.getID()] = saturation_avg;
    (*process_data_.element_porosity)[element_.getID()] = porosity_avg;

    Eigen::Map<KV>(&(*process_data_.element_stresses)[element_.getID() *
                                                      KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(element_, is_axially_symmetric_, p_L,
                         *process_data_.pressure_interpolated);
    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(element_, is_axially_symmetric_, T,
                         *process_data_.temperature_interpolated);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
unsigned ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getNumberOfIntegrationPoints() const
{
    return integration_method_.getNumberOfPoints();
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
typename MaterialLib::Solids::MechanicsBase<
    DisplacementDim>::MaterialStateVariables const&
ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getMaterialStateVariablesAt(unsigned integration_point)
    const
{
    return *ip_data_[integration_point].material_state_variables;
}
}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
