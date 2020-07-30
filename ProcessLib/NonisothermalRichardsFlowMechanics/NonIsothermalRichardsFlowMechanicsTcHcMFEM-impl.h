/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 30, 2020, 2:56 PM
 */

#pragma once

#include <cmath>

#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/KelvinVector.h"
#include "NonIsothermalRichardsFlowMechanicsTcHcMFEM.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "WaterVaporProperty.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::
    assembleForStaggeredSchemeImplementation(double const t, double const dt,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_xdot,
                                             int const process_id,
                                             std::vector<double>& local_M_data,
                                             std::vector<double>& local_K_data,
                                             std::vector<double>& local_b_data)
{
    if (process_id == data_of_staggeredTcHcM_.T_process_id)
    {
        return assembleForEnergyBalanceEquations(t, dt, local_x, process_id,
                                                 local_M_data, local_K_data,
                                                 local_b_data);
    }

    if (process_id == data_of_staggeredTcHcM_.H_process_id)
    {
        return assembleForMassBalanceEquations(t, dt, local_x, local_xdot,
                                               process_id, local_M_data,
                                               local_K_data, local_b_data);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::
    assembleWithJacobianForStaggeredSchemeImplementation(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double /*dxdot_dx*/,
        const double /*dx_dx*/, int const /*process_id*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    assembleWithJacobianForMomentumEquations(t, dt, local_x, local_xdot,
                                             local_b_data, local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::
    assembleForEnergyBalanceEquations(double const t, double const dt,
                                      Eigen::VectorXd const& local_x,
                                      int const /* to be removed process_id*/,
                                      std::vector<double>& local_M_data,
                                      std::vector<double>& local_K_data,
                                      std::vector<double>& /*local_b_data*/)
{
    int const temperature_size = ShapeFunction::NPOINTS;
    int const pressure_size = ShapeFunction::NPOINTS;
    int const temperature_index =
        getVariableIndex(data_of_staggeredTcHcM_.T_process_id);
    int const pressure_index =
        getVariableIndex(data_of_staggeredTcHcM_.H_process_id);

    auto const nodal_p =
        local_x.template segment<pressure_size>(pressure_index);
    auto const nodal_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto local_M = MathLib::createZeroedMatrix<TemperatureMatrixTye>(
        local_M_data, temperature_size, temperature_size);
    auto local_K = MathLib::createZeroedMatrix<TemperatureMatrixTye>(
        local_K_data, temperature_size, temperature_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& medium =
        *(this->_process_data.media_map->getMedium(this->_element.getID()));
    auto const& solid_phase = medium.phase("Solid");
    auto const& gas_phase = medium.phase("Gas");
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    MPL::VariableArray vars;
    auto const& b = this->_process_data.specific_body_force;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto& ip_data = this->_ip_data[ip];
        auto const& w = ip_data.integration_weight;

        // T and H processes use the same order of interpolation
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_p, N_p, p_at_xi);
        double T_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_T, N_p, T_at_xi);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_at_xi;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            std::max(0.0, p_at_xi);
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
            std::max(0.0, -p_at_xi);

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);

        auto const rho_l =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        auto const specific_heat_capacity_fluid =
            liquid_phase.property(MaterialPropertyLib::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        auto const specific_heat_capacity_solid =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);
        auto const rho_s =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        auto& S_l = ip_data.saturation;
        S_l = medium.property(MPL::PropertyType::saturation)
                  .template value<double>(vars, pos, t, dt);

        const double vapor_density =
            data_of_staggeredTcHcM_.has_vapor_diffusion
                ? waterVaporDensity(T_at_xi, p_at_xi, rho_l)
                : 0.0;

        auto const specific_heat_capacity_gas =
            gas_phase.property(MaterialPropertyLib::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        // Assemble mass matrix. The vapor enthalpy is omitted.
        local_M.noalias() +=
            w *
            (rho_s * specific_heat_capacity_solid * (1 - porosity) +
             ((1 - S_l) * vapor_density * specific_heat_capacity_gas +
              S_l * rho_l * specific_heat_capacity_fluid) *
                 porosity) *
            N_p.transpose() * N_p;

        // Assemble Laplace matrix
        auto const& sigma_eff = ip_data.sigma_eff;

        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = S_l;
        const double bishop =
            medium.property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vars, pos, t, dt);

        // For stress dependent permeability.
        vars[static_cast<int>(MPL::Variable::stress)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                (sigma_eff - alpha_B * bishop * identity2 * p_at_xi).eval()));
        auto const intrinsic_permeability =
            MPL::formEigenTensor<DisplacementDim>(
                medium.property(MPL::PropertyType::permeability)
                    .value(vars, pos, t, dt));
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);

        double const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<double>(vars, pos, t, dt);

        GlobalDimMatrixType const hydraulic_conductivity =
            intrinsic_permeability * k_rel / viscosity;
        GlobalDimVectorType const grad_p = dNdx_p * nodal_p;
        GlobalDimVectorType const velocity_l =
            -hydraulic_conductivity * (grad_p - rho_l * b);

        auto const thermal_conductivity_solid =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .value(vars, pos, t, dt);

        auto const thermal_conductivity_fluid =
            liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(vars, pos, t, dt) *
                S_l +
            gas_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(vars, pos, t, dt) *
                (1.0 - S_l);

        auto const thermal_conductivity =
            MaterialPropertyLib::formEffectiveThermalConductivity<
                DisplacementDim>(thermal_conductivity_solid,
                                 thermal_conductivity_fluid, porosity);

        local_K.noalias() +=
            w * (dNdx_p.transpose() * thermal_conductivity * dNdx_p +
                 N_p.transpose() * velocity_l.transpose() * dNdx_p * rho_l *
                     specific_heat_capacity_fluid * porosity);

        if (data_of_staggeredTcHcM_.has_vapor_diffusion)
        {
            auto const f_Tv =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_diffusion_enhancement_factor)
                    .template value<double>(vars, pos, t, dt);
            auto const tortuosity =
                medium.property(MaterialPropertyLib::PropertyType::tortuosity)
                    .template value<double>(vars, pos, t, dt);

            GlobalDimVectorType const vapor_flux =
                f_Tv * DTv(T_at_xi, p_at_xi, rho_l, S_l, porosity, tortuosity) *
                    dNdx_p * nodal_T +
                Dpv(T_at_xi, p_at_xi, rho_l, S_l, porosity, tortuosity) *
                    grad_p;

            local_K.noalias() += w * N_p.transpose() * vapor_flux.transpose() *
                                 dNdx_p * specific_heat_capacity_gas;
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::
    assembleForMassBalanceEquations(double const t, double const dt,
                                    Eigen::VectorXd const& local_x,
                                    Eigen::VectorXd const& local_xdot,
                                    int const /*to be removed process_id*/,
                                    std::vector<double>& local_M_data,
                                    std::vector<double>& local_K_data,
                                    std::vector<double>& local_b_data)
{
    int const temperature_size = ShapeFunction::NPOINTS;
    int const pressure_size = ShapeFunction::NPOINTS;
    int const displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;

    int const temperature_index =
        getVariableIndex(data_of_staggeredTcHcM_.T_process_id);
    int const pressure_index =
        getVariableIndex(data_of_staggeredTcHcM_.H_process_id);
    int const displacement_index =
        getVariableIndex(data_of_staggeredTcHcM_.M_process_id);
    auto const nodal_p =
        local_x.template segment<pressure_size>(pressure_index);
    auto const nodal_T =
        local_x.template segment<temperature_size>(temperature_index);

    auto const T_dot =
        local_xdot.template segment<temperature_size>(temperature_index);
    auto const u_dot =
        local_xdot.template segment<displacement_size>(displacement_index);

    auto local_M = MathLib::createZeroedMatrix<TemperatureMatrixTye>(
        local_M_data, temperature_size, temperature_size);
    auto local_K = MathLib::createZeroedMatrix<TemperatureMatrixTye>(
        local_K_data, temperature_size, temperature_size);

    auto local_rhs = MathLib::createZeroedVector<PressureVectorTye>(
        local_b_data, pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& medium =
        *(this->_process_data.media_map->getMedium(this->_element.getID()));
    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    MPL::VariableArray vars;
    auto const& b = this->_process_data.specific_body_force;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(DisplacementDim, DisplacementDim));

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto& ip_data = this->_ip_data[ip];
        auto const& w = ip_data.integration_weight;

        // T and H processes use the same order of interpolation
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_p, N_p, p_at_xi);
        double T_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_T, N_p, T_at_xi);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_at_xi;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            std::max(0.0, p_at_xi);
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
            std::max(0.0, -p_at_xi);

        auto& S_l = ip_data.saturation;
        S_l = medium.property(MPL::PropertyType::saturation)
                  .template value<double>(vars, pos, t, dt);

        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);
        auto const storage =
            medium.property(MaterialPropertyLib::PropertyType::storage)
                .template value<double>(vars, pos, t, dt);

        auto const rho_l =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        auto const drho_l_dp =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::phase_pressure, pos, t,
                    dt);

        double const dSw_dpc =
            medium.property(MaterialPropertyLib::PropertyType::saturation)
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::capillary_pressure,
                    pos, t, dt);

        // Assemble mass matrix
        double mass_coefficient =
            porosity * (S_l * (storage + drho_l_dp / rho_l) - dSw_dpc);

        double rho_v_over_rho_w = 0.0;
        if (data_of_staggeredTcHcM_.has_vapor_diffusion)
        {
            rho_v_over_rho_w =
                waterVaporDensity(T_at_xi, p_at_xi, rho_l) / rho_l;

            mass_coefficient +=
                porosity * rho_v_over_rho_w *
                (dSw_dpc + (1 - S_l) / (rho_l * T_at_xi *
                                        MaterialLib::PhysicalConstant::
                                            SpecificGasConstant::WaterVapour));
        }

        local_M.noalias() += w * mass_coefficient * N_p.transpose() * N_p;

        // Assemble Laplace matrix
        auto const& sigma_eff = ip_data.sigma_eff;

        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = S_l;
        const double bishop =
            medium.property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vars, pos, t, dt);

        // For stress dependent permeability.
        vars[static_cast<int>(MPL::Variable::stress)].emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                (sigma_eff - alpha_B * bishop * identity2 * p_at_xi).eval()));
        auto const intrinsic_permeability =
            MPL::formEigenTensor<DisplacementDim>(
                medium.property(MPL::PropertyType::permeability)
                    .value(vars, pos, t, dt));
        auto const viscosity =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);

        double const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<double>(vars, pos, t, dt);

        GlobalDimMatrixType const hydraulic_conductivity =
            intrinsic_permeability * k_rel / viscosity;
        if (data_of_staggeredTcHcM_.has_vapor_diffusion)
        {
            auto const tortuosity =
                medium.property(MaterialPropertyLib::PropertyType::tortuosity)
                    .template value<double>(vars, pos, t, dt);
            local_K.noalias() +=
                w * dNdx_p.transpose() *
                (hydraulic_conductivity +
                 (Dpv(T_at_xi, p_at_xi, rho_l, S_l, porosity, tortuosity) /
                  rho_l) *
                     I) *
                dNdx_p;
        }
        else
        {
            local_K.noalias() +=
                w * dNdx_p.transpose() * hydraulic_conductivity * dNdx_p;
        }

        // Assemble RHS
        double T_dot_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(T_dot, N_p, T_dot_at_xi);

        local_rhs.noalias() +=
            w * rho_l * dNdx_p.transpose() * hydraulic_conductivity * b;

        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->_element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->_is_axially_symmetric);
        auto const eps_dot = B * u_dot;
        double const dv_dt = Invariants::trace(eps_dot);

        local_rhs.noalias() -= S_l * alpha_B * dv_dt * N_p * w;

        // Add the thermal expansion term is the point is in the fully
        // saturated zone.
        if (p_at_xi > 0.0)
        {
            auto const solid_thermal_expansion =
                solid_phase.property(MPL::PropertyType::thermal_expansivity)
                    .template value<double>(vars, pos, t, dt);
            const double drho_l_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t, dt);
            const double eff_thermal_expansion =
                3.0 * (alpha_B - porosity) * solid_thermal_expansion -
                porosity * drho_l_dT / rho_l;
            local_rhs.noalias() +=
                eff_thermal_expansion * T_dot_at_xi * w * N_p;
        }

        if (data_of_staggeredTcHcM_.has_vapor_diffusion)
        {
            double const coeff =
                (1.0 - S_l) * alpha_B * rho_v_over_rho_w * dv_dt -
                porosity * (1 - S_l) *
                    dwaterVaporDensitydT(T_at_xi, p_at_xi, rho_l) *
                    T_dot_at_xi / rho_l;
            local_rhs.noalias() += coeff * N_p * w;

            auto const f_Tv =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_diffusion_enhancement_factor)
                    .template value<double>(vars, pos, t, dt);

            auto const tortuosity =
                medium.property(MaterialPropertyLib::PropertyType::tortuosity)
                    .template value<double>(vars, pos, t, dt);
            double const D_Tv =
                DTv(T_at_xi, p_at_xi, rho_l, S_l, porosity, tortuosity) / rho_l;

            local_rhs.noalias() -=
                f_Tv * D_Tv * dNdx_p.transpose() * (dNdx_p * nodal_T) * w;
        }
    }

    if (this->_process_data.apply_mass_lumping)
    {
        for (int idx_ml = 0; idx_ml < local_M.cols(); idx_ml++)
        {
            double const mass_lump_val = local_M.col(idx_ml).sum();
            local_M.col(idx_ml).setZero();
            local_M(idx_ml, idx_ml) = mass_lump_val;
        }
    }  // end of mass lumping
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
void NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::
    assembleWithJacobianForMomentumEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data)
{
    int const temperature_size = ShapeFunction::NPOINTS;
    int const pressure_size = ShapeFunction::NPOINTS;
    int const displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    int const temperature_index =
        getVariableIndex(data_of_staggeredTcHcM_.T_process_id);
    int const pressure_index =
        getVariableIndex(data_of_staggeredTcHcM_.H_process_id);
    int const displacement_index =
        getVariableIndex(data_of_staggeredTcHcM_.M_process_id);

    auto const nodal_p =
        local_x.template segment<pressure_size>(pressure_index);
    auto const nodal_dp =
        dt * local_xdot.template segment<pressure_size>(pressure_index);

    auto const nodal_T =
        local_x.template segment<temperature_size>(temperature_index);
    auto const nodal_dT =
        dt * local_xdot.template segment<temperature_size>(temperature_index);

    auto const nodal_u =
        local_x.template segment<displacement_size>(displacement_index);

    auto local_Jac = MathLib::createZeroedMatrix<DisplacementMatrixTye>(
        local_Jac_data, displacement_size, displacement_size);
    auto local_rhs = MathLib::createZeroedVector<DisplacementVectorTye>(
        local_b_data, displacement_size);

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material =
        *this->_process_data.solid_materials[0];

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const& medium =
        *(this->_process_data.media_map->getMedium(this->_element.getID()));
    auto const& solid_phase = medium.phase("Solid");
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    MPL::VariableArray vars;
    auto const& b = this->_process_data.specific_body_force;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;

    unsigned const n_integration_points =
        this->_integration_method.getNumberOfPoints();

    for (unsigned ip(0); ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip);

        auto& ip_data = this->_ip_data[ip];
        auto const& w = ip_data.integration_weight;

        // T and H processes use the same order of interpolation
        auto const& N_p = ip_data.N_p;
        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;

        double p_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_p, N_p, p_at_xi);
        double T_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_T, N_p, T_at_xi);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_at_xi;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            std::max(0.0, p_at_xi);
        vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
            std::max(0.0, -p_at_xi);

        auto& S_l = ip_data.saturation;
        S_l = medium.property(MPL::PropertyType::saturation)
                  .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(MPL::Variable::liquid_saturation)] = S_l;

        auto const rho_l =
            liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        auto const rho_s =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        auto const alpha_B =
            medium.property(MPL::PropertyType::biot_coefficient)
                .template value<double>(vars, pos, t, dt);

        const double bishop =
            medium.property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vars, pos, t, dt);
        auto const porosity =
            medium.property(MaterialPropertyLib::PropertyType::porosity)
                .template value<double>(vars, pos, t, dt);

        // Assemble Jacobian

        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->_element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->_is_axially_symmetric);
        auto& eps = ip_data.eps;
        auto const& eps_prev = ip_data.eps_prev;
        eps.noalias() = B * nodal_u;

        // assume isotropic thermal expansion
        double dT_at_xi = 0.;
        NumLib::shapeFunctionInterpolate(nodal_dT, N_p, dT_at_xi);
        auto const solid_thermal_expansion =
            solid_phase.property(MPL::PropertyType::thermal_expansivity)
                .template value<double>(vars, pos, t, dt) *
            dT_at_xi;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& sigma_eff_prev = ip_data.sigma_eff_prev;
        auto& state = ip_data.material_state_variables;
        auto sigma_eff0 = sigma_eff_prev;

        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            double dp_at_xi = 0.;
            NumLib::shapeFunctionInterpolate(nodal_dp, N_p, dp_at_xi);
            vars[static_cast<int>(MPL::Variable::capillary_pressure)] =
                std::max(0.0, -p_at_xi + dp_at_xi);

            double const S_l_0 = medium.property(MPL::PropertyType::saturation)
                                     .template value<double>(vars, pos, t, dt);
            NumLib::shapeFunctionInterpolate(nodal_dp, N_p, dp_at_xi);
            vars[static_cast<int>(MPL::Variable::liquid_saturation_rate)] =
                (S_l - S_l_0) / dt;

            double const dswelling_stress =
                solid_phase.property(MPL::PropertyType::swelling_stress_rate)
                    .template value<double>(vars, pos, t, dt);
            sigma_eff0 -= dswelling_stress * Invariants::identity2;
        }

        auto&& solution = solid_material.integrateStress(
            t, pos, dt, eps_prev,
            eps - solid_thermal_expansion * Invariants::identity2, sigma_eff0,
            *state, T_at_xi);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);

        local_Jac.noalias() += B.transpose() * C * B * w;

        // Assemble RHS

        // porosity * vapor_density * (1 - S_l) is omitted due to its relative
        // tiny value.
        double const rho = rho_s * (1 - porosity) + porosity * rho_l * S_l;

        auto const& N_u_op = ip_data.N_u_op;
        local_rhs.noalias() -=
            (B.transpose() *
                 (sigma_eff - alpha_B * bishop * identity2 * p_at_xi) -
             N_u_op.transpose() * rho * b) *
            w;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          typename IntegrationMethod, int DisplacementDim>
int NonIsothermalRichardsFlowMechanicsAssemblerTcHcM<
    ShapeFunctionDisplacement, ShapeFunction, IntegrationMethod,
    DisplacementDim>::getVariableIndex(const int process_id) const
{
    // Index of displacement
    if (process_id == data_of_staggeredTcHcM_.M_process_id)
    {
        return data_of_staggeredTcHcM_.M_process_id * ShapeFunction::NPOINTS;
    }

    // Index of pressure or temperature
    switch (process_id)
    {
        case 0:
            return 0;
        case 1:
            return data_of_staggeredTcHcM_.M_process_id == 0
                       ? ShapeFunctionDisplacement::NPOINTS * DisplacementDim
                       : ShapeFunction::NPOINTS;
        case 2:
            return ShapeFunctionDisplacement::NPOINTS * DisplacementDim +
                   ShapeFunction::NPOINTS;
    }
    return 0;  // to avoid warning.
}

}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
