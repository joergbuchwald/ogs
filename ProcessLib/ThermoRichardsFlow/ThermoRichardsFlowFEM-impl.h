/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cassert>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "WaterVaporProperty.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    ThermoRichardsFlowLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoRichardsFlowProcessData& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& medium = _process_data.media_map->getMedium(_element.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = shape_matrices[ip];
        x_position.setIntegrationPoint(ip);
        _ip_data.emplace_back();
        auto& ip_data = _ip_data[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm.integralMeasure * sm.detJ;

        ip_data.N_p = sm.N;
        ip_data.dNdx_p = sm.dNdx;

        // Initial porosity. Could be read from intergration point data or mesh.
        ip_data.porosity =
            medium->property(MPL::porosity)
                .template initialValue<double>(
                    x_position,
                    std::numeric_limits<
                        double>::quiet_NaN() /* t independent */);
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::size_t ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element {:d} is different "
            "from the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "saturation_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::saturation);
    }
    if (name == "porosity_ip")
    {
        return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                         &IpData::porosity);
    }
    return 0;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::setInitialConditionsConcrete(std::vector<double> const& local_x,
                                             double const t, bool const /*use_monolithic_scheme*/,
                                             int const /*process_id*/)
{
    assert(local_x.size() == temperature_size + pressure_size);

    auto p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    MPL::VariableArray variables;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& N_p = _ip_data[ip].N_p;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap_ip;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = -p_cap_ip;

        // Note: temperature dependent saturation model is not considered so
        // far.
        _ip_data[ip].saturation_prev =
            medium->property(MPL::PropertyType::saturation)
                .template value<double>(
                    variables, x_position, t,
                    std::numeric_limits<double>::quiet_NaN());
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assembleWithJacobian(double const t, double const dt,
                                     std::vector<double> const& local_x,
                                     std::vector<double> const& local_xdot,
                                     const double /*dxdot_dx*/,
                                     const double /*dx_dx*/,
                                     std::vector<double>& /*local_M_data*/,
                                     std::vector<double>& /*local_K_data*/,
                                     std::vector<double>& local_rhs_data,
                                     std::vector<double>& local_Jac_data)
{
    auto const local_matrix_dim = pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
    auto const p_L_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_xdot.data() + pressure_index, pressure_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<local_matrix_dim,
                                                        local_matrix_dim>>(
        local_Jac_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<local_matrix_dim>>(
        local_rhs_data, local_matrix_dim);

    typename ShapeMatricesType::NodalMatrixType M_TT =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType K_TT =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType K_Tp =
        ShapeMatricesType::NodalMatrixType::Zero(temperature_size,
                                                 pressure_size);
    typename ShapeMatricesType::NodalMatrixType M_pT =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size,
                                                 temperature_size);
    typename ShapeMatricesType::NodalMatrixType laplace_p =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_p =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    typename ShapeMatricesType::NodalMatrixType storage_p_a_S =
        ShapeMatricesType::NodalMatrixType::Zero(pressure_size, pressure_size);

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

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

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);

        auto storage_correction = 0.0;
        if (medium->hasProperty(MPL::PropertyType::storage_correction))
        {
            storage_correction =
                medium->property(MPL::PropertyType::storage_correction)
                    .template value<double>(variables, x_position, t, dt);
        }
        auto thermal_expansivity_correction = 0.0;
        if (medium->hasProperty(MPL::PropertyType::storage_correction))
        {
            thermal_expansivity_correction =
                medium->property(MPL::PropertyType::thermal_expansivity_correction)
                    .template value<double>(variables, x_position, t, dt);
        }

        //bulk_modulus correct name for bulk modulus of solid skeleton
        auto const K_S =
            solid_phase.property(MPL::PropertyType::bulk_modulus)
                .template value<double>(variables, x_position, t, dt);
        auto const beta_SR = (1 - alpha) / K_S;
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        auto const K_LR =
            liquid_phase.property(MPL::PropertyType::bulk_modulus)
                .template value<double>(variables, x_position, t, dt);

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        auto const& b = _process_data.specific_body_force;

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

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update

            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                _ip_data[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);
 
        auto const K_intrinsic = MPL::formEigenTensor<GlobalDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const Ki_over_mu = K_intrinsic / mu;
        GlobalDimMatrixType const rho_Ki_over_mu = rho_LR * Ki_over_mu;

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

        //double const dthermal_strain =
        //    solid_linear_thermal_expansion_coefficient.trace() * T_dot_ip * dt;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
         variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] = p_FR;
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        //double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;

        //
        // pressure equation, pressure part.
        //
        laplace_p.noalias() +=
            dNdx_p.transpose() * k_rel * rho_Ki_over_mu * dNdx_p * w;

        double const a0 = (alpha > phi) ? 0.0 : (alpha - phi) * beta_SR;
        double const specific_storage_a_p = S_L * (phi / K_LR + S_L * a0 +
                storage_correction);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi / K_LR + 2 * S_L * a0 + storage_correction);
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
        // Note:  d (rho_l * beta _T)/dp * dotT
        // Add the thermal expansion term is the point is in the fully
        // saturated zone.
        if (p_cap_ip <= 0.0)  // p_l>0.0
        {
            double const fluid_volumetric_thermal_expansion_coefficient =
                MPL::getLiquidThermalExpansivity(liquid_phase, variables,
                                                 rho_LR, x_position, t, dt);
            const double eff_thermal_expansion =
                (alpha - phi) *
                    solid_linear_thermal_expansion_coefficient.trace() +
                phi * fluid_volumetric_thermal_expansion_coefficient +
                thermal_expansivity_correction;
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
                    GlobalDim>(thermal_conductivity_solid,
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

        if (_process_data.has_water_vaporization)
        {
            double const p_L_ip = -p_cap_ip;
            double const rho_wv = waterVaporDensity(T_ip, p_L_ip, rho_LR);
            double const storage_coefficient_by_water_vapor =
                phi * rho_wv *
                (dS_L_dp_cap +
                 (1 - S_L) / (rho_LR * T_ip *
                              MaterialLib::PhysicalConstant::
                                  SpecificGasConstant::WaterVapour));
            storage_p_a_p.noalias() +=
                N_p.transpose() * storage_coefficient_by_water_vapor * N_p * w;

            double const vapor_expansion =
                phi * (1 - S_L) * dwaterVaporDensitydT(T_ip, p_L_ip, rho_LR);
            M_pT.noalias() += N_p.transpose() * vapor_expansion * N_p * w;

            auto const f_Tv =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_diffusion_enhancement_factor)
                    .template value<double>(variables, x_position, t, dt);

            auto const tortuosity =
                medium->property(MaterialPropertyLib::PropertyType::tortuosity)
                    .template value<double>(variables, x_position, t, dt);
            double const f_Tv_D_Tv =
                f_Tv * DTv(T_ip, p_L_ip, rho_LR, S_L, phi, tortuosity);

            local_Jac
                .template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index)
                .noalias() += dNdx_p.transpose() * f_Tv_D_Tv * dNdx_p * w;

            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() -= f_Tv_D_Tv * dNdx_p.transpose() * (dNdx_p * T) * w;

            laplace_p.noalias() +=
                dNdx_p.transpose() *
                Dpv(T_ip, p_L_ip, rho_LR, S_L, phi, tortuosity) * dNdx_p * w;
        }
    }

    if (_process_data.apply_mass_lumping)
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

    //
    // -- Residual
    //
    // temperature equation
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        M_TT * T_dot + K_TT * T;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L + (storage_p_a_p + storage_p_a_S) * p_L_dot +
        M_pT * T_dot;
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

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        cache_matrix.col(ip).noalias() = _ip_data[ip].v_darcy;
    }

    return cache;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getSaturation() const
{
    std::vector<double> result;
    getIntPtSaturation(0, {}, {}, result);
    return result;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::saturation, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getPorosity() const
{
    std::vector<double> result;
    getIntPtPorosity(0, {}, {}, result);
    return result;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                     &IpData::porosity, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
std::vector<double> const&
ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDryDensitySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
{
    return ProcessLib::getIntegrationPointScalarData(
        _ip_data, &IpData::dry_density_solid, cache);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod,
                                      GlobalDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                std::vector<double> const& local_xdot,
                                double const /*t*/, double const /*dt*/,
                                bool const /*use_monolithic_scheme*/,
                                int const /*process_id*/)
{
    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

    //auto const& medium = _process_data.media_map->getMedium(_element.getID());
    MPL::VariableArray variables;
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    //auto const& solid_phase = medium->phase("Solid");

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = _ip_data[ip].N_p;
        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N_p, T_ip);
        double T_dot_ip;
        NumLib::shapeFunctionInterpolate(T_dot, N_p, T_dot_ip);
        variables[static_cast<int>(MPL::Variable::temperature)] = T_ip;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        /*Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));
        */
        //double const dthermal_strain =
        //    solid_linear_thermal_expansion_coefficient.trace() * T_dot_ip * dt;
    }
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermoRichardsFlowLocalAssembler<ShapeFunction, IntegrationMethod,
                                      GlobalDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     std::vector<double> const& local_x,
                                     std::vector<double> const& local_x_dot)
{
    auto const T = Eigen::Map<typename ShapeMatricesType::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);
    auto const T_dot =
        Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x_dot.data() + temperature_index,
                                     temperature_size);

    auto const p_L = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x.data() + pressure_index, pressure_size);

    auto p_L_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<pressure_size> const>(
        local_x_dot.data() + pressure_index, pressure_size);

    auto const& medium = _process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_p = _ip_data[ip].N_p;

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

        auto& S_L = _ip_data[ip].saturation;
        auto const S_L_prev = _ip_data[ip].saturation_prev;
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

        auto const K_SR =
            solid_phase.property(MPL::PropertyType::bulk_modulus)
                .template value<double>(variables, x_position, t, dt);
        auto const beta_SR = (1 - alpha) / K_SR;
        variables[static_cast<int>(MPL::Variable::grain_compressibility)] =
            beta_SR;

        variables[static_cast<int>(MPL::Variable::effective_pore_pressure)] =
            -chi_S_L * p_cap_ip;
        variables_prev[static_cast<int>(
            MPL::Variable::effective_pore_pressure)] =
            -chi_S_L_prev * (p_cap_ip - p_cap_dot_ip * dt);

        auto& phi = _ip_data[ip].porosity;
        {  // Porosity update
            variables_prev[static_cast<int>(MPL::Variable::porosity)] =
                _ip_data[ip].porosity_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables[static_cast<int>(MPL::Variable::porosity)] = phi;
        }

        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);
        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        auto const K_intrinsic = MPL::formEigenTensor<GlobalDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables[static_cast<int>(MPL::Variable::solid_grain_pressure)] = p_FR;
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        _ip_data[ip].dry_density_solid = (1 - phi) * rho_SR;

        // Consider anisotropic thermal expansion.
        // Read in 3x3 tensor. 2D case also requires expansion coeff. for z-
        // component.
        /*Eigen::Matrix<double, 3,
                      3> const solid_linear_thermal_expansion_coefficient =
            MaterialPropertyLib::formEigenTensor<3>(
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_expansivity)
                    .value(variables, x_position, t, dt));
        */
        //double const dthermal_strain =
        //    solid_linear_thermal_expansion_coefficient.trace() * T_dot_ip * dt;

        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        _ip_data[ip].v_darcy.noalias() =
            -K_over_mu * dNdx_p * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;

    (*_process_data.element_saturation)[_element.getID()] = saturation_avg;
    (*_process_data.element_porosity)[_element.getID()] = porosity_avg;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
unsigned ThermoRichardsFlowLocalAssembler<
    ShapeFunction, IntegrationMethod, GlobalDim>::getNumberOfIntegrationPoints()
    const
{
    return _integration_method.getNumberOfPoints();
}

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
