/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SaturationDependentHeatConduction.h"

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
SaturationDependentHeatConduction::SaturationDependentHeatConduction(
    std::string name, double const K_dry, double const K_wet)
    : K_dry_(K_dry), K_wet_(K_wet)
{
    name_ = std::move(name);
}

void SaturationDependentHeatConduction::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'SaturationDependentHeatConduction' is "
            "implemented on the 'medium' scale only.");
    }
}

PropertyDataType SaturationDependentHeatConduction::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    return K_dry_ * (1 - S_L) + K_wet_ * S_L;
}

PropertyDataType SaturationDependentHeatConduction::dValue(
    VariableArray const& /*variable_array*/, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)variable;
    assert((variable == Variable::liquid_saturation) &&
           "SaturationDependentHeatConduction::dvalue is implemented for "
           "derivatives with "
           "respect to liquid saturation only.");

    return K_wet_ - K_dry_;
}
}  // namespace MaterialPropertyLib
