/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 13, 2020, 1:00 PM
 */

#include "GetLiquidCompressibility.h"

#include "MaterialLib/MPL/Phase.h"

namespace MaterialPropertyLib
{
class Phase;

double getLiquidCompressibility(Phase const& phase,
                                   VariableArray const& vars,
                                   const double density,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt)
{
    // The thermal expansivity is explicitly given in the project file.
    if (phase.hasProperty(
            MaterialPropertyLib::PropertyType::compressibility))
    {
        return phase
            .property(MaterialPropertyLib::PropertyType::compressibility)
            .template value<double>(vars, pos, t, dt);
    }
    else if (phase.hasProperty(
            MaterialPropertyLib::PropertyType::bulk_modulus))
    {
        return 1 / phase
            .property(MaterialPropertyLib::PropertyType::bulk_modulus)
            .template value<double>(vars, pos, t, dt);
    }
    else
    {

    // The thermal expansivity calculated by the density model directly.
    return (density == 0.0)
               ? 0.0
               : phase.property(MaterialPropertyLib::PropertyType::density)
                         .template dValue<double>(
                             vars, MaterialPropertyLib::Variable::phase_pressure,
                             pos, t, dt) /
                     density;
    }
}
}  // namespace MaterialPropertyLib
