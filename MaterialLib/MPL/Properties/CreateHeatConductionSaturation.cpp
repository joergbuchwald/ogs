/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "HeatConductionSaturation.h"

namespace MaterialPropertyLib
{
std::unique_ptr<HeatConductionSaturationDependent>
createHeatConductionSaturation(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "HeatConductionSaturationDependent");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create saturation dependent thermal_conductivity property {:s}.",
         property_name);

    auto const K_dry =
        //! \ogs_file_param{properties__property__HeatConductionSaturationDependent__value_dry}
        config.getConfigParameter<double>("value_dry");

    auto const K_wet =
        //! \ogs_file_param{properties__property__HeatConductionSaturationDependent__value_wet}
        config.getConfigParameter<double>("value_wet");

    return std::make_unique<
        MaterialPropertyLib::HeatConductionSaturationDependent>(
        std::move(property_name), K_dry, K_wet);
}
}  // namespace MaterialPropertyLib
