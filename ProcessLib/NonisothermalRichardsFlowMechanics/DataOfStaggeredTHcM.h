/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 11, 2020, 9:40 AM
 */

#pragma once
namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
struct DataOfStaggeredTHcM
{
    const int TH_process_id;
    const int M_process_id;
    const bool has_vapor_diffusion;
};
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
