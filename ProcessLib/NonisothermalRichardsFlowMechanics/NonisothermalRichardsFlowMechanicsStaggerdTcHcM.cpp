/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 29, 2020, 9:32 AM
 */

#include "NonisothermalRichardsFlowMechanicsStaggerdTcHcM.h"

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Deformation/SolidMaterialInternalToSecondaryVariables.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/RichardsMechanics/CreateLocalAssemblers.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace NonisothermalRichardsFlowMechanics
{
template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<2>;
template class NonisothermalRichardsFlowMechanicsStaggerdTcHcM<3>;
}  // namespace NonisothermalRichardsFlowMechanics
}  // namespace ProcessLib
