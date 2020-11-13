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

#pragma once

#include "MaterialLib/MPL/VariableType.h"  // for VariableArray
#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * It gets the liquid compressibility coefficient.
 *
 * If the the liquid compressibility coefficient is given in the project file via
 * the media property of compressibility, e.g
 * \verbatim
 *     <property>
 *       <name>compressibility</name>
 *       <type>Constant</type>
 *       <value>1e-10</value>
 *     </property>
 * \endverbatim
 * or
 * \verbatim
 *     <property>
 *       <name>bulk_modulus</name>
 *       <type>Constant</type>
 *       <value>1e10</value>
 *     </property>
 * \endverbatim
 * it returns the value of the given property. Otherwise it returns the value
 * computed from the density model by the following formula
 * \f[
 *      (\frac{\partial \rho}{\partial p})/\rho
 * \f]
 * where \f$\rho\f$ is the density, \f$p\f$ is the phase pressure.
 */
double getLiquidCompressibility(Phase const& phase,
                                   VariableArray const& vars,
                                   const double density,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t, double const dt);
}  // namespace MaterialPropertyLib
