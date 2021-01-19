/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "EquilibriumReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void EquilibriumReactant::print(std::ostream& os,
                                std::size_t const chemical_system_id) const
{
    os << name << " " << saturation_index << " "
       << (*molality)[chemical_system_id] << "\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
