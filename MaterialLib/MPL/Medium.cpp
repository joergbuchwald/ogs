/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Medium.h"

#include "BaseLib/Algorithm.h"
#include "Properties/Properties.h"

namespace MaterialPropertyLib
{
Medium::Medium(std::vector<std::unique_ptr<Phase>>&& phases,
               std::unique_ptr<PropertyArray>&& properties)
    : phases_(std::move(phases))
{
    if (properties)
    {
        overwriteExistingProperties(properties_, *properties, this);
    }
}

Phase const& Medium::phase(std::size_t const index) const
{
    return *phases_[index];
}

Phase const& Medium::phase(std::string const& name) const
{
    return *BaseLib::findElementOrError(
        phases_.begin(), phases_.end(),
        [&name](std::unique_ptr<MaterialPropertyLib::Phase> const& phase) {
            return phase->name == name;
        },
        "Could not find phase name '" + name + "'.");
}

Property const& Medium::property(PropertyType const& p) const
{
    Property const* const property = properties_[p].get();
    if (property == nullptr)
    {
        OGS_FATAL("Trying to access undefined property '{:s}' of {:s}",
                  property_enum_to_string[p], description());
    }
    return *properties_[p];
}

bool Medium::hasProperty(PropertyType const& p) const
{
    return properties_[p] != nullptr;
}

std::size_t Medium::numberOfPhases() const
{
    return phases_.size();
}

std::string Medium::description()
{
    return "medium";
}
}  // namespace MaterialPropertyLib
