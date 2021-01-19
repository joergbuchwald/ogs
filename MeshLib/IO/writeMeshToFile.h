/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <filesystem>
#include <string>

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
int writeMeshToFile(MeshLib::Mesh const& mesh,
                    std::filesystem::path const& file_path);
}  // namespace MeshLib::IO