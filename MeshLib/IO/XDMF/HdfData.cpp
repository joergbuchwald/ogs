/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "HdfData.h"

#include <hdf5.h>

#include <map>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "partition.h"

namespace MeshLib::IO
{
static hid_t meshPropertyType2HdfType(MeshPropertyDataType const ogs_data_type)
{
    std::map<MeshPropertyDataType const, hid_t> ogs_to_hdf_type = {
        {MeshPropertyDataType::float64, H5T_IEEE_F64LE},
        {MeshPropertyDataType::float32, H5T_IEEE_F32LE},
        {MeshPropertyDataType::int32, H5T_STD_I32LE},
        {MeshPropertyDataType::int64, H5T_STD_I64LE},
        {MeshPropertyDataType::uint32, H5T_STD_U32LE},
        {MeshPropertyDataType::uint64, H5T_STD_U64LE},
        {MeshPropertyDataType::int8, H5T_STD_I8LE},
        {MeshPropertyDataType::uint8, H5T_STD_U8LE}};
    try
    {
        return ogs_to_hdf_type.at(ogs_data_type);
    }
    catch (std::exception const& e)
    {
        OGS_FATAL("No known HDF5 type for OGS type. {:s}", e.what());
    }
}

HdfData::HdfData(void const* data_start, std::size_t const size_partitioned_dim,
                 std::size_t const size_tuple, std::string const& name,
                 MeshPropertyDataType const mesh_property_data_type)
    : data_start(data_start),
      data_space{size_partitioned_dim, size_tuple},
      name(name)
{
    auto const& partition_info = getPartitionInfo(size_partitioned_dim);
    DBUG("HdfData: The partition of dataset {:s} has dimension {:d} and offset {:d}.",
         name, size_partitioned_dim, partition_info.first);
    auto const& offset_partitioned_dim = partition_info.first;
    offsets = {offset_partitioned_dim, 0};
    file_space = {partition_info.second, size_tuple};
    data_type = meshPropertyType2HdfType(mesh_property_data_type);
}
}  // namespace MeshLib::IO