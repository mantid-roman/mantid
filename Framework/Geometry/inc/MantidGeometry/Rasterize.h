// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_GEOMETRY_RASTERIZE_H_
#define MANTID_GEOMETRY_RASTERIZE_H_

#include "MantidGeometry/DllConfig.h"
#include "MantidKernel/V3D.h"
#include <boost/shared_ptr.hpp>

namespace Mantid {
namespace Geometry {

class IObject;
class CSGObject;

/// Holds the information used for doing numerical integrations of object in the
/// beam
struct MANTID_GEOMETRY_DLL Raster {
public:
  void reserve(size_t numVolumeElements);

  std::vector<double> l1;            ///< Cached L1 distances
  std::vector<double> volume;        ///< Cached element volumes
  std::vector<Kernel::V3D> position; ///< Cached element positions
  double totalvolume;                ///< Volume of the object
};

namespace Rasterize {

MANTID_GEOMETRY_DLL Raster calculate(const Kernel::V3D &beamDirection,
                                     const Geometry::CSGObject &shape,
                                     const double cubeSizeInMetre);

MANTID_GEOMETRY_DLL Raster calculate(const Kernel::V3D &beamDirection,
                                     const boost::shared_ptr<const IObject> shape,
                                     const double cubeSizeInMetre);

MANTID_GEOMETRY_DLL Raster calculateCylinder(const Kernel::V3D &beamDirection,
                                             const Geometry::CSGObject &shape,
                                             const size_t numSlices,
                                             const size_t numAnnuli);

MANTID_GEOMETRY_DLL Raster calculateCylinder(
    const Kernel::V3D &beamDirection, const boost::shared_ptr<const IObject> shape,
    const size_t numSlices, const size_t numAnnuli);

} // namespace Rasterize
} // namespace Geometry
} // namespace Mantid

#endif /* MANTID_GEOMETRY_RASTERIZE_H_ */
