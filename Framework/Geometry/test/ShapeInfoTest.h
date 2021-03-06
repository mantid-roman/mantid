// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidGeometry/Rendering/ShapeInfo.h"
#include "MantidKernel/V3D.h"
#include <cxxtest/TestSuite.h>

using Mantid::Kernel::V3D;
using namespace Mantid::Geometry::detail;

class ShapeInfoTest : public CxxTest::TestSuite {
public:
  void testConstructEmptyInitiazesEverythingZero() {
    ShapeInfo shapeInfo;
    TS_ASSERT(shapeInfo.points().size() == 0);
    TS_ASSERT(shapeInfo.height() == 0);
    TS_ASSERT(shapeInfo.radius() == 0);
    TS_ASSERT(shapeInfo.shape() == ShapeInfo::GeometryShape::NOSHAPE);
  }

  void testSetSphere() {
    ShapeInfo shapeInfo;
    V3D center(0, 0, 0);
    double radius = 10;
    shapeInfo.setSphere(center, radius);

    TS_ASSERT_EQUALS(shapeInfo.shape(), ShapeInfo::GeometryShape::SPHERE);
    TS_ASSERT_EQUALS(shapeInfo.radius(), radius);
    TS_ASSERT_EQUALS(shapeInfo.height(), 0);
    TS_ASSERT_EQUALS(shapeInfo.points().size(), 1);
    TS_ASSERT_EQUALS(shapeInfo.points()[0], center);
  }

  void testSetCuboid() {
    ShapeInfo shapeInfo;
    V3D p1(0, 0, 0);
    V3D p2(0, 0, 1);
    V3D p3(0, 1, 0);
    V3D p4(0, 1, 1);

    shapeInfo.setCuboid(p1, p2, p3, p4);

    TS_ASSERT_EQUALS(shapeInfo.shape(), ShapeInfo::GeometryShape::CUBOID);
    TS_ASSERT_EQUALS(shapeInfo.radius(), 0);
    TS_ASSERT_EQUALS(shapeInfo.height(), 0);
    TS_ASSERT_EQUALS(shapeInfo.points().size(), 4);
    TS_ASSERT_EQUALS(shapeInfo.points(), (std::vector<V3D>{p1, p2, p3, p4}));
  }

  void testSetHexahedron() {
    ShapeInfo shapeInfo;
    V3D p1(0, 0, 0);
    V3D p2(0, 0, 1);
    V3D p3(0, 1, 0);
    V3D p4(0, 1, 1);
    V3D p5(1, 0, 0);
    V3D p6(1, 0, 1);
    V3D p7(1, 1, 0);
    V3D p8(1, 1, 1);

    shapeInfo.setHexahedron(p1, p2, p3, p4, p5, p6, p7, p8);

    TS_ASSERT_EQUALS(shapeInfo.shape(), ShapeInfo::GeometryShape::HEXAHEDRON);
    TS_ASSERT_EQUALS(shapeInfo.radius(), 0);
    TS_ASSERT_EQUALS(shapeInfo.height(), 0);
    TS_ASSERT_EQUALS(shapeInfo.points().size(), 8);
    TS_ASSERT_EQUALS(shapeInfo.points(),
                     (std::vector<V3D>{p1, p2, p3, p4, p5, p6, p7, p8}));
  }

  void testSetCone() {
    ShapeInfo shapeInfo;
    V3D center(0, 0, 0);
    V3D axis(1, 0, 0);
    double radius = 10;
    double height = 5;
    shapeInfo.setCone(center, axis, radius, height);

    TS_ASSERT_EQUALS(shapeInfo.shape(), ShapeInfo::GeometryShape::CONE);
    TS_ASSERT_EQUALS(shapeInfo.radius(), radius);
    TS_ASSERT_EQUALS(shapeInfo.height(), height);
    TS_ASSERT_EQUALS(shapeInfo.points().size(), 2);
    TS_ASSERT_EQUALS(shapeInfo.points()[0], center);
    TS_ASSERT_EQUALS(shapeInfo.points()[1], axis);
  }

  void testSetCylinder() {
    ShapeInfo shapeInfo;
    V3D center(0, 0, 0);
    V3D axis(1, 0, 0);
    double radius = 10;
    double height = 5;
    shapeInfo.setCylinder(center, axis, radius, height);

    TS_ASSERT_EQUALS(shapeInfo.shape(), ShapeInfo::GeometryShape::CYLINDER);
    TS_ASSERT_EQUALS(shapeInfo.radius(), radius);
    TS_ASSERT_EQUALS(shapeInfo.height(), height);
    TS_ASSERT_EQUALS(shapeInfo.points().size(), 2);
    TS_ASSERT_EQUALS(shapeInfo.points()[0], center);
    TS_ASSERT_EQUALS(shapeInfo.points()[1], axis);
  }

  void testGetObjectGeometry() {
    ShapeInfo shapeInfo;
    V3D center(0, 0, 0);
    double radius = 10;
    shapeInfo.setSphere(center, radius);

    ShapeInfo::GeometryShape tshape;
    std::vector<V3D> tpoints;
    double theight;
    double tradius;

    shapeInfo.getObjectGeometry(tshape, tpoints, tradius, theight);
    TS_ASSERT_EQUALS(tradius, radius);
    TS_ASSERT(theight == 0);
    TS_ASSERT(tpoints.size() == 1);
    TS_ASSERT_EQUALS(tpoints[0], center);
    TS_ASSERT_EQUALS(tshape, ShapeInfo::GeometryShape::SPHERE);
  }

  void testCuboidGeometry() {
    ShapeInfo shapeInfo;
    const V3D p1(0, 0, 0);
    const V3D p2(0, 0, 1);
    const V3D p3(0, 1, 0);
    const V3D p4(0, 1, 1);

    shapeInfo.setCuboid(p1, p2, p3, p4);
    const auto geometry = shapeInfo.cuboidGeometry();
    TS_ASSERT_EQUALS(geometry.leftFrontBottom, p1)
    TS_ASSERT_EQUALS(geometry.leftFrontTop, p2)
    TS_ASSERT_EQUALS(geometry.leftBackBottom, p3)
    TS_ASSERT_EQUALS(geometry.rightFrontBottom, p4)
  }

  void testHexahedronGeometry() {
    ShapeInfo shapeInfo;
    V3D p1(0, 0, 0);
    V3D p2(0, 0, 1);
    V3D p3(0, 1, 0);
    V3D p4(0, 1, 1);
    V3D p5(1, 0, 0);
    V3D p6(1, 0, 1);
    V3D p7(1, 1, 0);
    V3D p8(1, 1, 1);

    shapeInfo.setHexahedron(p1, p2, p3, p4, p5, p6, p7, p8);
    const auto geometry = shapeInfo.hexahedronGeometry();
    TS_ASSERT_EQUALS(geometry.leftBackBottom, p1)
    TS_ASSERT_EQUALS(geometry.leftFrontBottom, p2)
    TS_ASSERT_EQUALS(geometry.rightFrontBottom, p3)
    TS_ASSERT_EQUALS(geometry.rightBackBottom, p4)
    TS_ASSERT_EQUALS(geometry.leftBackTop, p5)
    TS_ASSERT_EQUALS(geometry.leftFrontTop, p6)
    TS_ASSERT_EQUALS(geometry.rightFrontTop, p7)
    TS_ASSERT_EQUALS(geometry.rightBackTop, p8)
  }

  void testSphereGeometry() {
    ShapeInfo shapeInfo;
    const V3D center(0, 0, 0);
    constexpr double radius = 10;
    shapeInfo.setSphere(center, radius);
    const auto geometry = shapeInfo.sphereGeometry();
    TS_ASSERT_EQUALS(geometry.centre, center)
    TS_ASSERT_EQUALS(geometry.radius, radius)
  }

  void testCylinderGeometry() {
    ShapeInfo shapeInfo;
    const V3D center(0, 0, 0);
    const V3D axis(1, 0, 0);
    constexpr double radius = 10;
    constexpr double height = 5;
    shapeInfo.setCylinder(center, axis, radius, height);
    const auto geometry = shapeInfo.cylinderGeometry();
    TS_ASSERT_EQUALS(geometry.centreOfBottomBase, center)
    TS_ASSERT_EQUALS(geometry.axis, axis)
    TS_ASSERT_EQUALS(geometry.radius, radius)
    TS_ASSERT_EQUALS(geometry.height, height)
  }

  void testConeGeometry() {
    ShapeInfo shapeInfo;
    const V3D center(0, 0, 0);
    const V3D axis(1, 0, 0);
    constexpr double radius = 10;
    constexpr double height = 5;
    shapeInfo.setCone(center, axis, radius, height);
    const auto geometry = shapeInfo.coneGeometry();
    TS_ASSERT_EQUALS(geometry.centre, center)
    TS_ASSERT_EQUALS(geometry.axis, axis)
    TS_ASSERT_EQUALS(geometry.radius, radius)
    TS_ASSERT_EQUALS(geometry.height, height)
  }

  void testCopyConstructor() {
    ShapeInfo shapeInfo;
    V3D center(0, 2, 1);
    double radius = 10;
    shapeInfo.setSphere(center, radius);

    ShapeInfo shapeInfoCopy(shapeInfo);

    TS_ASSERT_EQUALS(shapeInfo.shape(), shapeInfoCopy.shape());
    TS_ASSERT_EQUALS(shapeInfo.radius(), shapeInfoCopy.radius());
    TS_ASSERT_EQUALS(shapeInfo.height(), shapeInfoCopy.height());
    TS_ASSERT_EQUALS(shapeInfo.points(), shapeInfoCopy.points());
  }

  void testEquality() {
    ShapeInfo shapeInfo;
    V3D center(0, 2, 1);
    double radius = 10;
    shapeInfo.setSphere(center, radius);
    ShapeInfo shapeInfo2, shapeInfo3;
    shapeInfo2.setSphere(center, radius);
    shapeInfo3.setCuboid(V3D(), V3D(), V3D(), V3D());
    TS_ASSERT_EQUALS(shapeInfo2, shapeInfo);
    TS_ASSERT_DIFFERS(shapeInfo3, shapeInfo);
  }
};
