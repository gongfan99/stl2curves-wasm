#pragma once
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Standard_Version.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#if OCC_VERSION_HEX < 0x080000
#include <TColStd_HArray1OfReal.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#endif

namespace stl2curves::detail {
using RealArray = NCollection_Array1<double>;
using IntArray = NCollection_Array1<int>;
using PointArray = NCollection_Array1<gp_Pnt>;
using PointGrid = NCollection_Array2<gp_Pnt>;
using ShapeList = NCollection_List<TopoDS_Shape>;
using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
using ShapeAncestors = NCollection_IndexedDataMap<TopoDS_Shape, ShapeList, TopTools_ShapeMapHasher>;
#if OCC_VERSION_HEX >= 0x080000
using PointHandleArray = NCollection_HArray1<gp_Pnt>;
using RealHandleArray = NCollection_HArray1<double>;
#else
// Older kernels expose distinct RTTI classes for handle-managed arrays.
using PointHandleArray = TColgp_HArray1OfPnt;
using RealHandleArray = TColStd_HArray1OfReal;
#endif
} // namespace stl2curves::detail
