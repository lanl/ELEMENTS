#ifndef KOKKOS_ALIAS_H
#define KOKKOS_ALIAS_H

#include <stdlib.h> 
#include "materials.h"
#include <Kokkos_Core.hpp>

//MACROS to make the code less scary
#define kmalloc(size) ( Kokkos::kokkos_malloc<CudaMemSpace>(size) )
#define kfree(pnt)        (  Kokkos::kokkos_free(pnt) ) 
//#define parallel_for      ( Kokkos::parallel_for )
//#define parallel_reduce      ( Kokkos::parallel_reduce )
//#define deep_copy           ( Kokkos::deep_copy )
#define ProfileRegionStart  ( Kokkos::Profiling::pushRegion )
#define ProfileRegionEnd  ( Kokkos::Profiling::popRegion )
#define HostMirror        ( Kokkos::create_mirror_view )

using real_t = double;
using u_int  = unsigned int;

using CudaMemSpace    = Kokkos::CudaSpace;
using UVMMemSpace     = Kokkos::CudaUVMSpace;
using ExecSpace       = Kokkos::Cuda;
//using range_policy    = Kokkos::RangePolicy<ExecSpace>;
using Layout          = Kokkos::LayoutLeft;
using mdrange_policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
using mdrange_policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

using RMatrix1D    = Kokkos::View<real_t *,Layout,ExecSpace>;
using RMatrix2D    = Kokkos::View<real_t **,Layout,ExecSpace>;
using RMatrix3D    = Kokkos::View<real_t ***,Layout,ExecSpace>;
using RMatrix4D    = Kokkos::View<real_t ****,Layout,ExecSpace>;
using RMatrix5D    = Kokkos::View<real_t *****,Layout,ExecSpace>;
using IMatrix1D    = Kokkos::View<int *,Layout,ExecSpace>;
using IMatrix2D    = Kokkos::View<int **,Layout,ExecSpace>;
using IMatrix3D    = Kokkos::View<int ***,Layout,ExecSpace>;
using IMatrix4D    = Kokkos::View<int ****,Layout,ExecSpace>;
using IMatrix5D    = Kokkos::View<int *****,Layout,ExecSpace>;
using SVar         = Kokkos::View<size_t,Layout,ExecSpace>;
using SArray1D     = Kokkos::View<size_t *,Layout,ExecSpace>;
using SArray2D     = Kokkos::View<size_t **,Layout,ExecSpace>;
using SArray3D     = Kokkos::View<size_t ***,Layout,ExecSpace>;
using SArray4D     = Kokkos::View<size_t ****,Layout,ExecSpace>;
using SArray5D     = Kokkos::View<size_t *****,Layout,ExecSpace>;

using SHArray1D     = Kokkos::View<size_t *,Layout,Kokkos::HostSpace>;

using Material1D     = Kokkos::View<material_models*,Layout,ExecSpace>;
using MaterialHost1D = Kokkos::View<material_models*,Layout,Kokkos::HostSpace>;

#endif
