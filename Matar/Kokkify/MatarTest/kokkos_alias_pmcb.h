#ifndef KOKKOS_ALIAS_PMCB_H
#define KOKKOS_ALIAS_PMCB_H

#include <Kokkos_Core.hpp>

using real_t = double;
using u_int  = unsigned int;

using ExecSpace = Kokkos::OpenMP;
using Layout    = Kokkos::LayoutRight;

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







#endif
