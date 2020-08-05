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

#endif
