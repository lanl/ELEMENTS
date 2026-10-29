#ifndef ELEMENTS_PRECISION_UTILS_H
#define ELEMENTS_PRECISION_UTILS_H

// ELEMENTS stores and computes all floating-point data as MATAR's real_t,
// whose precision (double or float) is fixed at configure time by
// ELEMENTS_REAL. Fixed tolerances such as 1e-12 are meaningless in single
// precision, so wrap them in real_tol().

#include "matar.h"

// Machine epsilon of the working precision (device-safe)
KOKKOS_INLINE_FUNCTION constexpr real_t real_eps()
{
    return Kokkos::Experimental::epsilon_v<real_t>;
}

// A tolerance that is exactly `tol` when real_t is double, but never drops
// below `ulps` machine epsilons, so it stays usable when real_t is float.
KOKKOS_INLINE_FUNCTION constexpr real_t real_tol(const real_t tol, const real_t ulps = 100)
{
    return tol > ulps * real_eps() ? tol : ulps * real_eps();
}

#endif // ELEMENTS_PRECISION_UTILS_H
