#ifndef ELEMENTS_PRECISION_UTILS_H
#define ELEMENTS_PRECISION_UTILS_H

// ELEMENTS stores and computes all floating-point data as MATAR's real_t,
// whose precision (double or float) is fixed at configure time by
// ELEMENTS_REAL. Fixed tolerances such as 1e-12 are meaningless in single
// precision, so wrap them in real_tol(). For a tolerance on a squared quantity
// (e.g. a squared distance), apply real_tol to the unsquared value and square
// the result: `const real_t tol = real_tol(d); ... < tol*tol`.
// For machine epsilon use mtr::prec::epsilon<real_t>().
//
// Casts to real_t are always written C-style, (real_t)value, never as
// real_t(value).

#include <type_traits>
#include "matar.h"

// A tolerance that is exactly `tol` when real_t is double, and never below
// `ulps` machine epsilons otherwise, so it stays usable when real_t is float.
KOKKOS_INLINE_FUNCTION real_t real_tol(const real_t tol, const real_t ulps = 100)
{
    if constexpr (std::is_same_v<real_t, double>) {
        return tol;
    } else {
        return mtr::prec::fmax(tol, ulps * mtr::prec::epsilon<real_t>());
    }
}

#endif // ELEMENTS_PRECISION_UTILS_H
