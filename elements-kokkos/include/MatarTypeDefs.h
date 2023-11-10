#include "common.h"
#include "matar.h"

#ifdef MATAR_WITH_KOKKOS
typedef CArrayKokkos<Real> MatarRealCArray;
typedef CArrayKokkos<SizeType> MatarUIntCArray;
typedef CArrayKokkos<int> MatarIntCArray;
#else
typedef CArray<Real> MatarRealCArray;
typedef CArray<SizeType> MatarUIntCArray;
typedef CArray<int> MatarIntCArray;
#endif // MATAR_WITH_KOKKOS
