#include "matar.h"

#ifdef MATAR_WITH_KOKKOS
typedef CArrayKokkos<real_t> MatarRealCArray;
typedef CArrayKokkos<size_t> MatarUIntCArray;
typedef CArrayKokkos<int> MatarIntCArray;
#else
typedef CArray<real_t> MatarRealCArray;
typedef CArray<size_t> MatarUIntCArray;
typedef CArray<int> MatarIntCArray;
#endif // MATAR_WITH_KOKKOS
