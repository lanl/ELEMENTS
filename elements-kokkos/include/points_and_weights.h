#pragma once


#include "MatarTypeDefs.h"
#include "common.h"


/* Equispaced and Chebyshev point distributions */
template <typename NumType>
void equispaced_points(SizeType N, NumType &zl, NumType &zr, NumType *z);

template <typename NumType>
void chebyshev_points(SizeType N, NumType &zl, NumType &zr, NumType *z);


/* Hard-coded Gauss-Lobatto and Gauss-Legendre quadrature points and weights */
void lobatto_points(MatarRealCArray &lob_nodes_1D, const int &num);
void lobatto_weights(MatarRealCArray &lob_weights_1D, const int &num);

void legendre_points(MatarRealCArray &leg_nodes_1D, const int &num);
void legendre_weights(MatarRealCArray &leg_weights_1D, const int &num);
