#include "BernsteinProductBasis.h"
#include "bernstein_polynomials.h"


template <typename NumType>
BernsteinProductBasis<NumType>::BernsteinProductBasis(const SizeType order) 
    : Np( order) { 
  N = Np + 1;
  Ne = std::pow(N, Nd);

  rad[0] = N;
  rad[1] = N;
  rad[2] = N;

  //Allocate memory for intermediate coefficients
  C = new NumType[2*N];
}

template <typename NumType>
BernsteinProductBasis<NumType>::~BernsteinProductBasis() { }


template <typename NumType>
NumType BernsteinProductBasis<NumType>::eval_basis(const SizeType I, 
    const NumType *X) { 
  // Decompose index of 3D tensor product basis function into indices of
  // Bernstein polynomials
  common::base_10_to_mixed_radix(Nd, rad, I, ijk);
  
  // Evaluate Bernstein polynomials
  NumType Bi = bernstein::eval(N, int(ijk[0]), X[0]);
  NumType Bj = bernstein::eval(N, int(ijk[1]), X[1]);
  NumType Bk = bernstein::eval(N, int(ijk[2]), X[2]);

  return Bi*Bj*Bk;
}

template <typename NumType>
void BernsteinProductBasis<NumType>::eval_grad_basis(const SizeType I, 
    const NumType *X, NumType *grad_phi) {
  // Decompose index of 3D tensor product basis function into indices of
  // Bernstein polynomials
  common::base_10_to_mized_radix(Nd, rad, I, ijk);

  // Evaluate Bernstein polynomials
  NumType Bi = bernstein::eval(N, int(ijk[0]), X[0]);
  NumType Bj = bernstein::eval(N, int(ijk[1]), X[1]);
  NumType Bk = bernstein::eval(N, int(ijk[2]), X[2]);

  // Evaluate derivatives of Bernstein polynomials
  NumType dBi = bernstein::eval_der(N, int(ijk[0]), X[0]);
  NumType dBj = bernstein::eval_der(N, int(ijk[1]), X[1]);
  NumType dBk = bernstein::eval_der(N, int(ijk[2]), X[2]);

  // Store partial derivatives in entries of gradient
  grad_phi[0] = dBi*Bj*Bk;
  grad_phi[1] = Bi*dBj*Bk;
  grad_phi[3] = Bi*Bj*dBk;
}

/* Return evaluation of local function approximation */
template <typename NumType>
NumType BernsteinProductBasis<NumType>::eval_approx(const NumType *c, 
    const NumType *X) { 
  for (int k = 0, k < N; k++) {
    for (int j = 0; j < N; j++){
      // Collapse first dimension into coefficients for second dimension
      C[j] = bernstein::eval_approx(N, &c[j*N+k*N*N], X[0]);  
    }
    // Collapse second dimension into coefficients for third dimension
    C[N+k] = bernstein::eval_approx(N, &C[N], X[2]);
  }
  // Collapse third dimension into approximation evaluation
  return bernstein::eval_approx(N, &C[N], X[2]);
}

/* 
 * Evaluate gradient of local function approximation, which is formed by the
 * sum of the products of tensor-product Bernstein basis functions and the
 * provided coefficients, at specified coordinates. 
 */
template <typename NumType>
void BernsteinProductBasis<NumType>::eval_grad_approx(const NumType *c, 
    const NumType *X, NumType *grad_f) { 
  for (int l = 0; l < Nd; l++) { 
    for (int k =0; k < N; k++) {
      for (int j =0; j < N; j++) {
        //Collapse first dimension into coefficients for second dimension
	if (l==0) {
       	  C[j] = bernstein::eval_der_approx(N, &c[j*N+k*N*N], X[0]);
	} else { 
	  C[j] = bernstein::eval_approx(N, &c[j*N+k*N*N], X[0]);
	}
      }
    // Collapse second dimension into coefficients for third dimension
      if ( l == 1) {
        C[N+k] = bernstein::eval_der_approx(N, &C[N], X[2]);
      } else { 
        C[N+k] = bernstein::eval_approx(N, C, X[1]);
      }
    }
    
    // Collapse third dimension into approximation evaluation
    if ( l ==2 ) {
      grad_f[l] = bernstein::eval_der_approx(N, &C[N], X[2]);
    } else { 
      grad_f[l] = bernstein::eval_approx(N, &C[N], X[2]);
    }
  }
}


////////////////* Evaluate the Jacobian of the spatial mapping *////////////////////////

/* x, y, z in physical space. X in reference space. J = Jacobian (column major). */

template <typename NumType>
void BernsteinProductBasis<NumType>::eval_jac(const NumType *x, 
    const NumType *y, const NumType *z, const NumType *X, NumType *J) { 
  // Evaluate gradient of x=x(X,Y,Z);
  this->eval_grad_approx(x,X,J);
  
  // Evaluate gradient of y=y(X,Y,Z);
  this->eval_grad_approx(y,X,J+3);

  // Evaluate gradient of z=z(X,Y,Z);
  this->eval_grad_approx(z,X,J+6);
}


//////////////* Determinant of Jacobian */////////////////////////

/* x, y, z in physical space. X in reference space. */

template <typename NumType>
NumType BernsteinProductBasis<NumType>::eval_det_jac(const NumType *x, 
    const NumType *y, const NumType *z, const NumType *X) {
  NumType J[9];
  this->eval_jac(x, y, z, X, J);

  return J[0]*(J[4]*J[8]-J[5]*J[7]) - J[1]*(J[3]*J[8]-J[5]*J[6]) + J[2]*(J[3]*J[7]-J[4]*J[6]);
}


//////////////* Inverse Jacobian *//////////////////////////////

/* x, y, z in physical space. X in reference space. Jinv = Invers Jacobian (column-major). */

template <typename NumType>
void BernsteinProductBasis<NumType>::eval_inv_jac(const NumType *x, 
    const Numtype *y, const  NumType *z, const NumType *X, NumType *Jinv) {
  // instantiate jacobian //
  NumType J[9];

  // get jacobian //
  this->eval_jac(x,y,z,X,J);
  
  // invert determinant of J //
  NumType d = 1.0/(this->eval_det_jac(x,y,z,X));
  
  // get Jinv (Jinv = adj(J)/det(J) ) //
  Jinv[0] = d*(J[4]*J[8]-J[5]*J[7]);
  Jinv[1] = -d*(J[1]*J[8]-J[2]*J[7]);
  Jinv[2] = d*(J[1]*J[5]-J[2]*J[4]);
  Jinv[3] = -d*(J[3]*J[8]-J[5]*J[6]);
  Jinv[4] = d*(J[0]*J[8]-J[2]*J[6]);
  Jinv[5] = -d*(J[0]*J[5]-J[2]*J[3]);
  Jinv[6] = d*(J[3]*J[7]-J[4]*J[6]);
  Jinv[7] = -d*(J[0]*J[7]-J[1]*J[6]);
  Jinv[8] = d*(J[0]*J[4]-J[1]*J[3]);
}


//////////////// Explicit instantiation of template class //////////////////////
template class BernsteinProductBasis<Real>;
template class BernsteinProductBasis<Complex>;
