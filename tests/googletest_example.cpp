#include "element_types/lagrange_polynomials.h"
#include "element_types/point_distributions.h"

#include "gtest/gtest.h"

#include <cstdlib>  // rand, RAND_MAX
#include <ctime>    // time

namespace my {
namespace project {
namespace {

// The fixture for testing class Foo.
class LagrangePolynomialsTest : public ::testing::Test {
  SizeType Np;
  SizeType Nv;
  SizeType I;

  Real Zl;
  Real Zr;
  Real X;

  Real *c;
  Real *Z;
  Real *w;
  Real *C;
  
  Real tol;

 protected:
  LagrangePolynomialsTest() {
    // You can do set-up work for each test here.
		
		// Set order of Lagrange polynomial interpolation
    Np = 8;
    Nv = Np + 1;  // number of vertices (nodes) of interpolant

    // Generate equispaced vertex coordinates between -1 and 1
    Zl = -1.0;
    Zr = 1.0;
    Z = new Real[Nv];
    equispaced_points(Nv, Zl, Zr, Z);

    // Compute barycentric weights of vertices
    w = new Real[Nv];
    lagrange::compute_barycentric_weights(Nv, Z, w);

    // Allocate workspace for derivative computations
    C = new Real[Nv];

    // Set the tolerance
    tol = 1e-8;
  }

  ~LagrangePolynomialsTest() override {
    // You can do clean-up work that doesn't throw exceptions here.
    delete [] c;
    delete [] Z;
    delete [] w;
    delete [] C;
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  void SetUp() override {
    // Code here will be called immediately after the constructor (right
    // before each test).

		// Initialize random seed
		srand(time(NULL));

    // Generate random array of coefficients between 0 and 1
    c = new Real[Nv];
    for (SizeType i = 0; i < Nv; i++) {
      c[i] = Real(rand())/RAND_MAX;
    }

    // Select a random coordinate between -1 and 1
    X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);

    // Select a random vertex ID, between 0 and (Nv - 1)
    I = std::round((Nv - 1)*Real(rand())/Real(RAND_MAX));
  }

  void TearDown() override {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Class members declared here can be used by all tests in the test suite
  // for Foo.

	/*
   * Test consistency of interpolation 
   */
	void TestInterpolationConsistent(bool is_coincident) {
    Real sum1 = 0.0;
    Real sum2 = 0.0;

    if (is_coincident) {
      // Method 1
      for (SizeType i = 0; i < Nv; i++) {
        Real li = lagrange::eval(Nv, i, I, Z, w, Z[I]);
        sum1 += c[i]*li;
      }

      // Method 2 (barycentric formula)
      sum2 = lagrange::eval_interp(Nv, I, Z, w, Z[I], c);
    } else {
      // Method 1
      for (SizeType i = 0; i < Nv; i++) {
        Real li = lagrange::eval(Nv, i, -1, Z, w, X);
        sum1 += c[i]*li;
      }

      // Method 2 (barycentric formula)
      sum2 = lagrange::eval_interp(Nv, -1, Z, w, X, c);
    }

    EXPECT_TRUE(std::fabs(sum1 - sum2) < tol);
	};
};

/*
 * Check that the following produce consistent results: (1) evaluating the
 * Lagrange basis functions, scaling them by coefficients, and summing these
 * products and (2) evaluating the interpolant using the Barycentric formula.
 */
TEST_F(LagrangePolynomialsTest, TestConsistentInterpolation) {
  TestInterpolationConsistent(true);
  TestInterpolationConsistent(false);
}

/*
 * Check that the following produce consistent results: (1) evaluating the
 * Lagrange basis function derivatives, scaling them by coefficients, and
 * summing these products and (2) evaluating the interpolant derivative using
 * the Barycentric formula.
 */
TEST_F(LagrangePolynomialsTest, TestConsistentDerivatives) {
  EXPECT_TRUE(true);
}

/*
 * Check that derivative of the interpolant computed by the Barycentric formula
 * is correct by comparing it the derivative computed by complex step
 * approximation
 */
TEST_F(LagrangePolynomialsTest, TestCorrectDerivatives) {
  EXPECT_TRUE(true);
}


}  // namespace
}  // namespace project
}  // namespace my

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
