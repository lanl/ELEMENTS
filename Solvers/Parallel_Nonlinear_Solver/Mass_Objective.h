#ifndef MASS_OBJECTIVE_TOPOPT_H
#define MASS_OBJECTIVE_TOPOPT_H

#include "utilities.h"
#include "../Solver.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_ParameterList.hpp"

template<class Real>
class MassObjective_TopOpt : public ROL::Objective_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;  

private:
  ROL::Ptr<FEM<Real> > FEM_;

  //forward declare nonlinear parallel solver object to use linear system solve method in objective

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  bool constant_element_density;

  Objective_TopOpt(ROL::Ptr<FEM<Real> > FEM) 
    : FEM_(FEM), useLC_(true) {}

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
 
    // Apply Jacobian
    vector KU(up->size(),0.0);
    if ( useLC_ ) {
      FEM_->build_force(KU);
    }
    else {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp);
    }
    // Compliance
    Real c = 0.0;
    for (uint i=0; i<up->size(); i++) {
      c += (*up)[i]*KU[i];
    }
    return c;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<vector> gp = getVector(g);
    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    vector KU(up->size(),0.0);
    if ( useLC_ ) {
      FEM_->build_force(KU);
      // Apply jacobian to u
      for (uint i=0; i<up->size(); i++) {
        (*gp)[i] = KU[i];
      }
    }
    else {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp);
      // Apply jacobian to u
      for (uint i=0; i<up->size(); i++) {
        (*gp)[i] = 2.0*KU[i];
      }
    }
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<vector> gp = getVector(g);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    g.zero();
    if ( !useLC_ ) {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_adjoint_jacobian(*gp,U,*zp,U);
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector KV(vp->size(),0.0);
      vector V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_jacobian(KV,V,*zp);
      for (uint i=0; i<vp->size(); i++) {
        (*hvp)[i] = 2.0*KV[i];
      }
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap hv
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector KU(up->size(),0.0);
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp,*vp);
      for (uint i=0; i<up->size(); i++) {
        (*hvp)[i] = 2.0*KU[i];
      }
    }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
 
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      std::vector<Real> U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      std::vector<Real> V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,V);
      for (uint i=0; i<hvp->size(); i++) {
        (*hvp)[i] *= 2.0;
      }
    }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
    
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      vector V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,*vp,U);
    }
  }
};
