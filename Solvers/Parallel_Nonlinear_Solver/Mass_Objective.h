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
#include "Parallel_Nonlinear_Solver.h"

class MassObjective_TopOpt : public ROL::Objective_SimOpt<real_t> {
  
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::MultiVector<real_t, LO, GO, Node> MV;
  typedef ROL::Vector<real_t> V;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,Node> ROL_MV;  

private:

  ROL::Ptr<Parallel_Nonlinear_Solver> FEM_;
  ROL::Ptr<MV> Element_Masses;
  ROL::Ptr<ROL_MV> ROL_Element_Masses;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;

  MassObjective_TopOpt(ROL::Ptr<Parallel_Nonlinear_Solver> FEM, bool nodal_density_flag) 
    : FEM_(FEM), useLC_(true) {
      nodal_density_flag_ = nodal_density_flag;
      Element_Masses = ROL::makePtr<MV>(FEM_->element_map,1,true);
      ROL_Element_Masses = ROL::makePtr<ROL_MV>(Element_Masses);
  }

  real_t value( const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
    real_t c = 0.0;

    host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
   
    //local data count
    size_t nlocal_densities = design_densities.extent(0);

    //set per element mass
    if(nodal_density_flag_){

    }
    else{
      ROL::Ptr<ROL_MV> ROL_Global_Element_Volumes = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Volumes);
    }

    //sum per element results across all MPI ranks
    Elementwise::ReductionSum<Real> sumreduc;
    c = ROL_Element_Masses->reduce(sumreduc);

    return c;
  }

  void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    /*
    // Unwrap g
    ROL::Ptr<MV> gp = getVector(g);
    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);

    // Apply Jacobian
    MV KU(up->size(),0.0);
    if ( useLC_ ) {
      FEM_->build_force(KU);
      // Apply jacobian to u
      for (size_t i=0; i<up->size(); i++) {
        (*gp)[i] = KU[i];
      }
    }
    else {
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp);
      // Apply jacobian to u
      for (size_t i=0; i<up->size(); i++) {
        (*gp)[i] = 2.0*KU[i];
      }
    }
    */
  }

  void gradient_2( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    /*
    // Unwrap g
    ROL::Ptr<MV> gp = getVector(g);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);

    // Apply Jacobian
    g.zero();
    if ( !useLC_ ) {
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_adjoint_jacobian(*gp,U,*zp,U);
    }
    */
  }
  /*
  void hessVec_11( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      MV KV(vp->size(),0.0);
      MV V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_jacobian(KV,V,*zp);
      for (size_t i=0; i<vp->size(); i++) {
        (*hvp)[i] = 2.0*KV[i];
      }
    }
    
  }

  void hessVec_12( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    
    // Unwrap hv
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      MV KU(up->size(),0.0);
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp,*vp);
      for (size_t i=0; i<up->size(); i++) {
        (*hvp)[i] = 2.0*KU[i];
      }
    }
    
  }

  void hessVec_21( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
                     
    // Unwrap g
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
 
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      std::MV<real_t> U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      std::MV<real_t> V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,V);
      for (size_t i=0; i<hvp->size(); i++) {
        (*hvp)[i] *= 2.0;
      }
    }
    
  }

  void hessVec_22( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
                     
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
    
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      MV V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,*vp,U);
    }
    
  }
  */
};

#endif // end header guard
