#ifndef GEOMETRY_H
#define GEOMETRY_H

// This pulls in kokkos, matar, mesh, ref_elem stuff, and PT-Scotch
#include "ELEMENTS.h"

#include "matar.h"
#include "shapes.h"


/////////////////////////////////////////////////////////////////////////////
///
/// \fn jacobian
///
/// \brief Calculates the jacobian matrix in 2D, 3D, 4D, ...
///
/// \param jacobian     The jacobian matrix calculated in the routine, J[dims,dims]
/// \param node_coords  An array containing all node coords on the mesh
/// \param nodes_in_an_elem  The node indices in a single elem, its a 1D array
/// \param a_grad_basis The gradient of the basis at a single point, Grad[DOFs,dims]
///
/////////////////////////////////////////////////////////////////////////////
template <typename T1, typename T2, typename T3, typename T4>
KOKKOS_INLINE_FUNCTION
void jacobian(
    const T1 &jacobian,         // e.g., ViewCArrayKokkos <double>
    const T2 &node_coords,      // e.g., DCArrayKokkos    <double>
    const T3 &nodes_in_an_elem, // e.g., ViewCArrayKokkos <size_t>
    const T4 &a_grad_basis){    // e.g., ViewCArrayKokkos <double>

    const size_t dims = a_grad_basis.dims(1);
    const size_t num_dofs_in_elem = nodes_in_an_elem.size();

    // setting jacobian matrix to all zeros
    for(size_t i = 0; i < dims; i++)  // looping over dimension
    for(size_t j = 0; j < dims; j++){ // looping over dimension
        jacobian(i, j) = 0.0;
    } // end for 

    // Calculate Jacobian: J[i,j] = partial x_i/partial \xi_j
    for(size_t i = 0; i < dims; i++) // looping over dimension (partial)
    for(size_t j = 0; j < dims; j++) // looping over dimension (node position)
    for(size_t node_lid = 0; node_lid < num_dofs_in_elem; node_lid++){ 
        const size_t node_gid = nodes_in_an_elem(node_lid);  
        jacobian(i, j) += node_coords(node_gid, i)*a_grad_basis(node_lid, j);           
    } // end for  


} // end of jacobian function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_quadrature_point_connectivity
///
/// \brief Using mesh coordinates, mesh connectivity data structures, and  
///        reference element to build connectivity between quadrature points
///
/// This routine finds the corresponding local index for the face-adjacent
/// quadrature point on the surface.  This map is critical to e.g.
/// discontinous Galerkin or cell-centered finite volume methods that solve a 
/// Riemann problem on the surface at quadrature points.
///       nbr_qpt_lid = surf_qpt_qpt_map(surf_gid,face_lid,qpt_lid) 
/// where surf_gid=[0:num_surfs], face_lid=[0,1], and qpt_lid=[0:num_surf_qpts]
/// so
///       nbr_qpt_lid = surf_qpt_qpt_map(surf_gid,0,qpt_lid)     
///       qpt_lid     = surf_qpt_qpt_map(surf_gid,1,nbr_qpt_lid)  
///
///
/// \param Mesh The unstructured mesh object
/// \param RefSurface& The reference surface object
/// \param surf_qpt_qpt_map Array holding the quadrature point surface map
/// \param node_coords The mesh nodal coordinates
///
/// \return void
///
/////////////////////////////////////////////////////////////////////////////
void build_quadrature_point_connectivity(const swage::Mesh_t& Mesh,
                                         const elements::ReferenceSurface_t& RefSurf,
                                         CArrayKokkos<int>& surf_qpt_qpt_map,
                                         const DCArrayKokkos<double>& node_coords){

    const size_t elem_dims = RefSurf.elem_dims;
    const size_t num_surfs = surf_qpt_qpt_map.dims(0);
    const size_t num_nodes_in_elem = Mesh.num_nodes_in_elem;
    const size_t num_qpts_in_surf = RefSurf.qpt_basis.dims(1);

    FOR_ALL(surf_gid, 0, num_surfs, {
        
        // get the first elem id and face in this surf
        const size_t elem_gid = Mesh.elems_in_surf(surf_gid, 0);
        const size_t face_lid = Mesh.faces_in_surf(surf_gid, 0);

        const size_t num_elems_in_surf = Mesh.num_elems_in_surf(surf_gid);

        // get the neighbor, where on the bdys, we use the first elem info
        size_t nbr_elem_gid = elem_gid;
        size_t nbr_face_lid = face_lid;
        if(num_elems_in_surf==2){
            nbr_elem_gid = Mesh.elems_in_surf(surf_gid, 1); // second elem
            nbr_face_lid = Mesh.faces_in_surf(surf_gid, 1); // second elem
        }    

        ViewCArrayKokkos<size_t> nodes_in_elem(&Mesh.nodes_in_elem(elem_gid,0), num_nodes_in_elem);
        ViewCArrayKokkos<size_t> nodes_in_nbr_elem(&Mesh.nodes_in_elem(nbr_elem_gid,0), num_nodes_in_elem);

        // loop the quadrature points on side_lid 0 and match them to side_lid 1
        for(size_t qpt_lid=0; qpt_lid<num_qpts_in_surf; qpt_lid++){
            
            if (surf_qpt_qpt_map(surf_gid,0,qpt_lid)>=0) continue;  // this qpt was tagged

            ViewCArrayKokkos<double> a_basis(&RefSurf.qpt_basis(face_lid,qpt_lid,0), num_nodes_in_elem);

            double qpt_x[3];
            
            for (size_t dim=0; dim<elem_dims; dim++){
               qpt_x[dim] = 0.0;
            }

            for (size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                const size_t node_gid = nodes_in_elem(node_lid);
                for (size_t dim=0; dim<elem_dims; dim++){
                    qpt_x[dim] += a_basis(node_lid)*node_coords(node_gid,dim);
                }
            }            

            bool found=false;

            // loop over the quadrature points on the other side, slide_lid=1
            for (size_t nbr_qpt_lid=0; nbr_qpt_lid<num_qpts_in_surf; nbr_qpt_lid++){

                if (surf_qpt_qpt_map(surf_gid,1,nbr_qpt_lid)>=0) continue; // this nbr_qpt was tagged

                ViewCArrayKokkos<double> nbr_a_basis(&RefSurf.qpt_basis(nbr_face_lid,nbr_qpt_lid,0), num_nodes_in_elem);

                double nbr_qpt_x[3];
                
                for(size_t dim=0; dim<elem_dims; dim++){
                    nbr_qpt_x[dim] = 0.0;
                }

                for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                    const size_t nbr_node_gid = nodes_in_nbr_elem(node_lid);
                    for(size_t dim=0; dim<elem_dims; dim++){
                        nbr_qpt_x[dim] += nbr_a_basis(node_lid)*node_coords(nbr_node_gid,dim);
                    }
                } 

                double diff_sqrd = 0.0;
                for(size_t dim=0; dim<elem_dims; dim++){
                    diff_sqrd += (nbr_qpt_x[dim]-qpt_x[dim])*(nbr_qpt_x[dim]-qpt_x[dim]);
                }

                if(diff_sqrd<1.e-16){
                    // It is a match!  We found your friend!
                    surf_qpt_qpt_map(surf_gid,0,qpt_lid)     = nbr_qpt_lid;
                    surf_qpt_qpt_map(surf_gid,1,nbr_qpt_lid) = qpt_lid;
                    found = true;
                    break; // go to the next qpt on surf_gid 
                }

            } // end loop over neighbors qpts, that is side_lid 1

            if(found==false) Kokkos::abort("ERROR: Failed to find the matching quadrture point \n");

        } // end loop over qpts on side_lid 0
        
    }); // end parallel for

} // end function

#endif // GEOMETRY_H