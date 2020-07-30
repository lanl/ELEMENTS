#ifndef PSEUDO_MESH_H
#define PSEUDO_MESH_H

#include "matar.h"

using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;


class pseudo_mesh {

    public:
        int size1, size2, size3, size4, size5, size6;
        // Simple Matar
        CArrayKokkos  <real_t>  carray;
        CMatrixKokkos <real_t>  cmatrix;
        FArrayKokkos  <real_t>  farray;
        FMatrixKokkos <real_t>  fmatrix;
        // Ragged (Static)
        CArrayKokkos           <size_t> mystride;
        RaggedRightArrayKokkos <real_t> raggedright;
        //RaggedDownArrayKokkos  <real_t> raggeddown;

        // default constructor
        //KOKKOS_FUNCTION
        pseudo_mesh();

        // init function
        //KOKKOS_FUNCTION
        void init(int pnts1, int pnts2);

        //KOKKOS_FUNCTION
        //void init(int pnts1, int pnts2, TeamPolicy::member_type teamMember);

	    //~~~adding different init functions to test dimensions of FArray
	
        /*
	//3D
	KOKKOS_FUNCTION
	void init(int pnts1, int pnts2, int pnts3);

	    //4D mesh
	    KOKKOS_FUNCTION
	    void init(int pnts1, int pnts2, int pnts3, int pnts4);

	    //5D mesh
	    KOKKOS_FUNCTION
	    void init(int pnts1, int pnts2, int pnts3, int pnts4, int pnts5);

	//6D mesh
	KOKKOS_FUNCTION
	void init(int pnts1, int pnts2, int pnts3, int pnts4, int pnts5, int pnts6);
    */

        // destructor
        KOKKOS_FUNCTION
        ~pseudo_mesh() {};

};




#endif
