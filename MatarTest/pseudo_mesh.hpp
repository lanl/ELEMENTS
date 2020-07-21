#ifndef PSEUDO_MESH_H
#define PSEUDO_MESH_H

#include "matar.h"

using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;


class pseudo_mesh {

    public:
        int size1, size2, size3;
	    //int size4, size5, size6;
        //Change this name
        CArrayKokkos <size_t> mystride;
        RaggedDownArrayKokkos <double> var;
	    FMatrixKokkos <int> var1;

        // default constructor
        KOKKOS_FUNCTION
        pseudo_mesh();

        // init function
        KOKKOS_FUNCTION
        void init(int pnts1, int pnts2);

        KOKKOS_FUNCTION
        void init(int pnts1, int pnts2, TeamPolicy::member_type teamMember);

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
