#include "pseudo_mesh.hpp"




KOKKOS_FUNCTION
pseudo_mesh::pseudo_mesh() {};

KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2) {

    size1 = pnts1;
    size2 = pnts2;
    // Change this name
    mystride = CArrayKokkos <size_t> (size2);
    for (int i = 0; i < size2; i++) {
        mystride(i) = i+1;
    }
    var = RaggedDownArrayKokkos <double> (&mystride(0), size2);
    var1 = FMatrixKokkos<int>(size1,size2);
}

KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2, TeamPolicy::member_type teamMember) {

    size1 = pnts1;
    size2 = pnts2;
    // Change this name
    mystride = CArrayKokkos <size_t> (size2);
    Kokkos::parallel_for (Kokkos::TeamThreadRange (teamMember, size2), [=] (const int i) {
        mystride(i) = i+1;
        printf("%d\n", mystride(i));
        });
    teamMember.team_barrier();
    var = RaggedDownArrayKokkos <double> (&mystride(0), size2);
    var1 = FMatrixKokkos<int>(size1,size2);
}

/*
//3D
KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2, int pnts3) {
	size1 = pnts1;
	size2 = pnts2;
	size3 = pnts3;
	var = FArrayKokkos<double>(size1, size2, size3);
	var1 = CMatrixKokkos<int>(size1,size2,size3);
    var2 = FMatrixKokkos<double>(size1, size2, size3);
    var3 = CArrayKokkos<double>(size1, size2, size3);
}

//4D
KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2, int pnts3, int pnts4) {
	size1 = pnts1;
	size2 = pnts2;
	size3 = pnts3;
	size4 = pnts4;
	var = FArrayKokkos<double>(size1,size2,size3,size4);
	var1 = CMatrixKokkos<int>(size1,size2,size3,size4);
    var2 = FMatrixKokkos<double>(size1, size2, size3, size4);
    var3 = CArrayKokkos<double>(size1, size2, size3, size4);
}

//5D
KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2, int pnts3, int pnts4, int pnts5) {
	size1 = pnts1;
	size2 = pnts2;
	size3 = pnts3;
	size4 = pnts4;
	size5 = pnts5;
	var = FArrayKokkos<double>(size1,size2,size3,size4,size5);
	var1 = CMatrixKokkos<int>(size1,size2,size3,size4,size5);
    var2 = FMatrixKokkos<double>(size1, size2, size3, size4, size5);
    var3 = CArrayKokkos<double>(size1, size2, size3, size4, size5);
}
	
//6D
KOKKOS_FUNCTION
void pseudo_mesh::init(int pnts1, int pnts2, int pnts3, int pnts4, int pnts5, int pnts6) {
	size1 = pnts1;
	size2 = pnts2;
	size3 = pnts3;
	size4 = pnts4;
	size5 = pnts5;
	size6 = pnts6;
	var = FArrayKokkos<double>(size1,size2,size3,size4,size5,size6);
	var1 = CMatrixKokkos<int>(size1,size2,size3,size4,size5,size6);
    var2 = FMatrixKokkos<double>(size1, size2, size3, size4, size5, size6);
    var3 = CArrayKokkos<double>(size1, size2, size3, size4, size5, size6);
}
*/















