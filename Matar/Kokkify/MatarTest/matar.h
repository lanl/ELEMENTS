#ifndef MATAR_H
#define MATAR_H 


//==============================================================================
//   MATAR > 'MAT'rix & 'AR'rays
//
//     This code is for creating and viewing matricies and arrays using
//     standard notation e.g., A(i,j,k).  The indicies for a matrix
//     start at 1 for consistancy with mathematics.  The indicies for
//     an array start at 0 for consistancy with the C++ langauge.  The
//     entries in the arrays and matrices are stored in a 1D array that
//     is contiguous in memory, but are accessed using e.g., A(i,j,k).
//     The classes that have f_array in the name stride through the 1D  
//     array following the Fortran convention (i.e., first index varies
//     the quickest).  The classes with c_array in the name stride through
//     the 1D array following the C convention (i.e., last index varies
//     the quickest).
//
//   view_f_array(array_in, dim1, ...)
//   f_array_t(dim1,...)
//   view_f_matrix(matrix_in, dim1, ...)
//   f_matrix_t(dim1,...)
//
//   view_c_array(array_in, dim1, ...)
//   c_array_t(dim1,...)
//   view_c_matrix(matrix_in, dim1, ...)
//   c_matrix_t(dim1,...)
//
//
//==============================================================================
//
//  matar.cpp
//
//     g++ --std=c++14 matar.cpp
//
//  Created by Nathaniel Morgan on 5/18/19.
//
//==============================================================================

//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <Kokkos_Core.hpp>
#include <assert.h>
#include "kokkos_alias.h"

// To disable asserts, uncomment the following line
//#define NDEBUG


//==============================================================================
//   Fortran stride array related classes  (first index varies the quickest)
//==============================================================================

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class ViewFArray {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    T * this_array;
    
public:
    
    // default constructor
    ViewFArray ();

    //---1D array---
    ViewFArray(T *some_array, size_t some_dim1);
    T& operator()(size_t i);
    
    //--- 2D array ---
    
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2);
    T& operator()(size_t i, size_t j);
    
    //--- 3D array ---
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3);
    T& operator()(size_t i, size_t j, size_t k);
    
    //--- 4D array ---
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4);
    T& operator()(size_t i, size_t j, size_t k, size_t l);

    //--- 5D array ---
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
   
    //------6D -----

    ViewFArray (T *some_array,size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

    
}; // end of viewFArray

//class definitions for viewFArray

//~~~~constructors for viewFArray for 1D to 6D~~~~~~~

//no dimension
template <typename T>
ViewFArray<T>::ViewFArray(){}

//1D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1){
	dim1_ = some_dim1;
	this_array = some_array;
}

//2D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	this_array = some_array;
}

//3D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	this_array = some_array;
}

//4D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	this_array = some_array;
}

//5D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	this_array = some_array;
}

//6D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	dim6_ = some_dim6;
	this_array = some_array;
}

//~~~~~~operator () overload 
//for dimensions 1D to 6D
//indices for array are from 0...N-1

//1D
template <typename T>
T& ViewFArray<T>::operator()(size_t i)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 1D!");
	return this_array[i];
}

//2D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j) 
{
	assert( i < dim1_ && "i is out of bounds in ViewFArray 2D!");
	assert( j < dim2_ && "j is out of bounds in ViewFArray 2D!");
	return this_array[i + j*dim1_];
}

//3D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j,size_t k)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 3D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 3D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 3D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 4D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 4D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 4D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 4D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 5D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 5D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 5D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 5D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArray 5D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
T& ViewFArray<T>:: operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 6D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 6D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 6D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 6D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArray 6D!");
	assert(n < dim6_ && "n is out of bounds in ViewFArray 6D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

//~~~~~~end of class definitions for ViewFArray~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~Begin Class for ViewFArrayKokkos~~~~~~~~~~~~~~~~~~~~~`
template <typename T>
class ViewFArrayKokkos {

private: 
	size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
	T * this_array;

public:

	//default constructor
	KOKKOS_FUNCTION
	ViewFArrayKokkos();


	//~~~test~~~~
	//overload constructor with FArray
	//following overload syntax from RaggedDownArray
	//only doing 2D case
	//KOKKOS_FUNCTION
	//ViewFArrayKokkos(FArrayKokkos <size_t> &some_array, size_t dim1, size_t dim2);

	//overload constructors for 1D up to 6D
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1); //1D
	
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2); //2D
	
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3); //3D
	
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4); //4D
	
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5); //5D
	
	KOKKOS_FUNCTION
	ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6);

	//overload () operator to access elements like (i,j,...)
	//for arrays, indices are from 0,...,N-1
	
	KOKKOS_FUNCTION
	T& operator()(size_t i); //1D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j); //2D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k); //3D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k, size_t l); //4D

    KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m); //5D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n); //6D


}; //end class declarations for ViewFArrayKokkos

//~~~~~~~~~~~~~~~~~~~~ViewFArrayKokkos class definitions~~~~~~~~~~

//constructors from no dimensiton up to 6D

//no dimension
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos() {}

//1D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1) {
	dim1_ = dim1;
	this_array = some_array;
}

//2D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2){
	dim1_ = dim1;
	dim2_ = dim2;
	this_array = some_array;
}

/*
//2D with FArrayKokkos input
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(FArrayKokkos<size_t> &some_array, size_t dim1, size_t dim2) {
	dim1_ = dim1;
	dim2_ = dim2;
	this_array = some_array;
}
*/

//3D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	this_array = some_array;
}

//4D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	this_array = some_array;
}

//5D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	this_array = some_array;
}

//6D
template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	dim6_ = dim6;
	this_array = some_array;
}

//~~~~~~~~definitions for overload () operator~~~~~~~
//indices start are from 0,...,N-1

//1D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i){
	assert( i < dim1_ && "i is out of bounds in ViewFArrayKokkos 1D!");
	return this_array[i];
}

//2D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j){
	assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 2D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 2D!");
	return this_array[i + j*dim1_];
}

//3D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) {
	assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 3D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 3D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 3D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k , size_t l) {
	assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 4D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 4D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 4D!");
	assert(l < dim4_ && "l is out of bounds in VIewFArraykokkos 4D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m){
	assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 5D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 5D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 5D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArrayKokkos 5D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArrayKokkos 5D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) {
	assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 6D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 6D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 6D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArrayKokkos 6D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArrayKokkos 6D!");
	assert(n < dim6_ && "n is out of bounds in ViewFArrayKokkos 6D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

//~~~~~~~~~~~~~~~~~~~~~~end of ViewFArrayKokkos declarations~~~~~~~~~~~~


// create a 1D vector that is accessed as an array(i,...,n), where i=[0:N-1]
template <typename T>
class FArray {
    
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_, length_;
    T * this_array;
    
public:
    
    // default constructor
   FArray ();
   
    //overload constructors from 1D to 6D
     
   FArray(size_t some_dim1);
   FArray(size_t some_dim1, size_t some_dim2);
   FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3);
   FArray(size_t some_dim1, 
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4);
    
   FArray(size_t some_dim1,
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4,
          size_t some_dim5);

   FArray(size_t some_dim1,
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4,
          size_t some_dim5,
          size_t some_dim6);

    // overload operator() to access data as array(i,....,n);
    T& operator()(size_t i);
    T& operator()(size_t i, size_t j);
    T& operator()(size_t i, size_t j, size_t k);
    T& operator()(size_t i, size_t j, size_t k, size_t l);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

    //overload = operator
    FArray& operator=(const FArray& temp);

    // deconstructor
    ~FArray ( );
    
}; // end of f_array_t

//write out functions

//constructors
template <typename T>
FArray<T>::FArray(){}

template <typename T>
FArray<T>::FArray(size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array = new T[length_];
}

template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_array = new T[length_];
}

//3D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_array = new T[length_];
}

//4D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_array = new T[length_];
}

//5D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_array = new T[length_];
}

//6D
template <typename T>
FArray<T>::FArray(size_t some_dim1,size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_array = new T[length_];
}

//overload operator () for 1D to 6D
//indices are from [0:N-1]

//1D
template <typename T>
T& FArray<T>::operator()(size_t i)
{
    assert( i < dim1_ && "i is out of bounds in FArray 1D!");
    return this_array[i];
}

//2D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j)
{
    assert( i < dim1_ && "i is out of bounds in FArray 2D!");
    assert( j < dim2_ && "j is out of bounds in FArray 2D!");
    return this_array[i + j*dim1_];
}

//3D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k)
{
    assert( i < dim1_ && "i is out of bounds in FArray 3D!");
    assert( j < dim2_ && "j is out of bounds in Farray 3D!");
    assert( k < dim3_ && "k is out of bounds in FArray 3D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
    assert( i < dim1_ && "i is out of bounds in FArray 4D!");
    assert( j < dim2_ && "j is out of bounds in FArray 4D!");
    assert( k < dim3_ && "k is out of bounds in FArray 4D!");
    assert( l < dim4_ && "l is out of bounds in FArray 4D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
    assert( i < dim1_ && "i is out of bounds in FArray 5D!");
    assert( j < dim2_ && "j is out of bounds in FArray 5D!");
    assert( k < dim3_ && "k is out of bounds in FArray 5D!");
    assert( l < dim4_ && "l is out of bounds in FArray 5D!");
    assert( m < dim5_ && "m is out of bounds in FArray 5D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{

    assert( i < dim1_ && "i is out of bounds in FArray 5D!");
    assert( j < dim2_ && "j is out of bounds in FArray 5D!");
    assert( k < dim3_ && "k is out of bounds in FArray 5D!");
    assert( l < dim4_ && "l is out of bounds in FArray 5D!");
    assert( m < dim5_ && "m is out of bounds in FArray 5D!");
    assert( n < dim6_ && "n is out of bounds in FArray 6D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

// = operator
//THIS = FArray <> TEMP(n,m,...)
template <typename T>
FArray<T>& FArray<T>::operator= (const FArray& temp)
{
	if(this != & temp) {
	  dim1_ = temp.dim1_;
	  dim2_ = temp.dim2_;
	  dim3_ = temp.dim3_;
	  dim4_ = temp.dim4_;
	  dim5_ = temp.dim5_;
	  dim6_ = temp.dim6_;
	  length_ = temp.length_;
	  this_array = new T[length_];
	}
} 

//delete FArray
template <typename T>
FArray<T>::~FArray(){
    delete [] this_array;
}

//~~~~~~~~~~FArrayKokkos~~~~~~~~~~~~
template <typename T>
class FArrayKokkos {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_, length_;
    T * this_array;

public:

    //default constructor
    KOKKOS_FUNCTION
    FArrayKokkos();

    //overloaded constructors from 1D - 6D
    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1);
    
    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1, size_t some_dim2);

    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3);

    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4);

    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5); 

    KOKKOS_FUNCTION
    FArrayKokkos(size_t some_dim1, size_t sone_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6);

    //overload operator() to acces data
    //from 1D to 6D
    
    KOKKOS_FUNCTION
    T& operator()(size_t i);
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j);

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k);

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l);

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m);

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

    //overload = operator
    KOKKOS_FUNCTION
    FArrayKokkos& operator= (const FArrayKokkos &temp);

    //destructor
    KOKKOS_FUNCTION
    ~FArrayKokkos();    

}; //end of FArrayKokkos declarations

//~~~~~begin FArrayKokkos definitions

//constructors 1D to 6D

//no size
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos() {}

//1D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1){
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array = new T[length_];
}

//2D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_array = new T[length_];
}

//3D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_array = new T[length_];
}

//4D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_array = new T[length_];
}

//5D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_array = new T[length_];
}

//6D
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_array = new T[length_];
}

//definitions of overload operator()
//for 1D to 6D
//note the indices for array all start at 0

//1D
template<typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()( size_t i) {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 1D!");
    return this_array[i];
}

//2D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j){
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 2D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 2D!");
    return this_array[i +j*dim1_];
}

//3D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 3D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 3D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 3D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 4D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 4D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 4D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 4D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 5D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 5D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 5D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 5D!");
    assert( m < dim5_ && "m is out of bounds in FArrayKokkos 5D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 6D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 6D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 6D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 6D!");
    assert( m < dim5_ && "m is out of bounds in FArrayKokkos 6D!");
    assert( n < dim6_ && "n is out of bounds in FArrayKokkos 6D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

//overload = operator
//for object assingment THIS = FArrayKokkos<> TEMP(n,m,,,,)
template <typename T>
FArrayKokkos<T>& FArrayKokkos<T>::operator= (const FArrayKokkos& temp){
	if( this != &temp){
	  dim1_ = temp.dim1_;
	  dim2_ = temp.dim2_;
	  dim3_ = temp.dim3_;
	  dim4_ = temp.dim4_;
	  dim5_ = temp.dim5_;
	  dim6_ = temp.dim6_;
	  length_ = temp.length_;
	  this_array = new T[length_];
	}
    return *this;
}

//destructor
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::~FArrayKokkos() {
    delete[] this_array;
}
//~~~~~~~~~~~~~~~~END OF FArrayKokkos~~~~~~~~

//==============================================================================
//   C stride matrix related classes  (last index varies the quickest)
//==============================================================================

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class ViewCMatrix {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
     T * this_matrix;
		    
public:
		    
    // default constructor
    ViewCMatrix();
		    
		    
    //--- 1D array ---	   	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,size_t some_dim1);
    T& operator() (size_t i);
		    
    //--- 2D array ---	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix, size_t some_dim1, size_t some_dim2);
		    
    T& operator() (size_t i, size_t j);
		    
    //--- 3D array ---	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3);
    T& operator() (size_t i, size_t j, size_t k);
		    
    //--- 4D array ---
		    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3,
		size_t some_dim4);
		    
    T& operator() (size_t i, size_t j, size_t k, size_t l);

		    
    //--- 5D array ---
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3,
		size_t some_dim4,
		size_t some_dim5);
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m);

    //--- 6D array ---		    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		   size_t some_dim1,
		   size_t some_dim2,
		   size_t some_dim3,
		   size_t some_dim4,
		   size_t some_dim5,
		   size_t some_dim6);
		    
   T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

		    
}; // end of ViewCMatrix

//~~~~~~~ ViewCMatrix Class Definitions~~~~~~

//constructors from no dim to 6D

//no dim
template <typename T>
ViewCMatrix<T>::ViewCMatrix(){}

//1D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix,size_t some_dim1) {
	dim1_ = some_dim1;
	this_matrix = some_matrix;
}

//2D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	this_matrix = some_matrix;
}

//3D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	this_matrix = some_matrix;
}

//4D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	this_matrix = some_matrix;
}

//5D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	this_matrix = some_matrix;
}

//6D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	dim6_ = some_dim6;
	this_matrix = some_matrix;
}

//~~~~~~~overload () operator for 1D to 6D
//recall indices are from 1,...,N

//1D
template <typename T>
T& ViewCMatrix<T>:: operator() (size_t i)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 1D!");
	return this_matrix[i-1];
}

//2D
template <typename T>
T& ViewCMatrix<T>::operator() (size_t i, size_t j)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 2D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 2D!");
	return this_matrix[(i-1)*dim2_ + (j-1)];
}

//3D
template <typename T>
T& ViewCMatrix<T>::operator () (size_t i, size_t j, size_t k)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 3D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 3D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 3D!");
	return this_matrix[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 4D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 4D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 4D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 4D!");
	return this_matrix[(l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k,size_t l, size_t m)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 5D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 5D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 5D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 5D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewCMatrix 5D!");
	return this_matrix[(m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 6D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 6D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 6D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 6D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewCMatrix 6D!");
	assert(n < dim6_+1 && "n is out of bounds for ViewCMatrix 6D!");
	return this_matrix[(n-1)+ (m-1)*dim6_ + (l-1)*dim5_*dim6_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim6_*dim4_*dim3_*dim2_];
}

//~~~~~~end of ViewCMatrix Class Definitions~~~~~~~

//~~~~~~~~~begin ViewCMatrixKokkos Class Declarations~~~~~~
template <typename T>
class ViewCMatrixKokkos {

private:
	size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
	T *this_matrix;

public:
	//default constructor
	KOKKOS_FUNCTION
	ViewCMatrixKokkos();

	//overload constructors for 1D up to 6D
	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1); //1D

	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2); //2D

	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3); //3D

	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4); //4D

	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5); //5D

	KOKKOS_FUNCTION
	ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6); //6D

	//overload () operator for 1D up to 6D
	//Note indices are from 1,...,N
	KOKKOS_FUNCTION
	T& operator()(size_t i); //1D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j); //2D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j , size_t k); //3D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k , size_t l); //4D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m); //5D

	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n); //6D

}; //~~~~~end ViewCMatrixKokkos class declarations~~~~~~~

//~~~~~~~~~~~~begin ViewCMatrixKokkos class declarations~~~~~~~~~~

//begin with constructors from no dimension up to 6D

//no dimensions
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(){}

//1D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T *some_matrix, size_t dim1){
	dim1_ = dim1;
	this_matrix = some_matrix;
}

//2D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2) {
	dim1_ = dim1;
	dim2_ = dim1;
	this_matrix = some_matrix;
}

//3D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	this_matrix = some_matrix;
}

//4D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4){
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	this_matrix = some_matrix;
}

//5D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T *some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5) {
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	this_matrix = some_matrix;
}

//6D
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T * some_matrix, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6){
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	dim6_ = dim6;
	this_matrix = some_matrix;
}

//~~~~overload () operator for 1D up to 6D
//Note: Matrix indices are from 1,...,N
//1D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i){
	assert(i < dim1_+1 && "i is out of bounds in ViewCMatrixKokkos1D!");
	return this_matrix[i-1];
}

//2D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j) {
	assert(i < dim1_+1 && "i is out of bounds in ViewCMatrixKokkos 2D!");
	assert(j < dim2_+1 && "j is out of bounds in ViewCMatrixKokkos 2D!");
	return this_matrix[(i-1)*dim2_ + (j-1)];
}

//3D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k){
	assert(i < dim1_+1 && "i is out of bounds in ViewCMatrixKokkos 3D!");
	assert(j < dim2_+1 && "j is out of bounds in ViewCMatrixKokkos 3D!");
	assert(k < dim3_+1 && "k is out of bounds in ViewCMatrixKokkos 3D!");
	return this_matrix[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j , size_t k, size_t l) { 
	assert(i < dim1_+1 && "i is out of bounds in ViewCMatrixKokkos in 4D!");
	assert(j < dim2_+1 && "j is out of bounds in ViewCMatrixKokkos 4D!");
	assert(k < dim3_+1 && "k is out of bounds in ViewCMatrixKokkos 4D!");
	assert(l < dim4_+1 && "l is out of bounds in ViewCMatrixKokkos 4D!");
	return this_matrix[ (l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) {
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrixKokkos 5D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrixKokkos 5D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrixKokkos in 5D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrixKokkos in 5D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewcMatrixKokkos in 5D!");
	return this_matrix[(m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) {
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrixKokkos in 6D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrixKokkos in 6D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrixKokkos in 6D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrixKokkos in 6D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewCMatrixKokkos in 6D!");
	assert(n < dim6_+1 && "n is out of bounds for ViewCMatrixKokkos in 6D!");
	return this_matrix[(n-1) + (m-1)*dim6_ + (l-1)*dim6_*dim5_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim6_*dim5_*dim4_*dim3_*dim2_];
}

//~~~~~~~~end of ViewCMatrixKokkos class definitions ~~~~~~~~~~~


/*
 * ViewFMatrix
*/

// view a 1D vector as an matrix(i,...,n), where i = [1:N]
template <typename T>
class ViewFMatrix {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T * this_matrix_;
    
public:
    
    // Default constructor
    ViewFMatrix ();
    
    //--- 1D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1);
    
    // Overload  operator() to access data as matrix(i,j),
    // where i = [1:N], j = [1:N]
    T& operator()(size_t i) const;
    
    //--- 2D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2);
    
    // Overload operator() to access data as matrix(i,j),
    //  where i=[1:N], j=[1:N]
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3);
    
    // Overload operator() to access data as matrix(i,j,k),
    // where i = [1:N], j = [1:N], k = [1:N]
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4);
        
    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [0:n-1], j = [1:N], k = [1:N], l = [1:N]
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D matrix ---
    
    // Overloaded constructor
    ViewFMatrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5);
        
    // Overload operator() to access data as matrix(i,j,k,l,m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    T& operator() (size_t i, 
                   size_t j, 
                   size_t k, 
                   size_t l, 
                   size_t m) const;
    
    //--- 6D matrix ---
    
    // Overloaded constructor
    ViewFMatrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5,
                 size_t some_dim6);
        
    // Overload operator() to access data as matrix(i,j,k,l,m,n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    T& operator()(size_t i, 
                  size_t j, 
                  size_t k, 
                  size_t l, 
                  size_t m, 
                  size_t n) const;

    size_t size() const;
}; // end of ViewFMatrix

// Default constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix() {}
      
//--- 1D matrix ---
        
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i,j),
// where i = [1:N], j = [1:N]
 template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 1D");  // die if >= dim1
        
    return this_matrix_[(i - 1)];
}
        
//--- 2D matrix ---
        
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i,j),
//  where i=[1:N], j=[1:N]
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j) const
{
       
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 2D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 2D");  // die if >= dim2
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}
        
//--- 3D matrix ---
     
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix (T *some_matrix,
                             size_t some_dim1,
                             size_t some_dim2,
                             size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i,j,k),
// where i = [1:N], j = [1:N], k = [1:N]
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 3D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 3D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 3D");  // die if >= dim3
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}
        
//--- 4D matrix ---
       
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i, j, k, l),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N]
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 4D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 4D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 4D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 4D");  // die if >= dim4
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}
        
//--- 5D matrix ---
        
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i,j,k,l,m),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l, 
                                     size_t m) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 5D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 5D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 5D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 5D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 5D");  // die if >= dim5
       
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

        
//--- 6D matrix ---
       
// Overloaded constructor
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5,
                            size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = some_matrix;
}
        
// Overload operator() to access data as matrix(i,j,k,l,m,n),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 6D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 6D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 6D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 6D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 6D");  // die if >= dim5
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in c_matrix 6D");  // die if >= dim6
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
inline size_t ViewFMatrix<T>::size() const {
    return length_;
}

/*
 * ViewFMatrixKokkos
*/

// view a 1D vector as an matrix(i,...,n), where i = [1:N]
template <typename T>
class ViewFMatrixKokkos {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T * this_matrix_;
    
public:
    
    // Default constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos ();
    
    //--- 1D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos(T *some_matrix,
                      size_t some_dim1);
    
    // Overload  operator() to access data as matrix(i,j),
    // where i = [1:N], j = [1:N]
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;
    
    //--- 2D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos(T *some_matrix,
                      size_t some_dim1,
                      size_t some_dim2);
    
    // Overload operator() to access data as matrix(i,j),
    //  where i=[1:N], j=[1:N]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos(T *some_matrix,
                      size_t some_dim1,
                      size_t some_dim2,
                      size_t some_dim3);
    
    // Overload operator() to access data as matrix(i,j,k),
    // where i = [1:N], j = [1:N], k = [1:N]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos(T *some_matrix,
                      size_t some_dim1,
                      size_t some_dim2,
                      size_t some_dim3,
                      size_t some_dim4);
        
    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [0:n-1], j = [1:N], k = [1:N], l = [1:N]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos (T *some_matrix,
                       size_t some_dim1,
                       size_t some_dim2,
                       size_t some_dim3,
                       size_t some_dim4,
                       size_t some_dim5);
        
    // Overload operator() to access data as matrix(i,j,k,l,m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, 
                   size_t j, 
                   size_t k, 
                   size_t l, 
                   size_t m) const;
    
    //--- 6D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewFMatrixKokkos (T *some_matrix,
                       size_t some_dim1,
                       size_t some_dim2,
                       size_t some_dim3,
                       size_t some_dim4,
                       size_t some_dim5,
                       size_t some_dim6);
        
    // Overload operator() to access data as matrix(i,j,k,l,m,n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    KOKKOS_FUNCTION
    T& operator()(size_t i, 
                  size_t j, 
                  size_t k, 
                  size_t l, 
                  size_t m, 
                  size_t n) const;
    
    KOKKOS_FUNCTION
    size_t size() const;
    
}; // end of ViewFMatrixKokkos

// Default constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos() {}
    
//--- 1D matrix ---
    
// Overloaded constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T *some_matrix,
                                        size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i,j),
// where i = [1:N], j = [1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 1D");  // die if >= dim1
    
    return this_matrix_[(i - 1)];
}
    
//--- 2D matrix ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T *some_matrix,
                                        size_t some_dim1,
                                        size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i,j),
//  where i=[1:N], j=[1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, 
                                           size_t j) const
{
   
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 2D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 2D");  // die if >= dim2
    
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}
    
//--- 3D matrix ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos (T *some_matrix,
                                         size_t some_dim1,
                                         size_t some_dim2,
                                         size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i,j,k),
// where i = [1:N], j = [1:N], k = [1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, 
                                           size_t j, 
                                           size_t k) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 3D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 3D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 3D");  // die if >= dim3
    
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

    
//--- 4D matrix ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T *some_matrix,
                                        size_t some_dim1,
                                        size_t some_dim2,
                                        size_t some_dim3,
                                        size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i, j, k, l),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, 
                                           size_t j, 
                                           size_t k, 
                                           size_t l) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 4D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 4D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 4D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 4D");  // die if >= dim4
    
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

    
//--- 5D matrix ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i,j,k,l,m),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l, 
                                     size_t m) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 5D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 5D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 5D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 5D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 5D");  // die if >= dim5
   
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}
    
//--- 6D matrix ---
   
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5,
                            size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = some_matrix;
}
    
// Overload operator() to access data as matrix(i,j,k,l,m,n),
// where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                    size_t m, size_t n) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 6D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 6D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 6D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 6D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 6D");  // die if >= dim5
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in c_matrix 6D");  // die if >= dim6
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
KOKKOS_FUNCTION
size_t ViewFMatrixKokkos<T>::size() const {
    return length_;
}

/*
 * FMatrix 
 */

template <typename T>
class FMatrix {
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T* this_matrix_;

public:
    // Default constructor
    FMatrix ();

    // --- 1D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1);

    // Overload operator() to access data as matrix(i), where i = [1:N]
    T& operator() (size_t i) const;

    // --- 2D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as matrix(i, j),
    // where i = [1:N], j = [1:N]
    T& operator() (size_t i, size_t j) const;

    // --- 3D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as matrix(i, j, k),
    // where i = [1:N], j = [1:N], k = [1:N]
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as matrix(i, j, k, l, m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N],
    // m = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as matrix(i, j, k, l, m, n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], 
    // m = [1:N], n = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    FMatrix& operator=(const FMatrix& temp);

    size_t size() const;

    // Deconstructor
    ~FMatrix ();

}; // End of FMatrix

template <typename T>
FMatrix<T>::FMatrix() {}

// --- 1D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1) {
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = new T[length_];
}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i) const {
    assert(i >= 1 && i <= dim1_);
    return this_matrix_[i - 1];
}

// --- 2D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = new T[length_];
}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}

// --- 3D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = new T[length_];
}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

// --- 4D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = new T[length_];
}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

// --- 5D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = new T[length_];
}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_) 
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

// --- 6D matrix ---

template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = new T[length_];

}

template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    assert(n >= 1 && n <= dim6_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)  
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)  
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
inline FMatrix<T>& FMatrix<T>::operator= (const FMatrix& temp)
{
    // Do nothing if assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = new T[length_];
    }
    
    return *this;
}

template <typename T>
inline size_t FMatrix<T>::size() const {
    return length_;
}

template <typename T>
FMatrix<T>::~FMatrix() {
    delete[] this_matrix_;
}

/*
 * FMatrixKokkos 
 */

template <typename T>
class FMatrixKokkos {
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T* this_matrix_;

public:
    // Default constructor
    KOKKOS_FUNCTION
    FMatrixKokkos ();

    // --- 1D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1);

    // Overload operator() to access data as matrix(i), where i = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i) const;

    // --- 2D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as matrix(i, j),
    // where i = [1:N], j = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j) const;

    // --- 3D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as matrix(i, j, k),
    // where i = [1:N], j = [1:N], k = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as matrix(i, j, k, l, m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N],
    // m = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D matrix ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    FMatrixKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as matrix(i, j, k, l, m, n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], 
    // m = [1:N], n = [1:N]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    KOKKOS_FUNCTION
    FMatrixKokkos& operator= (const FMatrixKokkos& temp);

    KOKKOS_FUNCTION
    size_t size() const;

    // Deconstructor
    KOKKOS_FUNCTION
    ~FMatrixKokkos ();

}; // End of FMatrixKokkos

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos() {}

// --- 1D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1) {
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i) const {
    assert(i >= 1 && i <= dim1_);
    return this_matrix_[i - 1];
}

// --- 2D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}

// --- 3D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

// --- 4D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

// --- 5D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_) 
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

// --- 6D matrix ---

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = new T[length_];
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    assert(n >= 1 && n <= dim6_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)  
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)  
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>& FMatrixKokkos<T>::operator= (const FMatrixKokkos& temp)
{
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = new T[length_];
    }
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
size_t FMatrixKokkos<T>::size() const {
    return length_;
}

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::~FMatrixKokkos() {
    delete[] this_matrix_;
}

//==============================================================================
//   C stride array related classes  (last index varies the quickest)
//==============================================================================

/*
 * ViewCArray
 */

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class ViewCArray {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T * this_array_;
    
public:
    
    // Default constructor
    ViewCArray ();
    
    //--- 1D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1);
    
    // Overload  operator() to access data as array(i,j),
    // where i = [0:N-1], j = [0:N-1]
    T& operator()(size_t i) const;
    
    //--- 2D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2);
    
    // Overload operator() to access data as array(i,j),
    //  where i=[0:N-1], j=[0:N-1]
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3);
    
    // Overload operator() to access data as array(i,j,k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4);
        
    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:n-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D array ---
    
    // Overloaded constructor
    ViewCArray (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5);
        
    // Overload operator() to access data as array(i,j,k,l,m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;
    
    //--- 6D array ---
    
    // Overloaded constructor
    ViewCArray (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6);
        
    // Overload operator() to access data as array(i,j,k,l,m,n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;
    
    size_t size() const;
    
}; // end of ViewCArray

// Default constructor
template <typename T>
ViewCArray<T>::ViewCArray() {}
    
//--- 1D array ---
    
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j),
// where i = [0:N-1], j = [0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 1D");  // die if >= dim1
    
    return this_array_[i];
}
   
    
//--- 2D array ---
    
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j),
//  where i=[0:N-1], j=[0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j) const
{
   
    assert(i < dim1_ && "i is out of bounds in c_array 2D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 2D");  // die if >= dim2
    
    return this_array_[j + (i * dim2_)];
}
    
    
//--- 3D array ---
    
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray (T *some_array,
                           size_t some_dim1,
                           size_t some_dim2,
                           size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 3D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 3D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 3D");  // die if >= dim3
    
    return this_array_[k + (j * dim3_) 
                         + (i * dim3_ * dim2_)];
}

    
//--- 4D array ---
    
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i, j, k, l),
// where i = [0:n-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 4D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 4D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 4D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 4D");  // die if >= dim4
    
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_) 
                         + (i * dim4_ * dim3_ * dim2_)];
}

    
//--- 5D array ---
    
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4,
                          size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k,l,m),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l, 
                                    size_t m) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 5D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 5D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 5D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 5D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 5D");  // die if >= dim5
    
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_)
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

//--- 6D array ---
   
// Overloaded constructor
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4,
                          size_t some_dim5,
                          size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k,l,m,n),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 6D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 6D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 6D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 6D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 6D");  // die if >= dim5
    assert(n < dim6_ && "n is out of bounds in c_array 6D");  // die if >= dim6
    
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_) 
                         + (k * dim6_ * dim5_ * dim4_)
                         + (j * dim6_ * dim5_ * dim4_ * dim3_) 
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}
    
template <typename T>
inline size_t ViewCArray<T>::size() const {
    return length_;
} // end of ViewCArray

/*
 * ViewCArrayKokkos
 */

// view a 1D vector as an array(i,...,n), where i=[0:N-1]
template <typename T>
class ViewCArrayKokkos {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_;  // Length of 1D array
    T * this_array_;
    
public:
    // Default constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos ();
    
    //--- 1D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos(T *some_array,
                     size_t some_dim1);
    
    // Overload  operator() to access data as array(i,j),
    // where i = [0:N-1], j = [0:N-1]
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;
    
    //--- 2D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos(T *some_array,
                     size_t some_dim1,
                     size_t some_dim2);
    
    // Overload operator() to access data as array(i,j),
    //  where i=[0:N-1], j=[0:N-1]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos(T *some_array,
                     size_t some_dim1,
                     size_t some_dim2,
                     size_t some_dim3);
    
    // Overload operator() to access data as array(i,j,k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos(T *some_array,
                     size_t some_dim1,
                     size_t some_dim2,
                     size_t some_dim3,
                     size_t some_dim4);
        
    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:n-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos (T *some_array,
                      size_t some_dim1,
                      size_t some_dim2,
                      size_t some_dim3,
                      size_t some_dim4,
                      size_t some_dim5);
        
    // Overload operator() to access data as array(i,j,k,l,m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;
    
    //--- 6D array ---
    
    // Overloaded constructor
    KOKKOS_FUNCTION
    ViewCArrayKokkos (T *some_array,
                      size_t some_dim1,
                      size_t some_dim2,
                      size_t some_dim3,
                      size_t some_dim4,
                      size_t some_dim5,
                      size_t some_dim6);
        
    // Overload operator() to access data as array(i,j,k,l,m,n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;
    
    KOKKOS_FUNCTION
    size_t size() const;
    
}; // end of ViewCArrayKokkos

// Default constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos() {}
    
//--- 1D array ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T *some_array,
                                      size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j),
// where i = [0:N-1], j = [0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 1D");  // die if >= dim1
    
    return this_array_[i];
}
    
//--- 2D array ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T *some_array,
                                      size_t some_dim1,
                                      size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j),
//  where i=[0:N-1], j=[0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, 
                                          size_t j) const
{
   
    assert(i < dim1_ && "i is out of bounds in c_array 2D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 2D");  // die if >= dim2
    
    return this_array_[j + (i * dim2_)];
}
    
//--- 3D array ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos (T *some_array,
                                       size_t some_dim1,
                                       size_t some_dim2,
                                       size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, 
                                          size_t j, 
                                          size_t k) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 3D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 3D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 3D");  // die if >= dim3
    
    return this_array_[k + (j * dim3_) 
                         + (i * dim3_ * dim2_)];
}

    
//--- 4D array ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T *some_array,
                                      size_t some_dim1,
                                      size_t some_dim2,
                                      size_t some_dim3,
                                      size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i, j, k, l),
// where i = [0:n-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, 
                                          size_t j, 
                                          size_t k, 
                                          size_t l) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 4D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 4D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 4D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 4D");  // die if >= dim4
    
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_) 
                         + (i * dim4_ * dim3_ * dim2_)];
}

    
//--- 5D array ---
    
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T *some_array,
                                      size_t some_dim1,
                                      size_t some_dim2,
                                      size_t some_dim3,
                                      size_t some_dim4,
                                      size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k,l,m),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, 
                                          size_t j, 
                                          size_t k, 
                                          size_t l, 
                                          size_t m) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 5D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 5D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 5D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 5D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 5D");  // die if >= dim5
    
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_)
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

//--- 6D array ---
   
// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T *some_array,
                                      size_t some_dim1,
                                      size_t some_dim2,
                                      size_t some_dim3,
                                      size_t some_dim4,
                                      size_t some_dim5,
                                      size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = some_array;
}
    
// Overload operator() to access data as array(i,j,k,l,m,n),
// where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                   size_t m, size_t n) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 6D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 6D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 6D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 6D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 6D");  // die if >= dim5
    assert(n < dim6_ && "n is out of bounds in c_array 6D");  // die if >= dim6
    
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_) 
                         + (k * dim6_ * dim5_ * dim4_)
                         + (j * dim6_ * dim5_ * dim4_ * dim3_) 
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}
    
template <typename T>
KOKKOS_FUNCTION
size_t ViewCArrayKokkos<T>::size() const {
    return length_;
} // end of ViewCArrayKokkos

/*
 * CArray 
 */

template <typename T>
class CArray {
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T* this_array_;

public:
    // Default constructor
    CArray ();

    // --- 1D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1);

    // Overload operator() to access data as array(i), where i = [0:N-1]
    T& operator() (size_t i) const;

    // --- 2D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as array(i, j),
    // where i = [0:N-1], j = [0:N-1]
    T& operator() (size_t i, size_t j) const;

    // --- 3D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as array(i, j, k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as array(i, j, k, l, m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1],
    // m = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as array(i, j, k, l, m, n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], 
    // m = [0:N-1], n = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    CArray& operator= (const CArray& temp); 

    size_t size() const;

    // Deconstructor
    ~CArray ();

}; // End of CArray

template <typename T>
CArray<T>::CArray() {}

// --- 1D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1) {
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = new T[length_];
}

template <typename T>
inline T& CArray<T>::operator() (size_t i) const {
    assert(i < dim1_);
    return this_array_[i];
}

// --- 2D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = new T[length_];
}

template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j) const {
    assert(i < dim1_);
    assert(j < dim2_);
    return this_array_[j + (i * dim2_)];
}

// --- 3D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = new T[length_];
}

template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    return this_array_[k + (j * dim3_) + (i * dim3_ * dim2_)];
}

// --- 4D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_array_ = new T[length_];
}

template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_)  
                         + (i * dim4_ * dim3_ * dim2_)];
}

// --- 5D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_array_ = new T[length_];
}

template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_) 
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

// --- 6D array ---

template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = new T[length_];

}

template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    assert(n < dim6_);
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_)  
                         + (k * dim6_ * dim5_ * dim4_) 
                         + (j * dim6_ * dim5_ * dim4_ * dim3_)  
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];

}

template <typename T>
inline CArray<T>& CArray<T>::operator= (const CArray& temp)
{
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_array_ = new T[length_];
    }
    return *this;
}

template <typename T>
inline size_t CArray<T>::size() const {
    return length_;
}

template <typename T>
CArray<T>::~CArray() {
    delete[] this_array_;
}

/*
 * CArrayKokkos 
 */

template <typename T>
class CArrayKokkos {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    TArray1D this_array_; 

public:
    // Default constructor
    CArrayKokkos ();

    // --- 1D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1);

    // Overload operator() to access data as array(i), where i = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i) const;

    // --- 2D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as array(i, j),
    // where i = [0:N-1], j = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j) const;

    // --- 3D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as array(i, j, k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as array(i, j, k, l, m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1],
    // m = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D array ---
    
    // Overloaded constructor
    CArrayKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as array(i, j, k, l, m, n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], 
    // m = [0:N-1], n = [0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    CArrayKokkos& operator= (const CArrayKokkos& temp);

    // Method that returns size
    KOKKOS_FUNCTION
    size_t size() const;

    // Deconstructor
    KOKKOS_FUNCTION
    ~CArrayKokkos ();
}; // End of CArrayKokkos

template <typename T>
CArrayKokkos<T>::CArrayKokkos() {}

// --- 1D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i) const {
    assert(i < dim1_);
    return this_array_[i];
}

// --- 2D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i, size_t j) const {
    assert(i < dim1_);
    assert(j < dim2_);
    return this_array_[j + (i * dim2_)];
}

// --- 3D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    return this_array_[k + (j * dim3_) + (i * dim3_ * dim2_)];
}

// --- 4D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_)  
                         + (i * dim4_ * dim3_ * dim2_)];
}

// --- 5D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_) 
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

// --- 6D array ---

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    assert(n < dim6_);
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_)  
                         + (k * dim6_ * dim5_ * dim4_) 
                         + (j * dim6_ * dim5_ * dim4_ * dim3_)  
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}

template <typename T>
CArrayKokkos<T>& CArrayKokkos<T>::operator= (const CArrayKokkos& temp)
{

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_array_ = TArray1D("this_array_", length_);
    }
    
    return *this;
}

//return size
template <typename T>
KOKKOS_FUNCTION
size_t CArrayKokkos<T>::size() const {
	return length_;
}

template <typename T>
KOKKOS_FUNCTION
CArrayKokkos<T>::~CArrayKokkos() {
}

//~~~~~~~~~~~~~~~Begin CMatrix Class~~~~~~~~~~~~~~~~~~~~~~~`

// create a 1D vector that is accessed as an array(i,...,n), where i=[1:N-]
template <typename T>
class CMatrix {
        
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_,length_;
    T * this_matrix;
            
public:
        
       // default constructor
       CMatrix();
       CMatrix(size_t some_dim1);
       CMatrix(size_t some_dim1, size_t some_dim2);
       CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3);
       CMatrix(size_t some_dim1,
           size_t some_dim2,
           size_t some_dim3,
           size_t some_dim4);
       CMatrix(size_t some_dim1,
           size_t some_dim2,
           size_t some_dim3,
           size_t some_dim4,
           size_t some_dim5);
       CMatrix (size_t some_dim1,
            size_t some_dim2,
            size_t some_dim3,
            size_t some_dim4,
            size_t some_dim5,
            size_t some_dim6);
           
    //overload operators to access data
       T& operator()(size_t i);
       T& operator()(size_t i, size_t j);
       T& operator()(size_t i, size_t j, size_t k);
       T& operator()(size_t i, size_t j, size_t k, size_t l);
       T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
       T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

       //overload = operator
	CMatrix& operator= (const CMatrix &temp);
            
       // deconstructor
       ~CMatrix( );
        
}; // end of c_matrix_t


//implement functions
//First, constructors from 1D to 6D

//no dimenstion
template <typename T>
CMatrix<T>::CMatrix() {}

//1D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix = new T[length_];
}

//2D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_matrix = new T[length_];
}

//3D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_matrix = new T[length_];
}

//4D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_matrix= new T[length_];
}   

//5D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_matrix = new T[length_];
}

//6D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_matrix = new T[length_];
}

//~~~~operators() for 1D to 6D
//NOTE: the indices start at 1 for matrices
//i [1,N], j [1,N],...
template <typename T>
T& CMatrix<T>::operator()(size_t i)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 1D!");
    return this_matrix[i-1];
}

//2D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 2D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 2D!");
    return this_matrix[(j-1) + (i-1)*dim2_];
}

//3D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 3D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 3D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 3D!");
    return this_matrix[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 4D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 4D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 4D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 4D!");
    return this_matrix[ (l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l,size_t m)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 5D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 5D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 5D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 5D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrix 5D!");
    return this_matrix[(m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 6D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 6D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 6D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 6D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrix 6D!");
    assert( n < dim6_+1 && "n is out of bounds in CMatrix 6D!");
    return this_matrix[ (n-1) + (m-1)*dim6_ + (l-1)*dim6_*dim5_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim6_*dim5_*dim4_*dim3_*dim2_];
}

//overload = operator
//THIS = CMatrix<> temp
template <typename T>
CMatrix<T> &CMatrix<T>::operator= (const CMatrix &temp) {
	if(this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix = new T[length_];
	}
}

// Destructor
template <typename T>
CMatrix<T>::~CMatrix(){
    delete[] this_matrix;
}

//~~~~~KOKKOS CMATRIX~~~~~
template <typename T>
class CMatrixKokkos {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
private:
    size_t dim1_, dim2_, dim3_, dim4_,dim5_, dim6_, length_;
    TArray1D this_matrix_; 

public:
    //default constructor
    CMatrixKokkos();

    //overload constructor from 1D to 6D
    CMatrixKokkos(size_t some_dim1);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3);    

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Over load operator () to access data 
    //Note for matrix class, indices start at 1
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    // Overload = operator
    CMatrixKokkos& operator=(const CMatrixKokkos &temp);

    // Destructor
    KOKKOS_FUNCTION
    ~CMatrixKokkos();

}; // End of CMatrisKokkos declarations

//CMatrixKokkos Definitions

//Starting with constructors

// No dimension
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos() {}

//1D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1) { 

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//2D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2) { 

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//3D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//4D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//5D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//6D
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

//operators to access data
//From 1D to 6D
//Recall indices start at 1

//1D
template<typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i) const {
    assert( i < dim1_+1 && "i is out of bounds in CMatrixKokkos 1D!");
    return this_matrix_[i-1];
}

//2D
template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i < dim1_+1 && "i is out of bounds in CMatrixKokkos 2D!");
    assert(j < dim2_+1 && "j is out of bounds in CMatrixKokkos 2D!");
    return this_matrix_[ (j-1) + (i-1)*dim2_];
}

//3D
template<typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert( i < dim1_+1 && "i is out of bounds in CMatrixKokkos 3D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrixKokkos 3D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrixKokkos 3D!");
    return this_matrix_[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert( i < dim1_+1 && "i is out of bounds in CMatrixKokkos 4D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrixKokkos 4D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrixKokkos 4D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrixKokkos 4D!");
    return this_matrix_[ (l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const {
    assert( i < dim1_+1 && "i is out of bounds in CMatrixKokkos 5D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrixKokkos 5D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrixKokkos 5D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrixKokkos 5D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrixKokkos 5D!");
    return this_matrix_[ (m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const {
    assert( i < dim1_+1 && "i is out of bounds in CMatrixKokkos 6D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrixKokkos 6D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrixKokkos 6D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrixKokkos 6D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrixKokkos 6D!");
    assert( n < dim6_+1 && "n is out of bounds in CMatrixKokkos 6D!");
    return this_matrix_[(n-1) + (m-1)*dim6_ + (l-1)*dim6_*dim5_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim6_*dim5_*dim4_*dim3_*dim2_];
}

// Overload = operator
// for object assignment THIS = CMatrixKokkos <> temp
template <typename T>
CMatrixKokkos<T> & CMatrixKokkos<T>::operator= (const CMatrixKokkos &temp) {
	if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = TArray1D("this_matrix_", length_);
    }
	
    return *this;
}
    
// Deconstructor
template <typename T>
KOKKOS_FUNCTION
CMatrixKokkos<T>::~CMatrixKokkos() {
}

//==============================================================================
//   Ragged right array related classes 
//==============================================================================

template <typename T>
class RaggedRightArray {
private:
    size_t *start_index_;
    T * array_;
    
    size_t dim1_, length_;
    
public:
    // Default constructor
    RaggedRightArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    RaggedRightArray (CArray<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    RaggedRightArray (ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    RaggedRightArray (size_t *strides_array, size_t some_dim1);
    
    // A method to return the stride size
    size_t stride(size_t i) const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;

    RaggedRightArray& operator= (const RaggedRightArray &temp);

    // Destructor
    ~RaggedRightArray ( );
}; // End of RaggedRightArray

// Overloaded constructor
template <typename T>
RaggedRightArray<T>::RaggedRightArray (CArray<size_t> &strides_array){
    // The length of the stride array is some_dim1;
    //dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// Overloaded constructor
template <typename T>
RaggedRightArray<T>::RaggedRightArray (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    //dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// Overloaded constructor
template <typename T>
RaggedRightArray<T>::RaggedRightArray (size_t *strides_array, size_t dim1){
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// A method to return the stride size
template <typename T>
inline size_t RaggedRightArray<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_ + 1) && "i is greater than dim1_ in RaggedRightArray");

    return start_index_[(i + 1)] - start_index_[i];
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& RaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // get the 1D array index
    size_t start = start_index_[i];
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArray");  // die if >= stride
    
    return array_[j + start];
} // End operator()

template <typename T>
RaggedRightArray<T> & RaggedRightArray<T>::operator= (const RaggedRightArray &temp) {

    if( this != &temp) {
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        start_index_ = new size_t[dim1_ + 1];
        for (int j = 0; j < dim1_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
RaggedRightArray<T>::~RaggedRightArray () {
    delete[] array_;
    delete[] start_index_;
}

template <typename T>
class RaggedRightArrayKokkos {
private:
    size_t *start_index_;
    T * array_;
    
    size_t dim1_, length_;
    
public:
    // Default constructor
    KOKKOS_FUNCTION
    RaggedRightArrayKokkos ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    KOKKOS_FUNCTION
    RaggedRightArrayKokkos (CArray<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    KOKKOS_FUNCTION
    RaggedRightArrayKokkos (ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    KOKKOS_FUNCTION
    RaggedRightArrayKokkos (size_t *strides_array, size_t some_dim1);
    
    // A method to return the stride size
    KOKKOS_FUNCTION
    size_t stride(size_t i) const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    RaggedRightArrayKokkos& operator= (const RaggedRightArrayKokkos &temp);

    // Destructor
    KOKKOS_FUNCTION
    ~RaggedRightArrayKokkos ( );
}; // End of RaggedRightArray

// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos (CArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos (size_t *strides_array, size_t dim1) {
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[(i + 1)] = count;
    } // end for i
    
    array_ = new T[count];
} // End constructor

// A method to return the stride size
template <typename T>
KOKKOS_FUNCTION
size_t RaggedRightArrayKokkos<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_ + 1) && "i is greater than dim1_ in RaggedRightArray");

    return start_index_[(i + 1)] - start_index_[i];
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
KOKKOS_FUNCTION
T& RaggedRightArrayKokkos<T>::operator()(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArray");  // die if >= stride
    
    return array_[j + start];
} // End operator()

template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T> & RaggedRightArrayKokkos<T>::operator= (const RaggedRightArrayKokkos &temp) {

    if( this != &temp) {
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        start_index_ = new size_t[dim1_ + 1];
        for (int j = 0; j < dim1_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T>::~RaggedRightArrayKokkos () {
    delete[] array_;
    delete[] start_index_;
}

//~~~~~~~~~~~begin RaggedDownArray class declarations~~~~~~~~~~~~~~~~~~~`
template <typename T>
class RaggedDownArray { 
private:
    size_t *start_index_;
	T * array_;

	size_t dim2_;
    size_t length_;

public:
    //default constructor
    RaggedDownArray() ;

    //~~~~2D`~~~~
	//overload constructor with CArray
	RaggedDownArray(CArray<size_t> &strides_array);

	//overload with ViewCArray
	RaggedDownArray(ViewCArray <size_t> &strides_array);

	//overload with traditional array
	RaggedDownArray(size_t *strides_array, size_t dome_dim1);

	//method to return stride size
	size_t stride(size_t j);

	//overload () operator to access data as array (i,j)
	T& operator()(size_t i, size_t j);

    // method to return total size
    size_t size();

    RaggedDownArray& operator= (const RaggedDownArray &temp);

    //destructor
    ~RaggedDownArray();

}; //~~~~~end of RaggedDownArray class declarations~~~~~~~~	

//~~~~~~~~~~~~~~~begin RaggedDownArray class definitions~~~~~~~~~~~

// Overload constructor with CArray for the strides array
template <typename T>
RaggedDownArray<T>::RaggedDownArray() {}


template <typename T>
RaggedDownArray<T>::RaggedDownArray( CArray <size_t> &strides_array) {
    // Length of stride array
    //dim2_ = strides_array.size();

    // Create and initialize startding indices
    start_index_ = new size_t[dim2_+1]; //theres a plus 1, because 
    start_index_[0] = 0; //1D array starts at 0

		
	//length of strides
	dim2_ = strides_array.size();

    // Loop to find total length of 1D array
    size_t count = 0;
    for(size_t j = 0; j < dim2_ ; j++) { 
        count += strides_array(j);
        start_index_[j+1] = count;
    } 
    length_ = count;

    array_ = new T[count];

} // End constructor 

// Overload constructor with ViewCArray
template <typename T>
RaggedDownArray<T>::RaggedDownArray( ViewCArray <size_t> &strides_array) {
    // Length of strides
    //dim2_ = strides_array.size();

    //create array for holding start indices
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    size_t count = 0;
    // Loop over to get total length of 1D array
    for(size_t j = 0; j < dim2_ ;j++ ) {
        count += strides_array(j);
        start_index_[j+1] = count;
    }
    length_ = count;	
    array_ = new T[length_];

} // End constructor 

// Overload constructor with regualar array
template <typename T>
RaggedDownArray<T>::RaggedDownArray( size_t *strides_array, size_t dim2){
    // Length of stride array
    dim2_ = dim2;

    // Create and initialize starting index of entries
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    // Loop over to find length of 1D array
    // Represent ragged down array and set 1D index
    size_t count = 0;
    for(size_t j = 0; j < dim2_; j++) {
        count += strides_array[j];
        start_index_[j+1] = count;
	}

    length_ = count;	
    array_ = new T[length_];

} //end construnctor

// Check the stride size
template <typename T>
size_t RaggedDownArray<T>::stride(size_t j) {
    return start_index_[j+1] - start_index_[j];
}

template <typename T>
size_t RaggedDownArray<T>::size() {
    return length_;
}

// overload operator () to access data as an array(i,j)
// Note: i = 0:stride(j), j = 0:N-1
template <typename T>
T& RaggedDownArray<T>::operator()(size_t i, size_t j) {
    // Where is the array starting?
    // look at start index
    size_t start = start_index_[j]; 

    // Make sure we are within array bounds
    assert(i < stride(j) && "i is out of bounds in RaggedDownArray");
    assert(j < dim2_ && "j is out of dim2_ bounds in RaggedDownArray");
		
    return array_[i + start];

} // End () operator

template <typename T>
RaggedDownArray<T> & RaggedDownArray<T>::operator= (const RaggedDownArray &temp) {

    if( this != &temp) {
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        start_index_ = new size_t[dim2_ + 1];
        for (int j = 0; j < dim2_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
RaggedDownArray<T>::~RaggedDownArray() {
    delete[] array_;
    delete[] start_index_;

} // End destructor

//~~~~~~~~~~~~~~~~~~~end RaggedDownArray class definitions~~~~~~~


//~~~~~~~~~~~~~~~~~~~begin Complressed Sparse Column Kokkos Storage~~~~~
// csc storage
// array.value(i,j)
// array.row_index(i,j)

template <typename T>
class sparseColumnArrayKokkos {
private:
    size_t *start_index_;
    size_t *row_index_;
    T * array_;

    size_t dim2_;

public:
    // Default constructor
    KOKKOS_FUNCTION
    sparseColumnArrayKokkos(); 

    //---2D array---
		
    // Overload constructor for CArray
    KOKKOS_FUNCTION
    sparseColumnArrayKokkos( CArrayKokkos <size_t> &strides_array);

    // Overload constructor for a ViewCArray
    KOKKOS_FUNCTION
    sparseColumnArrayKokkos(ViewCArrayKokkos<size_t> &strides_array);

    // Overload constructor for a regular cpp array
    KOKKOS_FUNCTION
    sparseColumnArrayKokkos(size_t *strides_array, size_t some_dim1);

    // Method to return the strides size
    KOKKOS_FUNCTION
    size_t stride(size_t j) const;

    // Method to return row index; array.row_index(i,j)
    KOKKOS_FUNCTION
    size_t& row_index(size_t i, size_t j) const;

    // Method to access data as array.value(i,j),
    // where i = 0:N-1, j = stride(i)
    KOKKOS_FUNCTION
    T& value(size_t i, size_t j) const;

    // Destructor
    KOKKOS_FUNCTION
    ~sparseColumnArrayKokkos ();
}; //~~~~~~~~~~end of sparseColumnArray class declaration~~~~~~~~~

//~~~~~~~~~~~~~~~~begin sparseColumnArray class defitions

// Overload constructor with CArray
template <typename T>
KOKKOS_FUNCTION
sparseColumnArrayKokkos<T>::sparseColumnArrayKokkos(CArrayKokkos <size_t> &strides_array) {
    // Length of stride array
    // dim2_ = strides_array.size();

    // Create and initialize startind index 
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    // Loop over to get total length of 1D array
    size_t count = 0;
    for(size_t j = 0; j < dim2_ ; j++) {
        count += strides_array(j);
        start_index_[j+1] = count;
    }

    array_ = new T[count];
    row_index_ = new T[count];
} // End constructor with CArray

// Overload constructor with view array
template <typename T>
KOKKOS_FUNCTION
sparseColumnArrayKokkos<T>::sparseColumnArrayKokkos(ViewCArrayKokkos <size_t> &strides_array) {
    // Length of stride array
    //dim2_ = strides_array.size();

    // Create and initialize start indices
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    // Loop over total length of 1D array
    size_t count = 0;
    for( size_t j = 0; j < dim2_ ; j++) {
        count += strides_array(j);
        start_index_[j+1] = count;
    }

    array_ = new T[count];
    row_index_ = new T[count];
} // End constructor with ViewCArray

//overload constructor with a regular cpp array
template <typename T>
KOKKOS_FUNCTION
sparseColumnArrayKokkos<T>::sparseColumnArrayKokkos(size_t *strides_array, size_t dim2) {
    // Length of  stride array
    dim2_ = dim2;

    // Create and initialize staring index array
    start_index_ = new T[dim2_+1];
    start_index_[0] = 0;

    // Loop over total length of 1D
    size_t count = 0;
    for(size_t j = 0; j < dim2_; j++){
        count += strides_array[j];
        start_index_[j+1] = count;
    }

    array_ = new T[count];
    row_index_ = new T[count];
} // End constructor with regular array

// Method to return stride size
template <typename T>
KOKKOS_FUNCTION
size_t sparseColumnArrayKokkos<T>::stride(size_t j) const {
    return start_index_[j+1] - start_index_[j];
}

// Access data as array.row_index(i,j)
// i = stride(j), j = 0:N-1
template <typename T>
KOKKOS_FUNCTION
size_t& sparseColumnArrayKokkos<T>::row_index(size_t i, size_t j) const {
    // Get 1D array index
    size_t start = start_index_[j];

    // Assert
    assert(i < stride(j) && "i is out stride bounds in sparseColumnArrayKokkos");
    assert(j < dim2_ && "is out out of dim2 bounds in sparseColumnArrayKokkos");

    return row_index_[i+start];
} // End row index method

// Access data as array.value(i,j)
// where i = 0:stride(j), j = 0:N-1
template <typename T>
KOKKOS_FUNCTION
T& sparseColumnArrayKokkos<T>::value(size_t i, size_t j) const {
    // Get 1D array index
    size_t start = start_index_[j];

    // Assert indices are within bounds
	assert(i < stride(j) && "i is out of stride bounds in sparseColumnArrayKokkos");
    assert(j < dim2_ && "j is out of dim2 bounds in sparseColumnArrayKokkos");

    return array_[i + start];
} // End method

// Destructor
template <typename T>
KOKKOS_FUNCTION
sparseColumnArrayKokkos<T>::~sparseColumnArrayKokkos () {
    delete[] array_;
    delete[] start_index_;
    delete[] row_index_;
}

//~~~~~~~~~~~~~~~~end of sparseColumnArrayKokkos class definitions ~~~~~~~~~~~~

// A class to access the row indicies i of an Array(i,j), like CArray,
// that is sparsely populated but has dense sections along i
//     get i_global using (i_local, j_global)
/*
template <typename T>
class SparseRowAccessors {
private:
    size_t *stride_;
    T * dim1_index_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    SparseRowAccessors ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overloaded constructor
    SparseRowAccessors (size_t dim1, size_t dim2);
    
    // A method to return or set the stride size
    size_t& stride(size_t j) const;
    
    // A method to increase the stride size
    void push_back(size_t j) const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;
    
    // Destructor
    ~SparseRowAccessors ();
}; 

// Overloaded constructor
SparseRowAccessors::SparseRowAccessors (size_t dim1, size_t dim2) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1 * dim2;
    
    // Create memory on the heap for the values
    dim1_index_ = new T[length_];
    
    // Create memory for the stride size in each column
    stride_ = new size_t[dim2];
    
    // Initialize the stride to a dense array so stride_ = 0
    for (int j=0; j<dim2_; j++){
        stride_[j] = 0;
    }
    
    // Start index is always = i + j*dim1
} 

// A method to set the stride size for column j
size_t& SparseRowAccessors::stride(size_t j) const {
    return stride_[j];
}

// A method to increase the stride size for row i
void SparseRowAccessors::push_back(size_t j) const {
    stride_[j]++;
}

// Overload operator() to access data as array(i,j),
// where i=[0:N-1], j=[0:stride(i)]
inline T& SparseRowAccessors::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in SparseRowAccessors");  // die if >= dim1
    assert(j < dim2_ && "j is out of dim2 bounds in SparseRowAccessors");  // die if >= dim2
    assert(i < stride_[j] && "j is out of stride bounds in sparaseRowAccessors");  // die if >= stride
    
    return dim1_index_[i + j*dim1_];
} 

// Destructor
SparseRowAccessors::~SparseRowAccessors () {
    delete[] dim1_index_;
    delete[] stride_;
}
*/

// The DynamicRaggedRightArray is designed to include a buffer
template <typename T>
class DynamicRaggedRightArray {
private:
    size_t *stride_;
    T * array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedRightArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor
    DynamicRaggedRightArray (size_t dim1, size_t dim2);
    
    // A method to return or set the stride size
    size_t& stride(size_t i) const;
    
    // A method to increase the stride size
    void push_back(size_t i) const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;
    
    // Destructor
    ~DynamicRaggedRightArray ();
}; 

// Overloaded constructor
template <typename T>
DynamicRaggedRightArray<T>::DynamicRaggedRightArray (size_t dim1, size_t dim2) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;
    
    // Create memory on the heap for the values
    array_ = new T[dim1*dim2];
    
    // Create memory for the stride size in each row
    stride_ = new size_t[dim1];
    
    // Initialize the stride
    for (int i=0; i<dim1_; i++){
        stride_[i] = 0;
    }
    
    // Start index is always = j + i*dim2
} 

// A method to set the stride size for row i
template <typename T>
size_t& DynamicRaggedRightArray<T>::stride(size_t i) const {
    return stride_[i];
}


// A method to increase the stride size for row i
template <typename T>
void DynamicRaggedRightArray<T>::push_back(size_t i) const {
    stride_[i]++;
}

// Overload operator() to access data as array(i,j),
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& DynamicRaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in DynamicRaggedRight");  // die if >= dim1
    assert(j < dim2_ && "j is out of dim2 bounds in DynamicRaggedRight");  // die if >= dim2
    assert(j < stride_[i] && "j is out of stride bounds in DynamicRaggedRight");  // die if >= stride
    
    return array_[j + i*dim2_];
}

// Destructor
template <typename T>
DynamicRaggedRightArray<T>::~DynamicRaggedRightArray() {
    delete[] array_;
    delete[] stride_;
}

// The DynamicRaggedRightArray is designed to include a buffer
/*
template <typename T>
class DynamicRaggedRightArrayKokkos {
private:
    size_t *stride_;
    T * array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedRightArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // overload constructor
    DynamicRaggedRightArray (size_t dim1, size_t dim2);
    
    // A method to return or set the stride size
    size_t& stride(size_t i) const;
    
    // A method to increase the stride size
    void push_back(size_t i) const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;
    
    // Destructor
    ~DynamicRaggedRightArray ();
}; 

// Overloaded constructor
template <typename T>
DynamicRaggedRightArray<T>::DynamicRaggedRightArray (size_t dim1, size_t dim2) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;
    
    // Create memory on the heap for the values
    array_ = new T[dim1*dim2];
    
    // Create memory for the stride size in each row
    stride_ = new size_t[dim1];
    
    // Initialize the stride
    for (int i=0; i<dim1_; i++){
        stride_[i] = 0;
    }
    
    // Start index is always = j + i*dim2
} 

// A method to set the stride size for row i
template <typename T>
size_t& DynamicRaggedRightArray<T>::stride(size_t i) const {
    return stride_[i];
}


// A method to increase the stride size for row i
template <typename T>
void DynamicRaggedRightArray<T>::push_back(size_t i) const {
    stride_[i]++;
}

// Overload operator() to access data as array(i,j),
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& DynamicRaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in DynamicRaggedRight");  // die if >= dim1
    assert(j < dim2_ && "j is out of dim2 bounds in DynamicRaggedRight");  // die if >= dim2
    assert(j < stride_[i] && "j is out of stride bounds in DynamicRaggedRight");  // die if >= stride
    
    return array_[j + i*dim2_];
}

// Destructor
template <typename T>
DynamicRaggedRightArray<T>::~DynamicRaggedRightArray() {
    delete[] array_;
    delete[] stride_;
}
*/


template <typename T>
class SparseRowArray {
private:
    size_t *start_index_;
    size_t *column_index_;
    
    T * array_;
    
    size_t dim1_;
    
public:
    // Default constructor
    SparseRowArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    SparseRowArray (CArray<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    SparseRowArray (ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    SparseRowArray (size_t *strides_array, size_t some_dim1);
    
    // A method to return the stride size
    size_t stride(size_t i) const;
    
    // A method to return the column index as array.column_index(i,j)
    size_t column_index(size_t i, size_t j) const;
    
    // A method to access data as array.value(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& value(size_t i, size_t j) const;
    
    // Destructor
    ~SparseRowArray ();
}; 


// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (CArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 


// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 

// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (size_t *strides_array, size_t dim1) {
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 


// A method to return the stride size
template <typename T>
size_t SparseRowArray<T>::stride(size_t i) const {
    return start_index_[i+1] - start_index_[i];
}

// A method to return the column index
template <typename T>
size_t SparseRowArray<T>::column_index(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in SparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in SparseRowArray");  // die if >= stride
    
    return column_index_[j + start];
}

// Access data as array.value(i,j), 
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& SparseRowArray<T>::value(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in sparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in sparseRowArray");  // die if >= stride
    
    return array_[j + start];
} 

// Destructor
template <typename T>
SparseRowArray<T>::~SparseRowArray() {
    delete[] array_;
    delete[] start_index_;
    delete[] column_index_;
}

template <typename T>
class SparseRowArrayKokkos {
private:
    size_t *start_index_;
    size_t *column_index_;
    
    T * array_;
    
    size_t dim1_;
    
public:
    // Default constructor
    KOKKOS_FUNCTION
    SparseRowArrayKokkos ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    KOKKOS_FUNCTION
    SparseRowArrayKokkos (CArray<size_t> &strides_array);
    
    // overload constructor for a ViewCArray
    KOKKOS_FUNCTION
    SparseRowArrayKokkos (ViewCArray<size_t> &strides_array);
    
    // overloaded constructor for a traditional array
    KOKKOS_FUNCTION
    SparseRowArrayKokkos (size_t *strides_array, size_t some_dim1);
    
    
    // a method to return the stride size
    KOKKOS_FUNCTION
    size_t stride(size_t i) const;
    
    
    // a method to return the column index as array.column_index(i,j)
    KOKKOS_FUNCTION
    size_t column_index(size_t i, size_t j) const;
    
    
    // a method to access data as array.value(i,j)
    //  where i=[0:N-1], j=[stride(i)]
    KOKKOS_FUNCTION
    T& value(size_t i, size_t j) const;
    
    // destructor
    KOKKOS_FUNCTION
    ~SparseRowArrayKokkos ();
}; // End of SparseRowArrayKokkos

// overloaded constructor
template <typename T>
KOKKOS_FUNCTION
SparseRowArrayKokkos<T>::SparseRowArrayKokkos (CArray<size_t> &strides_array) {
    
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 

// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
SparseRowArrayKokkos<T>::SparseRowArrayKokkos (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 

// Overloaded constructor
template <typename T>
KOKKOS_FUNCTION
SparseRowArrayKokkos<T>::SparseRowArrayKokkos (size_t *strides_array, 
                                               size_t dim1) {
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[i+1] = count;
    } // end for i
    
    array_ = new T[count];
    column_index_ = new size_t[count];
} 

// A method to return the stride size
template <typename T>
KOKKOS_FUNCTION
size_t SparseRowArrayKokkos<T>::stride(size_t i) const {
    return start_index_[i+1] - start_index_[i];
}

// A method to return the column index
template <typename T>
KOKKOS_FUNCTION
size_t SparseRowArrayKokkos<T>::column_index(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in sparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in sparseRowArray");  // die if >= stride
    
    return column_index_[(j + start)];
}

// Access data as array.value(i,j),
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
KOKKOS_FUNCTION
T& SparseRowArrayKokkos<T>::value(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in sparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in sparseRowArray");  // die if >= stride
    
    return array_[(j + start)];
} 


// Destructor
template <typename T>
KOKKOS_FUNCTION
SparseRowArrayKokkos<T>::~SparseRowArrayKokkos() {
    delete[] array_;
    delete[] start_index_;
    delete[] column_index_;
}

//~~~~~~~~~~~begin RaggedDownArrayKokkos class declarations~~~~~~~~~~~~~~~~~~~`
template <typename T>
class RaggedDownArrayKokkos { 
private:
    size_t *start_index_;
	T * array_;

	size_t dim2_;
    size_t length_;

public:
    //default constructor
    KOKKOS_FUNCTION
    RaggedDownArrayKokkos() ;

    //~~~~2D`~~~~
	//overload constructor with CArrayKokkos
	KOKKOS_FUNCTION
	RaggedDownArrayKokkos(CArrayKokkos<size_t> &strides_array);

	//overload with ViewCArrayKokkos
	KOKKOS_FUNCTION
	RaggedDownArrayKokkos(ViewCArrayKokkos <size_t> &strides_array);

	//overload with traditional array
	KOKKOS_FUNCTION
	RaggedDownArrayKokkos(size_t *strides_array, size_t dome_dim1);

	//method to return stride size
	KOKKOS_FUNCTION
	size_t stride(size_t j);

	//overload () operator to access data as array (i,j)
	KOKKOS_FUNCTION
	T& operator()(size_t i, size_t j);

    // method to return total size
    KOKKOS_FUNCTION
    size_t size();

    KOKKOS_FUNCTION
    RaggedDownArrayKokkos& operator= (const RaggedDownArrayKokkos &temp);

    //destructor
    KOKKOS_FUNCTION
    ~RaggedDownArrayKokkos();

}; //~~~~~end of RaggedDownArrayKokkos class declarations~~~~~~~~	

//~~~~~~~~~~~~~~~begin RaggedDownArrayKokkos class definitions~~~~~~~~~~~

// Overload constructor with CArrayKokkos for the strides array
template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos() {}


template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos( CArrayKokkos <size_t> &strides_array) {
    // Length of stride array
    //dim2_ = strides_array.size();

    // Create and initialize startding indices
    start_index_ = new size_t[dim2_+1]; //theres a plus 1, because 
    start_index_[0] = 0; //1D array starts at 0

		
	//length of strides
	dim2_ = strides_array.size();

    // Loop to find total length of 1D array
    size_t count = 0;
    for(size_t j = 0; j < dim2_ ; j++) { 
        count += strides_array(j);
        start_index_[j+1] = count;
    } 
    length_ = count;

    array_ = new T[count];

} // End constructor 

// Overload constructor with ViewCArrayKokkos
template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos( ViewCArrayKokkos <size_t> &strides_array) {
    // Length of strides
    //dim2_ = strides_array.size();

    //create array for holding start indices
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    size_t count = 0;
    // Loop over to get total length of 1D array
    for(size_t j = 0; j < dim2_ ;j++ ) {
        count += strides_array(j);
        start_index_[j+1] = count;
    }
    length_ = count;	
    array_ = new T[length_];

} // End constructor 

// Overload constructor with regualar array
template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos( size_t *strides_array, size_t dim2){
    // Length of stride array
    dim2_ = dim2;

    // Create and initialize starting index of entries
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    // Loop over to find length of 1D array
    // Represent ragged down array and set 1D index
    size_t count = 0;
    for(size_t j = 0; j < dim2_; j++) {
        count += strides_array[j];
        start_index_[j+1] = count;
	}

    length_ = count;	
    array_ = new T[length_];

/*
	template <typename T>
	KOKKOS_FUNCTION
	RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos(size_t *strides_array, size_t dim2, Kokkos::TeamPolicy<Kokkos::Cuda>::member_type teamMember) {
		//length of stride array
		dim2_ = dim2;

		//create and initialize starting index of entries
		start_index_ = new size_t[dim2_+1];
		start_index_[0] = 0;
        for (int i = 0; i < dim2_; i++) {
            printf("hello %d %d\n", i, strides_array[i]);
        }

		//loop over to find length of 1D array
		//represent ragged down array and set 1D index
		//size_t count = 0;
        Kokkos::parallel_scan (Kokkos::ThreadVectorRange(teamMember, dim2_), [=] (const int j, size_t& count, const bool final) {
            const size_t idx_temp = strides_array[j];
            printf("%d %d\n", j, idx_temp);
	        if (final) {
                start_index_[j+1] = count; 
                //printf("%d %d\n", j, start_index_[j+1]);
            }
	        count += idx_temp;
        });
	    length_ = start_index_[dim2_];	
		array_ = new T[length_];

	} //end construnctor
*/

} //end construnctor

// Check the stride size
template <typename T>
KOKKOS_FUNCTION
size_t RaggedDownArrayKokkos<T>::stride(size_t j) {
    return start_index_[j+1] - start_index_[j];
}

template <typename T>
KOKKOS_FUNCTION
size_t RaggedDownArrayKokkos<T>::size() {
    return length_;
}

// overload operator () to access data as an array(i,j)
// Note: i = 0:stride(j), j = 0:N-1
template <typename T>
KOKKOS_FUNCTION
T& RaggedDownArrayKokkos<T>::operator()(size_t i, size_t j) {
    // Where is the array starting?
    // look at start index
    size_t start = start_index_[j]; 

    // Make sure we are within array bounds
    assert(i < stride(j) && "i is out of bounds in RaggedDownArrayKokkos");
    assert(j < dim2_ && "j is out of dim2_ bounds in RaggedDownArrayKokkos");
		
    return array_[i + start];

} // End () operator

template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T> & RaggedDownArrayKokkos<T>::operator= (const RaggedDownArrayKokkos &temp) {

    if( this != &temp) {
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        start_index_ = new size_t[dim2_ + 1];
        for (int j = 0; j < dim2_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::~RaggedDownArrayKokkos() {
    delete[] array_;
    delete[] start_index_;

} // End destructor


template <typename T>
class k_test_t {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
private:
    size_t dim1, dim2, dim3, dim4, dim5, dim6, length;
    //T * this_array;
    TArray1D this_array; 
public:
    // default constructor
    //KOKKOS_FUNCTION
    k_test_t ();
    
    //--- 1D array ---
   
    // overloaded constructor
    //KOKKOS_FUNCTION
    k_test_t (size_t some_dim1);

    //KOKKOS_FUNCTION
    k_test_t (size_t some_dim1, size_t some_dim2);

    // overload opperator() to access data as array(i) where i=[0:N-1]
    KOKKOS_FUNCTION
    T& operator() (size_t i) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j) const;

    // overload copy assignment operator
    //KOKKOS_FUNCTION
    k_test_t& operator= (const k_test_t& temp);

    // destructor
    KOKKOS_FUNCTION
    ~k_test_t ();
};


template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t() {}

template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t(size_t some_dim1) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1 = some_dim1;
    length = dim1;
    //this_array = new T[length];
    //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
    this_array = TArray1D("this_array", length);
}

template <typename T>
KOKKOS_FUNCTION
T& k_test_t<T>::operator() (size_t i) const {
    return this_array(i);
}

template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>::k_test_t(size_t some_dim1, size_t some_dim2) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1 = some_dim1;
    dim2 = some_dim2;
    length = dim1 * dim2;
    //this_array = new T[length];
    //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
    this_array = TArray1D("this_array", length);
}

template <typename T>
KOKKOS_FUNCTION
T& k_test_t<T>::operator() (size_t i, size_t j) const {
    return this_array(j + i*dim2);
}

// for the object assignment THIS = k_array_t <> TEMP(n,m,etc.) 
template <typename T>
//KOKKOS_FUNCTION
k_test_t<T>& k_test_t<T>::operator= (const k_test_t& temp) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // do nothing if the assignment is of form x = x
    if (this != &temp) { 
        // Assign old elements equal to the new ones
        dim1 = temp.dim1;
        dim2 = temp.dim2;
        dim3 = temp.dim3;
        dim4 = temp.dim4;
        dim5 = temp.dim5;
        dim6 = temp.dim6;
        length = temp.length;
        // Allocate new array
        //this_array = new T[length];
        //this_array = (T *) Kokkos::kokkos_malloc<Kokkos::CudaSpace>(length*sizeof(T));
        this_array = TArray1D("this_array", length);
    }
    // Return *this
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
k_test_t<T>::~k_test_t () {
    //delete[] this_array;
    //Kokkos::kokkos_free(this_array);
}


#endif //MATAR_H
