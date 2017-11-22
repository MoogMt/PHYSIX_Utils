#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

//=====
// MIN
//============================================================================
int min( int int1, int int2 );
double min ( double real1, double real2 );
double min(std::vector<double> vector);
//============================================================================

//=====
// MAX
//============================================================================
int max( int int1, int int2 );
double max ( double real1, double real2 );
double max(std::vector<double> vector);
//============================================================================

//==========
// AVERAGES
//============================================================================
double average(std::vector<int> data);
double average(std::vector<double> data);
//============================================================================

//======
// SUMS
//============================================================================
int sumFromO(int integer);
int sumBtw(int int1, int int2);
//============================================================================

//===================
// COMPUTE SEPARATOR
//============================================================================
int computeSep(int atom_index, int nb_atoms);
//============================================================================

//=========
// VECTORS
//============================================================================
//--------
// Switch
//---------------------------------------------------------------------
void switchV( std::vector<int> & vector , int index1 , int index2 );
void switchV( std::vector<double> & vector , int index1 , int index2 );
//----------------------------------------------------------------------
// Init
//--------------------------------------------------------
std::vector<int> initVector( int value );
std::vector<double> initVector( double value );
std::vector<std::string> initVector( std::string value );
//---------------------------------------------------------
// Make
//---------------------------------------------
std::vector<int> makeVec(int init, int final);
//---------------------------------------------
// Sort
//-----------------------------------------------------------------------------
std::vector<double> sortVector( std::vector<double> to_sort, bool increasing );
std::vector<double> sortVectorIncreasing( std::vector<double> unsorted );
std::vector<double> sortVectorDecreasing( std::vector<double> unsorted );
//-----------------------------------------------------------------------------
// Appending
//-----------------------------------------------------------------------------
void appendVector( std::vector<double> & vec_to, std::vector<double> vec_from );
//-----------------------------------------------------------------------------
// Average
//-----------------------------------------------------------------------------
double average( std::vector<double> & vector );
//-----------------------------------------------------------------------------
// Basic Operations
//-----------------------------------------------------------------------------
std::vector<double> square( std::vector<double> vector1 );
std::vector<double> squaroot( std::vector<double> vector1 );
std::vector<double> distance( std::vector<double> vector1 , double scalar );
std::vector<double> addVector( std::vector<double> vector1 , std::vector<double> vector2 );
std::vector<double> difference( std::vector<double> vector1 , std::vector<double> vector2 );
//-----------------------------------------------------------------------------
// Products
//-----------------------------------------------------------------------------
double scalarProduct( std::vector<double> vector1 , std::vector<double> vector2 );
std::vector<double> crossProduct( std::vector<double> vector1 , std::vector<double> vector2 );
//-----------------------------------------------------------------------------
// Norm
//-----------------------------------------------------------------------------
double norm( std::vector<double> vector );
//-----------------------------------------------------------------------------
// Distance From Plan
//-----------------------------------------------------------------------------
double getDistanceFromPlan( std::vector<double> vector1, std::vector<double> vector2 , std::vector<double> point_outside , std::vector<double> point_plan );
//============================================================================


//============
// CONVERSION
//============================================================================
double it2real (  std::istream_iterator<std::string> iterator );
//============================================================================

//=======
// ARRAY
//============================================================================
void zeros( int* vec , int nb_atoms );
int sum( int* vector , int nb_atoms );
void copy( int* try2 , int* try1 , int size );
//============================================================================

#endif 
