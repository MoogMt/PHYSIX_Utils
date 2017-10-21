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


// Min
//------------------------------------------
int min( int int1, int int2 );
double min ( double real1, double real2 );
double min(std::vector<double> vector);
//-------------------------------------------

// Max
//----------------------------------------
int max( int int1, int int2 );
double max ( double real1, double real2 );
double max(std::vector<double> vector);
//----------------------------------------

//----------
// Averages
//-------------------------------------------
double average(std::vector<int> data);
double average(std::vector<double> data);
//-------------------------------------------

// Sums
//---------------------------------
int sumFromO(int integer);
int sumBtw(int int1, int int2);
//---------------------------------

// Compute separator
//------------------------------------------------
int computeSep(int atom_index, int nb_atoms);
//------------------------------------------------

// Utils for vectors
//----------------------------------------------------------------------------
std::vector<int> makeVec(int init, int final);
std::vector<double> sortVector(std::vector<double> to_sort, bool increasing);
std::vector<double> sortVectorIncreasing(std::vector<double> unsorted);
std::vector<double> sortVectorDecreasing(std::vector<double> unsorted);
//----------------------------------------------------------------------------

// Conversion
//----------------------------------------------------------------------------
double it2real (  std::istream_iterator<std::string> iterator );
//----------------------------------------------------------------------------
  
#endif 
