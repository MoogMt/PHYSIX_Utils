#ifndef CONTACT_MATRIX_H
#define CONTACT_MATRIX_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

//------
// BINS
//----------------------
struct Bin
{
  double begin, end;
  int value;
};
//----------------------


// FUNCTIONS
//------------------------------------------------------------------------------
//Bins 
void fillBin( Bin & bin , std::vector<double> data );
Bin makeBin( int bin_min, int bin_max, std::vector<double> & data );
// Histograms
std::vector<bins> makeRegularHistogram( std::vector<double> data_x , std::vector<double> data_y , double x_min , double x_max , int number_bins )
//------------------------------------------------------------------------------

#endif 
