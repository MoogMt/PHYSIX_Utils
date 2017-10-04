#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "utils.h"

//------
// BINS
//----------------------
struct Bin
{
  double begin, end; // Inner and Outer bound of the bin
  int value;         // Number of elements inside the bin
};
//----------------------


// FUNCTIONS
//------------------------------------------------------------------------------
//Bins
Bin emptyBin();
double center(Bin bin);
bool overlap( Bin bin1, Bin bin2 );
void fillBin( Bin &bin , std::vector<double> data );
Bin makeBin( int bin_min, int bin_max, std::vector<double> data );
Bin addBinsMin( Bin bin1, Bin bin2);
Bin addBinsMax( Bin bin1, Bin bin2 );
// Histograms
// Regular histogram
std::vector<Bin> makeRegularHistogram( std::vector<double> data , double x_min , double x_max , int number_bins );
std::vector<Bin> makeRegularHistogram( std::vector<double> data_x , std::vector<double> data_y , int number_bins );
// Custom histogram
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<double> bins_limits);
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<Bin> bins);
std::vector<Bin> addHistograms( std::vector<Bin> hist1 , std::vector<Bin> hist2 );
// Write histogram
void writeHistogram( std::ofstream & file , std::vector<Bin> hist );
//------------------------------------------------------------------------------

#endif 
