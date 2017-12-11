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
#include <sstream>

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
struct BinReal
{
  double begin, end; // Inner and Outer bound of the bin
  double value;         // Number of elements inside the bin
};
//----------------------

//===========
// FUNCTIONS
//==============================================================================

//------
// Bins
//------------------------------------------------------------------------------
Bin emptyBin();
double center(Bin bin);
Bin makeBinReal( double begin , double end , int value );
bool overlap( Bin bin1, Bin bin2 );
void fillBin( Bin &bin , std::vector<double> data );
Bin makeBin( int bin_min, int bin_max, std::vector<double> data );
Bin addBinsMin( Bin bin1, Bin bin2);
Bin addBinsMax( Bin bin1, Bin bin2 );
//------------------------------------------------------------------------------

//----------------
// Histograms
//------------------------------------------------------------------------------
// Regular histogram
//------------------------------------------------------------------------------
std::vector<Bin> makeRegularHistogram( std::vector<double> data , double x_min , double x_max , int number_bins );
std::vector<Bin> makeRegularHistogram( std::vector<double> data_x , std::vector<double> data_y , int number_bins );
// Custom histogram
//------------------------------------------------------------------------------
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<double> bins_limits);
std::vector<Bin> makeHistograms( std::vector<double> data, std::vector<Bin> bins);
std::vector<Bin> addHistograms( std::vector<Bin> hist1 , std::vector<Bin> hist2 );
double getTotalValue( std::vector<Bin> hist );
void writeHistogram( std::ofstream & file , std::vector<Bin> hist );
void writeHistogram( std::string file_name , std::vector<Bin> hist );
//------------------------------------------------------------------------------

//-----------
// Bin Real
//------------------------------------------------------------------------------
BinReal emptyBinReal();
double center( BinReal bin );
bool overlap( BinReal bin1 , BinReal bin2 );
BinReal makeBinReal( double begin , double end , double value );
void fillBin( BinReal &bin , std::vector<double> data );
BinReal makeBinReal( double bin_min, double bin_max, std::vector<double> dat);
BinReal addBinsMin( BinReal bin1, BinReal bin2);
BinReal addBinsMax( BinReal bin1, BinReal bin2 );
//------------------------------------------------------------------------------

//-----------------
// Real Histogram
//------------------------------------------------------------------------------
std::vector<BinReal> normalizeHistogram( std::vector<Bin> hist );
void writeHistogram( std::ofstream & file , std::vector<BinReal> hist );
void writeHistogram( std::string file_name , std::vector<BinReal> hist );
bool checkSizeHists( std::vector<BinReal> hist , std::vector<BinReal> hist2 );
bool checkSizeHists( std::vector< std::vector<BinReal> > hist_list );
void writeBinReal( std::ofstream & file , BinReal bin, bool wcenter);
void writeBinRealCenter( std::ofstream & file, BinReal bin );
void writeHistBinCenter( std::ofstream & file , std::vector<BinReal> hist , int index );
void writeHistBin( std::ofstream & file , std::vector<BinReal> hist, int index , bool wcenter );
void writeHistograms( std::ofstream & file , std::vector< std::vector<BinReal> > hist_list );
//------------------------------------------------------------------------------

//-----------------------
// Integrate Histograms
//-------------------------------------------------------------------------------------
double average( std::vector<BinReal> & histogram );
double integrateHistogram( std::vector<BinReal> & histogram );
double integrateHistogram( std::vector<BinReal> & histogram , double start , double end );
//-------------------------------------------------------------------------------------

//----
// IO
//-------------------------------------------------------------------------------------
void readRegularHistogram( std::string file_name , std::vector<Bin> & histogram );
std::vector<Bin> readRegularHistogram( std::string file_name );
std::vector<BinReal> readRegularHistogramReal( std::string file_name );
//-------------------------------------------------------------------------------------

//======================================================================================

#endif 
