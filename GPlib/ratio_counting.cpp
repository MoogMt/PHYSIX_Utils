//---------------
// External Libs
//-------------------
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>
//---------------------

//----------------
// Internal Files
//-------------------------
#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "contact_matrix.h"
#include "xyz.h"
#include "histogram.h"
#include "pdb.h"
#include "sim.h"
#include "lut.h"
#include "molecules.h"
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main( void )
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //---------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream output( "ratio.dat" ,  std::ios::out );
  //--------------------------------------------------------------------------------
  
  //------
  // Data
  //-----------------------------------------
  std::vector<Bin> hist1CC;
  std::vector<Bin> hist2CC;
  std::vector<Bin> hist1CO;
  std::vector<Bin> hist2CO;
  std::vector<Bin> hist3CO;
  std::vector<Bin> hist4CO;
  std::vector<Bin> hist1OO;
  //-----------------------------------------

  // Reading Data
  //-------------------------------------------------------------------
  readRegularHistogram( "1CCnearest.dat" , hist1CC );
  readRegularHistogram( "1CCnearest.dat" , hist2CO );
  readRegularHistogram( "1CCnearest.dat" , hist3CO );
  readRegularHistogram( "1CCnearest.dat" , hist4CO );
  readRegularHistogram( "1CCnearest.dat" , hist1CC );
  readRegularHistogram( "1CCnearest.dat" , hist2CC );
  readRegularHistogram( "1CCnearest.dat" , hist1OO );
  //-------------------------------------------------------------------

  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  output.close();
  //----------------------
  
  return 0;
}
