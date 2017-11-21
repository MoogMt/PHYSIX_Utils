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
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream output4CO_1st( "4CO_1st.dat" ,  std::ios::out );
  std::ofstream output4CO_2nd( "4CO_2nd.dat" ,  std::ios::out );
  //--------------------------------------------------------------------------------
  
  //------
  // Data
  //-----------------------------------------
  std::vector<std::vector<BinReal> > hist_list;
  //-----------------------------------------

  //--------------
  // Reading Data
  //-------------------------------------------------------------------
  // 3CO
  //------------------------------------------------
  // 2500K
  hist_list.push_back( readRegularHistogramReal( "3500K/40GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "2500K/45GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "2500K/50GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "2500K/60GPa/4nearestCO.dat" ) );
  // 3000K
  hist_list.push_back( readRegularHistogramReal( "3000K/40GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3000K/45GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3000K/50GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3000K/60GPa/4nearestCO.dat" ) );
  // 3500K
  hist_list.push_back( readRegularHistogramReal( "3500K/40GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3500K/45GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3500K/50GPa/4nearestCO.dat" ) );
  hist_list.push_back( readRegularHistogramReal( "3500K/60GPa/4nearestCO.dat" ) );
  //-------------------------------------------------------------------i

  //----------
  // Pressure
  //-------------------------------------------------------------------
  double press[4] = { 40, 45 , 50, 60 };
  //-------------------------------------------------------------------

  //-----------------
  // Writting data
  //-------------------------------------------------------------------
  for ( int i=0 ; i < 4 ; i++ )
    {
      output4CO_1st << press[i] << " " << integrateHistogram( hist_list[i] , 0 , 1.75 ) << " " << integrateHistogram( hist_list[i+4] , 0 , 1.75 ) << " " << integrateHistogram( hist_list[i+8] , 0 , 1.75 )  <<  std::endl;
      output4CO_2nd << press[i] << " " << integrateHistogram( hist_list[i] , 1.75 , 3.00 ) << " " << integrateHistogram( hist_list[i+4] , 1.75 , 3.00 ) << " " << integrateHistogram( hist_list[i+8] , 1.75 , 3.00 )  <<  std::endl;
    }
  //-------------------------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  output4CO_1st.close();
  output4CO_2nd.close();
  //----------------------
  
  return 0;
}
