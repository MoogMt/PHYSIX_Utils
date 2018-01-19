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
  std::ifstream press_in("STRESS");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream press_out("pressure.dat");
  //-------------------------------------

  //------------------------
  // Checking if file exists
  //-------------------------------------------------------------
  if ( press_out )
    {
      std::cout << "File STRESS open" << std::endl;
    }
  else
    {
      std::cout << "File STRESS does not exists!" << std::endl;
      return 0;
    }
  //-------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------
  int step       = 1;    // Step counter
  int start_step = 5000; // Start step
  int end_step = 2000000;
  int comp_step  = 1;    // Frequency of computation
  //----------------------------------------------------

  //-------
  // Press 
  //----------------------------------------------------
  double pressure; std::vector<double> press_vec ;
  std::vector<Bin> press_hist;
  //----------------------------------------------------

  //-----------
  // Tolerance
  //----------------------------------------------------
  double tolerance = -1 ; 
  while ( tolerance > 1 || tolerance < 0 )
    {
      std::cout << "Tolerance: ";
      std::cin >> tolerance;
    }
  //----------------------------------------------------
    
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while(  readPressure( press_in, pressure ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step < end_step )
	{
	  press_vec.push_back( pressure );
	  std::cout << "Computing pressure, step: " << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  //------------
  // Histogram
  //---------------------------------------------------------------
  writeHistogram( press_out , normalizeHistogram( makeRegularHistogram( press_vec , min( press_vec ) , max( press_vec) , 200 ) ) );
  //---------------------------------------------------------------

  //----------------------------------------------
  // Print average and variance of the pressure 
  //---------------------------------------------------------------
  std::cout << "Block Average Pressure: " << blockaverage( press_vec , tolerance , 100, 100) << std::endl;
  std::cout << "Average Pressure: " << average( press_vec )  << std::endl;
  //---------------------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  press_in.close();
  press_out.close();
  //----------------------
  
  return 0;
}
