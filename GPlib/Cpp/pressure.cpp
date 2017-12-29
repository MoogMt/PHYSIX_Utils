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
  int step_press = 0;
  //----------------------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while(  readPressure( press_in, pressure ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  press_vec.push_back( pressure );
	  step_press++;
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  //------------
  // Histograms
  //---------------------------------------------------------------
  // Technical values
  double hist_start = min( press_vec ) ;
  double hist_end   = max( press_vec );
  int nb_box        = 750;
  //---------------------------------------------------------------
  writeHistogram( press_out , normalizeHistogram( makeRegularHistogram( press_vec , hist_start , hist_end , nb_box ) ) );
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
