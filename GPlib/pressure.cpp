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
  //--------------------------------------
  int step       = 1;    // Step counter
  int start_step = 2000; // Start step
  int end_step = 40000;
  int comp_step  = 1;    // Frequency of computation
  //--------------------------------------

  //---------
  // Pressure
  //-------------------------------------------------------------------
  int step_press = 200;
  int stride_press = 0;
  double gen_avg_press = 0;
  double gen_var_press = 0;
  double loc_avg_press = 0;
  double loc_var_press = 0;
  //-------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  double pressure = 0;
  while(  readPressure( press_in, pressure ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  gen_avg_press += pressure;
	  gen_var_press += pressure*pressure;
	  loc_avg_press += pressure;
	  loc_var_press += pressure*pressure;
	  if ( step_press % stride_press == 0 )
	    {
	      loc_avg_press /= (double)(stride_press);
	      loc_var_press = loc_var_press/(double)(stride_press) - loc_avg_press*loc_avg_press;
	      press_out << step << " " << gen_avg_press << " " << loc_var_press << " " << sqrt( loc_var_press ) << std::endl;
	      std::cout << step << " " << gen_avg_press << " " << loc_var_press << " " << sqrt( loc_var_press ) << std::endl;
	      loc_avg_press = 0;
	      loc_var_press = 0;
	    }
	  step_press++;
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  // Writting results
  gen_avg_press = gen_avg_press/(end_step-start_step);
  gen_var_press = gen_var_press/(end_step-start_step) - gen_avg_press*gen_avg_press;
  std::cout << gen_avg_press << " " << gen_var_press << " " << sqrt(gen_var_press) << std::endl;

  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  press_in.close();
  press_out.close();
  //----------------------
  
  return 0;
}
