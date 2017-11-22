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
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream diffusion("diffusion.dat");
  //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 2000; // Start step
  int end_step   = 120000;
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;     // Atoms in cell
  std::vector<double> r, r_0; // Atoms_old in cell
  double origin[3] = { 0 , 0 , 0};
  AllTypeLUT lut_list; // LUT for types
  //--------------------------------------------------

  //--------------------
  // Reading Cell File
  //-------------------------------------------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step ) r_0 = calculate_distance( atom_list , origin );
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> r = calculate_distance( atom_list , origin) ;
	  double r2_avg = abs( average( difference( r , r_0 ) ) ); r2_avg *= r2_avg;
	  diffusion << step << " " << r2_avg << std::endl;
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  
  //---------------
  //Closing fluxes
  //----------------------
  input.close();
  diffusion.close();
  //----------------------
      
  return 0;
}
