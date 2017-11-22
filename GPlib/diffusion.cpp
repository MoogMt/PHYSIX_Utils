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
  std::vector<double> r; // Atoms_old in cell
  AllTypeLUT lut_list; // LUT for types
  std::vector<double> x0, y0, z0;
  std::vector<double> x, y, z;
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

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //------------------------------------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = atom_list.x;
	  y0 = atom_list.y;
	  z0 = atom_list.z;
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( atom_list.x , x0 );
	  std::vector<double> y = difference( atom_list.y , y0 );
	  std::vector<double> z = difference( atom_list.z , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  diffusion << step << " " << average( r ) << std::endl;
	}      
      std::cout << step << std::endl;
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
