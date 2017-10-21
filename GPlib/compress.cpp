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
//-------------------------

//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //---------------------------------

  //--------
  // Output
  //-------------------------------------------------------------
  std::ofstream compStruc("compressedBox.xyz",  std::ios::out );
  //-------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------
  int step = 1;                            // Step counter
  double frac_a = 0.98; double frac_b = 0.98;   double frac_c = 0.98;
  //----------------------------------------------------------------

  //---------------
  // Initializers
  //----------------------------------------------
  std::vector<Atom> atom_list;   // Atoms in cell
  //----------------------------------------------

  //---------------
  // Reading cell
  //---------------------------------
  Cell box=readParamCell("cell.param");
  //---------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if ( step == 8000 )
	{
	  compressBox( { wrapPBC(atom_list,box) , box } , frac_a , frac_b , frac_c);
	  writePositions( compStruc , atom_list , "C" );
	  writePositions( compStruc , atom_list , "O" );
	}
      step++;
    } while( atom_list.size() != 0 );
  //----------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
  compStruc.close();
  //----------------------
  
  return 0;
}
