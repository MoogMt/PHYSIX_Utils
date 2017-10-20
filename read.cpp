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

  //----------------------
  // Physical parameters
  //---------------------------------------------------------------
  int step      = 1;  // Step counter
  int comp_step = 5;
  //----------------------------------------------------------------

  //---------------
  // Initializers
  //----------------------------------------------
  std::vector<Atom> atom_list;   // Atoms in cell
  std::vector<typeLUT> list_lut;
  //----------------------------------------------

  //---------------
  // Reading cell
  //---------------------------------
  Cell box=readParamCell("cell.param");
  //---------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZ( input , atom_list , list_lut ) )
    {
      if ( step % comp_step == 0 )
	{
	  ContactMatrix cm = makeContactMatrix( atom_list , list_lut , box );
	  std::cout << "step: " << step << std::endl;
	}
      step++;
     }
  //----------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
   //----------------------
  
  return 0;
}
