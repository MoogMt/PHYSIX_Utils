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
int main(void)
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //---------------------------------

  //----------------------
  // Physical parameters
  //--------------------------------------
  int step      = 1;  // Step counter
  int comp_step = 1;
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;               // Atoms in cell
  AllTypeLUT lut_list ; // LUT for types
  //---------------------------------------------------

  //--------------------
  // Reading Parameters
  //-------------------------------------------------------------------
  Cell         cell    = readParamCell ( "cell.param" );
  CutOffMatrix cut_off = readCutOff    ( "cut_off.dat" , lut_list );
  //-------------------------------------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZ( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0  ) 
	{
	  ContactMatrix cm =  makeContactMatrix ( atom_list, cell , cut_off , lut_list );
	  std::vector<MoleculeBasic> mols = makeMolecules( cm );
	  for ( int i=0 ; i < mols.size() ; i++ )
	    {
	      std::cout << "-------------------------------------------" << std::endl;
	      for ( int j=0 ; j < mols[i].names.size() ; j++ )
		{
		  std::cout << mols[i].names[j] << " " << mols[i].atom_index[j] << std::endl;
		}
	      std::cout << "-------------------------------------------" << std::endl;
	    }
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
