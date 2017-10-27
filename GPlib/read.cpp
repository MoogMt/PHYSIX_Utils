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
  int comp_step = 15;
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
      if ( step % comp_step == 0  ) //% comp_step == 0 )
	{
	  std::cout << "putain de " << step << std::endl;
	  ContactMatrix cm =  makeContactMatrix ( atom_list, cell , cut_off , lut_list );
	  std::cout << "salete" << step << std::endl;
	  std::vector<MoleculeBasic> mols = makeMolecules( cm );
	  std::cout << "de merde." << step << std::endl;
	  int count=0;
	  int count2=0;
	  /*std::cout << count << " "<<  count2 << " " << mols.size() << std::endl;
	  for ( int i=0 ; i < mols.size() ; i++ )
	    {
	      count += mols[i].names.size();
	      count2 += mols[i].atom_index.size();
	      }*/
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
