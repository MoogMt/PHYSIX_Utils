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
  std::ifstream input("input.xyz");
  //---------------------------------

  //--------
  // Output
  //-------------------------------------------------------------
  std::ofstream C_position( "C_position.cpmd" , std::ios::out );
  std::ofstream O_position( "O_position.cpmd" , std::ios::out );
  std::ofstream xyz( "position.xyz" , std::ios::out );
  //-------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------
  double frac_a = 0.918367347; double frac_b = 0.918367347; double frac_c =  0.918367347;
  //----------------------------------------------------------------

  //---------------
  // Initializers
  //----------------------------------------------
  AtomList atom_list;
  AllTypeLUT lut_list;
  //----------------------------------------------

  //---------------
  // Reading cell
  //---------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) ) return 1;
  //---------------------------------

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) ) return 1;
  //-------------------------------------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  if ( readStepXYZfast( input , atom_list , lut_list , false , false) )
    {
      wrapPBC( atom_list , cell) ;
      compressCell( atom_list , cell , frac_a , frac_b , frac_c );
      writeXYZ( xyz , atom_list );
      writePositions( C_position , atom_list , "C" );
      writePositions( O_position , atom_list , "O" );
    }
  //----------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
  C_position.close();
  O_position.close();
  xyz.close();
  //----------------------
  
  return 0;
}
