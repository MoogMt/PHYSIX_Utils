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
  std::ofstream output("velocity.dat");
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
  AtomList  atom_list;  // Atoms in cell
  AtomList  atom_list_prev;  // Atoms in cell
  AllTypeLUT lut_list; // LUT for types
  //--------------------------------------------------

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  
	  atom_list_prev = atom_list;
	  std::cout << step << std::endl;
	}
      else if ( step == start_step )
	{
	  atom_list_prev = atom_list;
	}
      step++;
     }
  //----------------------------------------------------
  
  //--------------------
  // Making histograms
  //----------------------------------------------------
  c2_angles_hist = makeRegularHistogram( c2_angles , hist_start , hist_end , nb_box );
  writeHistogram( c2_angles_out , normalizeHistogram( c2_angles_hist ) );
  c3_angles_hist = makeRegularHistogram( c3_angles , hist_start , hist_end , nb_box );
  writeHistogram( c3_angles_out , normalizeHistogram( c3_angles_hist ) );
  c4_angles_hist = makeRegularHistogram( c4_angles , hist_start , hist_end , nb_box );
  writeHistogram( c4_angles_out , normalizeHistogram( c4_angles_hist ) );
  o2_angles_hist = makeRegularHistogram( o2_angles , hist_start , hist_end , nb_box );
  writeHistogram( o2_angles_out , normalizeHistogram( o2_angles_hist ) );
  o3_angles_hist = makeRegularHistogram( o3_angles , hist_start , hist_end , nb_box );
  writeHistogram( o3_angles_out , normalizeHistogram( o3_angles_hist ) );
  //----------------------------------------------------

  // Print other values
  //----------------------------------------------------
  std::cout << "Other C Angles: " << std::endl;
  for ( int i=0 ; i < c_others_nb.size() ; i++ )
    {
      std::cout << "value: " << c_others_nb[i] << std::endl;
    }
  //----------------------------------------------------
  std::cout << "Other O Angles: " << std::endl;
  for ( int i=0 ; i < o_others_nb.size() ; i++ )
    {
      std::cout << "value: " << o_others_nb[i] << std::endl;
    }
  //----------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  c2_angles_out.close();
  c3_angles_out.close();
  c4_angles_out.close();
  o2_angles_out.close();
  o3_angles_out.close();
  //----------------------
      
  return 0;
}
