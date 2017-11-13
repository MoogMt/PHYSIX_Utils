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
  std::ofstream cangles("cangles.dat");
  //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 2000; // Start step
  int end_step = 22000;
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;  // Atoms in cell
  AllTypeLUT lut_list; // LUT for types
  ContactMatrix cm_connection;    // Contact Matrix
  ContactMatrix cm_distance;    // Contact Matrix
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
  //-------------------------------------------------------------------

  //------------------------------------
  // Histograms
  //------------------------------------
  // Technical values
  double hist_start =  0.0;
  double hist_end   =  180.0;
  int nb_box = 1800;
  //-----------------------------------
  std::vector<Bin> c_angles_hist;
  //------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  // Makes the contact matrix
	  makeContactMatrix( cm_connection , cm_distance , atom_list, cell , cut_off , lut_list );
	  // Making molecules
	  std::vector<Molecule> molecules = makeMolecules( cm_connection );
	  // Calculating angles
	  std::vector<double> c_angles;
	  for ( int i=0 ; i < molecules.size() ; i++ )
	    {
	      for ( int j=0 ; j < molecules[i].atom_index.size() ; j++ )
		{
		  if ( molecules[i].names[j] == "C" )
		    {
		      appendVector( c_angles , getAngleAtom( cm_distance, molecules[i] , molecules[i].atom_index[j] ) );
		    }
		}
	    }
	  std::vector<Bin> cangles_hist_temp = makeRegularHistogram( c_angles, hist_start , hist_end , nb_box);
	  // Building up the histogram
	  if ( step > start_step + 1 )
	    {
	      c_angles_hist = addHistograms( c_angles_hist , cangles_hist_temp );
	    }
	  else
	    {
	      c_angles_hist = cangles_hist_temp;
	    }
	  // Step making
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  //--------------------
  // Writting histograms
  //----------------------------------------------------
  writeHistogram( cangles , normalizeHistogram( c_angles_hist ) );
  //----------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  cangles.close();
  //----------------------
  
  return 0;
}
