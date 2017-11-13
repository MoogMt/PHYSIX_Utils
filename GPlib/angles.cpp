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
  std::ofstream oangles("oangles.dat");
  //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 2000; // Start step
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
  int nb_box = 1000;
  //-----------------------------------
  std::vector<double> c_angles;
  std::vector<double> o_angles;
  std::vector<Bin> c_angles_hist;
  std::vector<Bin> o_angles_hist;
  //------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 )
	{
	  // Makes the contact matrix
	  makeContactMatrix( cm_connection , cm_distance , atom_list, cell , cut_off , lut_list );
	  // Making molecules
	  std::vector<Molecule> molecules = makeMolecules( cm_connection );
	  for ( int j=0 ; j < molecules.size() ; j++ )
	    {
	      for ( int i=0 ; i < molecules[j].atom_index.size() ; i++ )
		{
		  if ( i < 32 )
		    {
		      appendVector( c_angles, getAngleAtom( cm_distance, molecules[j] , molecules[j].atom_index[i] ) ) ;
		    }
		  else
		    {
		      appendVector( o_angles , getAngleAtom( cm_distance, molecules[j] , molecules[j].atom_index[i] ) ) ;
		    }
		}
	    }
	  // Step making
	  std::cout << step << std::endl;
	}      
      step++;
     }
  //----------------------------------------------------

  //-------------------
  // Making Histograms
  //----------------------------------------------------
  c_angles_hist = makeRegularHistogram( c_angles, hist_start , hist_end , nb_box);
  o_angles_hist = makeRegularHistogram( o_angles, hist_start , hist_end , nb_box);
  //----------------------------------------------------

  //--------------------
  // Writting histograms
  //----------------------------------------------------
  writeHistogram( cangles , normalizeHistogram( c_angles_hist ) );
  writeHistogram( oangles , normalizeHistogram( o_angles_hist ) );
  //----------------------------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  cangles.close();
  oangles.close();
  //----------------------
  
  return 0;
}
