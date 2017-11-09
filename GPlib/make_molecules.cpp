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
  //---------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream molecules ("molecules.dat",  std::ios::out );
  //--------------------------------------------------------------------------------
  
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

  //-----------
  // Histogram
  //-----------------------------------------
  double hist_start = 0.5, hist_end = 96.5;
  int nb_box = 96;
  std::vector<Bin> hist;
  //-----------------------------------------
  
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
	  // Prints the bonds between atoms
	  for( int i=0 ; i < molecules.size() ; i++ )
	    {
	      for ( int j=0 ; j < molecules[i].bonds.size(); j++ )
		{
		  std::cout << molecules[i].bonds[j].atom1_index << " " << molecules[i].bonds[j].atom2_index << std::endl;
		}
	      std::cout << "------------" << std::endl;
	    }
	  std::cout << "================" << std::endl;
	  // Stock the size of the molecules in the box
	  std::vector<double> sizes;
	  for ( int i=0 ; i < molecules.size() ; i++ )
	    {
	      sizes.push_back( molecules[i].names.size() );
	    }
	  // Make an histogram of the sizes of the molecules in the box
	  if ( step == 1 )
	    {
	      hist = makeRegularHistogram( sizes , hist_start , hist_end , nb_box );
	    }
	  else
	    {
	      hist = addHistograms ( hist , makeRegularHistogram( sizes , hist_start , hist_end , nb_box ) );
	    }
	}      
      step++;
     }
  //----------------------------------------------------
  // Writes histogram
  writeHistogram( molecules , normalizeHistogram( hist ) );
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
  molecules.close();
   //----------------------
  
  return 0;
}
