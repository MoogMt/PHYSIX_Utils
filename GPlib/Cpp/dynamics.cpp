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
  std::ofstream graph_molecules ( "graph_molecules.dat" ,  std::ios::out );
  //-------------------------------------------------------------------------------

  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 5000; // Start step
  int end_step   = 1000000; // End step
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

  //--------------  
  // Atom Indexes
  //------------------------------
  std::vector<int> atom_indexesC;
  std::vector<int> atom_indexesO;
  //------------------------------
  
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

  //------------------------------
  std::vector<double> coordinanceC;
  std::vector<double> coordinanceO;
  std::vector<double> lifetime;
  std::vector<double> lifetime_C;
  std::vector<double> lifetime_C2;
  std::vector<double> lifetime_C3;
  std::vector<double> lifetime_C4;
  std::vector<double> lifetime_O;
  std::vector<double> lifetime_O1;
  std::vector<double> lifetime_O2;
  std::vector<double> changes_all;
  std::vector<double> changes_C;
  std::vector<double> changes_O;
  //------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      if ( step % comp_step == 0 && step >= start_step && step < end_step)
	{
	  makeContactMatrix( cm_distance , atom_list, cell , cut_off , lut_list );
	  if ( step == start_step )
	    {
	      nb_atoms = atom_list.names.size();
	      // Get index of species
	      atom_indexesC = getSpecieIndex( lut_list, "C" );
	      atom_indexesO = getSpecieIndex( lut_list, "O" );
	      // Initializing
	      lifetime.assign( nb_atoms , 0. );
	      coordinance.assign( nb_atoms, 0. );
	       for ( int i=0 ; i < atom_indexesC.size() ; i++ )
		{
		  coordinanceC[i] = cumSum( getAtomContact( cm , atom_indexesC[i] ) );
		}
	      for ( int i=0 ; i < atom_indexesO.size() ; i++ )
		{
		  coordinanceO[i] = cumSum( getAtomContact( cm , atom_indexesO[i] ) );
		}
	    }
	  else
	    {
	      int changes = 0; int c_changes=0; int o_changes=0;
	      for ( int i=0 ; i < atom_indexesC.size() ; i++ )
		{
		  double sum = cumSum( getAtomContact( cm , atom_indexesC[i] ) );
		  if( (sum-2.5)*(coordinanceC[i]-2.5) < 0  )
		    {
		      
		    }
		  else
		    {
		      
		    }
		}
	      for ( int i=0 ; i < atom_indexesO.size() ; i++ )
		{
		  double sum = cumSum( getAtomContact( cm , atom_indexesO[i] ) );
		}
	    }
	  std::cout << step << std::endl;
	}
      step++;
    }
  //----------------------------------------------------
  
  
  return 0;
}
