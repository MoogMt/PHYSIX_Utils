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
  //-------
  // DEBUG
  //-----------------
  bool debug=false;
  //-----------------
  
  //---------
  // Input
  //------------------------------------
  std::ifstream input("TRAJEC.xyz");
  //------------------------------------

  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream trackCC("trackCC.dat",  std::ios::out );
  std::ofstream CCxyz("CC-molecules.xyz",  std::ios::out  | std::ios::app );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------------------------
  double cut_off_radius = 1.6;             // Cut-Off for molecules
  int step = 1;                            // Step counter
  int comp_step=2;                         // The number of step you wait to compute CM
  int start_step = 5000, end_step = 31000; // Start and end step for datanalysis
  double hist_start  = 0.95;  double hist_end = 3.00; int nb_box = 300;
  //----------------------------------------------------------------------------------

  //---------
  // Initializers
  //----------------------------------------------------------------------------------
  std::vector<Atom> atom_list;                                    // Atoms in cell
  std::vector<int> atom_indexesC; std::vector<int> atom_indexesO; // Indexes of atoms
  Contact_Matrix contact_matrix_init;                             // Initial Contact Matrix
  //----------------------------------------------------------------------------------

  //---------------
  // Reading cell
  //---------------------------------
  Cell box=readParamCell("cell.param");
  //---------------------------------
  
  //-------------------
  // Reading XYZ file
  //---------------------------------------------------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if ( step == 1 )
	{
	  atom_indexesC = makeVec(0,31);
	  atom_indexesO = makeVec(32,atom_list.size());
	}
      if( step % comp_step == 0 && !(debug) && step >= start_step && step <=  end_step )
	{
	  //----------------
	  // Contact Matrix
	  //-------------------------------------------------------------------
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  //-------------------------------------------------------------------
	  for ( int i=0 ; i < atom_indexesC.size()-1 ; i++ )
	    {
	      for (int j=i+1 ; j < atom_indexesC.size() ; j++ )
		{
		  if ( getDistance(contact_matrix, atom_indexesC[i],  atom_indexesC[j] ) < 1.6 )
		    {
		      
		      trackCC << step << " " << i << " " << j << " " << getAtomNeighboursNb(contact_matrix, i , "O", cut_off_radius ) << " " << getAtomNeighboursNb(contact_matrix, j , "O", cut_off_radius ) << " " <<  getAtomNeighboursNb(contact_matrix, i , "C", cut_off_radius ) << " " << getAtomNeighboursNb(contact_matrix, j , "C", cut_off_radius ) << std::endl;
		      writeXYZ(CCxyz,atom_list);
		    }
		}
	    }
	  //--------------------------------------------------------------------
	  //-----------------------------------------
	  std::cout << "step " << step << std::endl;
	  //-----------------------------------------
	  }
      step++;
    } while( atom_list.size() != 0 );
  //-----------------------------------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
  trackCC.close();
  CCxyz.close();
  //----------------------
  
  return 0;
}
