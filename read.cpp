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
  std::ofstream outputC_3nn ("3nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_4nn ("4nearestC.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  std::ofstream outputO_2nn ("2nearestO.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------------------------
  Cell box = {9.0,9.0,9.0,90,90,90};    // Definition of simulation box
  double cut_off_radius = 1.6;          // Cut-Off for molecules
  int step = 1;                         // Step counter
  int comp_step=2;                      // The number of step you wait to compute CM
  std::vector<Atom> atom_list;          // Atoms in cell
  std::vector<int> atom_indexesC;       // Indexes of the atoms to study
  std::vector<int> atom_indexesO;       // Indexes of the atoms to study
  Contact_Matrix contact_matrix_init;   // Initial Contact Matrix
  std::vector<Bin> hist_3C;
  std::vector<Bin> hist_4C;
  std::vector<Bin> hist_2O;
  //----------------------------------------------------------------------------------

  //-------------------
  // Reading XYZ file
  //---------------------------------------------------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if ( step == 1 )
	{
	  atom_indexesC = makeVec(0,32);
	  atom_indexesO = makeVec(32,atom_list.size());
	  contact_matrix_init = makeContactMatrix( atom_list, box );
	  hist_3C = makeRegularHistogram( getNNearest( contact_matrix_init , 3 , atom_indexesC ) , 1.0 , 4.0 , 50 );
	  hist_4C = makeRegularHistogram( getNNearest( contact_matrix_init , 4 , atom_indexesC ) , 1.0 , 4.0 , 50 );
	  hist_2O = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesO ) , 1.0 , 4.0 , 50 );
	}
      if( step % comp_step == 0 && !(debug) && step != 1 )
	{
	  //----------------
	  // Contact Matrix
	  //-------------------------------------------------------------------
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  //-------------------------------------------------------------------
	  // Nearest Neighbours
	  //--------------------------------------------------------------------
	  hist_3C = addHistograms( hist_3C, makeRegularHistogram( getNNearest( contact_matrix_init , 3 , atom_indexesC ) , 1.0 , 4.0 , 50 ) );
	  hist_4C = addHistograms( hist_4C, makeRegularHistogram( getNNearest( contact_matrix_init , 4 , atom_indexesC ) , 1.0 , 4.0 , 50 ) );
	  hist_2O = addHistograms( hist_2O, makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesO ) , 1.0 , 4.0 , 50 ) );
	  //--------------------------------------------------------------------
	  // Print step
	  //-----------------------------------------
	  std::cout << "step " << step << std::endl;
	  //-----------------------------------------
	  }
      step++;
    } while( atom_list.size() != 0 );
  //-----------------------------------------------------------------------------

  writeHistogram( outputC_3nn , hist_3C );
  writeHistogram( outputC_4nn , hist_4C );
  writeHistogram( outputO_2nn , hist_2O );
  
  //Closing fluxes
  //----------------------
  input.close();
  outputC_3nn.close();
  outputC_4nn.close();
  outputO_2nn.close();
  //----------------------
  
  return 0;
}
