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
  std::ofstream outputC_1nn ("1nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_2nn ("2nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_3nn ("3nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_4nn ("4nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputO_1nn ("1nearestO.dat",  std::ios::out | std::ios::app );
  std::ofstream outputO_2nn ("2nearestO.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------------------------
  double cut_off_radius = 1.6;             // Cut-Off for molecules
  int step = 1;                            // Step counter
  int comp_step=2;                         // The number of step you wait to compute CM
  int start_step = 5000, end_step = 30000; // Start and end step for datanalysis
  double hist_start = 0.95;  double hist_end   = 3.50; int nb_box = 100;
  //----------------------------------------------------------------------------------

  //---------
  // Initializers
  //----------------------------------------------------------------------------------
  std::vector<Atom> atom_list;                                    // Atoms in cell
  std::vector<int> atom_indexesC; std::vector<int> atom_indexesO; // Indexes of atoms
  Contact_Matrix contact_matrix_init;                             // Initial Contact Matrix
  std::vector<Bin> hist_1C;  std::vector<Bin> hist_2C;   // Histograms
  std::vector<Bin> hist_3C;  std::vector<Bin> hist_4C;   // Histograms
  std::vector<Bin> hist_1O;  std::vector<Bin> hist_2O;   // Histograms
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
	  atom_indexesC = makeVec(0,32);
	  atom_indexesO = makeVec(32,atom_list.size());
	  contact_matrix_init = makeContactMatrix( atom_list, box );
	  hist_1C = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesC ) , hist_start , hist_end , nb_box );
	  hist_2C = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesC ) , hist_start , hist_end , nb_box );
	  hist_3C = makeRegularHistogram( getNNearest( contact_matrix_init , 3 , atom_indexesC ) , hist_start , hist_end , nb_box );
	  hist_4C = makeRegularHistogram( getNNearest( contact_matrix_init , 4 , atom_indexesC ) , hist_start , hist_end , nb_box );
	  hist_1O = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesO ) , hist_start , hist_end , nb_box );
	  hist_2O = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesO ) , hist_start , hist_end , nb_box );
	}
      if( step % comp_step == 0 && !(debug) && step >= 10000 && step <= 30000 )
	{
	  //----------------
	  // Contact Matrix
	  //-------------------------------------------------------------------
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  //-------------------------------------------------------------------
	  // Nearest Neighbours
	  //--------------------------------------------------------------------
	  hist_1C = addHistograms( hist_1C, makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesC ) , hist_start , hist_end , nb_box ) );
	  hist_2C = addHistograms( hist_2C, makeRegularHistogram( getNNearest( contact_matrix , 2 , atom_indexesC ) , hist_start , hist_end , nb_box ) );
	  hist_3C = addHistograms( hist_3C, makeRegularHistogram( getNNearest( contact_matrix , 3 , atom_indexesC ) , hist_start , hist_end , nb_box ) );
	  hist_4C = addHistograms( hist_4C, makeRegularHistogram( getNNearest( contact_matrix , 4 , atom_indexesC ) , hist_start , hist_end , nb_box ) );
	  hist_1O = addHistograms( hist_1O, makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesO ) , hist_start , hist_end , nb_box ) );
	  hist_2O = addHistograms( hist_2O, makeRegularHistogram( getNNearest( contact_matrix , 2 , atom_indexesO ) , hist_start , hist_end , nb_box ) );
	  //--------------------------------------------------------------------
	  // Print step
	  //-----------------------------------------
	  std::cout << "step " << step << std::endl;
	  //-----------------------------------------
	  }
      step++;
    } while( atom_list.size() != 0 );
  //-----------------------------------------------------------------------------
  writeHistogram( outputC_1nn , normalizeHistogram( hist_1C ) );
  writeHistogram( outputC_2nn , normalizeHistogram( hist_2C ) );
  writeHistogram( outputC_3nn , normalizeHistogram( hist_3C ) );
  writeHistogram( outputC_4nn , normalizeHistogram( hist_4C ) );
  writeHistogram( outputO_1nn , normalizeHistogram( hist_1O ) );
  writeHistogram( outputO_2nn , normalizeHistogram( hist_2O ) );
    
  //Closing fluxes
  //----------------------
  input.close();
  outputC_1nn.close();
  outputC_2nn.close();
  outputC_3nn.close();
  outputC_4nn.close();
  outputO_1nn.close();
  outputO_2nn.close();
  //----------------------
  
  return 0;
}
