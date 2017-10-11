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
  std::ofstream outputCC_1nn ("1nearestCC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputCC_2nn ("2nearestCC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputCO_1nn ("1nearestCO.dat",  std::ios::out | std::ios::app );
  std::ofstream outputCO_2nn ("2nearestCO.dat",  std::ios::out | std::ios::app );
  std::ofstream outputCO_3nn ("3nearestCO.dat",  std::ios::out | std::ios::app );
  std::ofstream outputCO_4nn ("4nearestCO.dat",  std::ios::out | std::ios::app );
  std::ofstream outputOC_1nn ("1nearestOC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputOC_2nn ("2nearestOC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputOO_1nn ("1nearestOO.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  
  //----------------------
  // Physical parameters
  //----------------------------------------------------------------------------------
  double cut_off_radius = 1.6;             // Cut-Off for molecules
  int step = 1;                            // Step counter
  int comp_step=2;                         // The number of step you wait to compute CM
  int start_step = 5000, end_step = 30000; // Start and end step for datanalysis
  double hist_start  = 0.95;  double hist_end    = 2.00; int nb_box = 200;
  double hist_start2 = 1.00;  double hist_end2   = 3.00;
  double hist_end3   = 2.5;
  //----------------------------------------------------------------------------------

  //---------
  // Initializers
  //----------------------------------------------------------------------------------
  std::vector<Atom> atom_list;                                    // Atoms in cell
  std::vector<int> atom_indexesC; std::vector<int> atom_indexesO; // Indexes of atoms
  Contact_Matrix contact_matrix_init;                             // Initial Contact Matrix
  // Histograms
  std::vector<Bin> hist_1CO;
  std::vector<Bin> hist_2CO;
  std::vector<Bin> hist_3CO;
  std::vector<Bin> hist_4CO;
  std::vector<Bin> hist_1CC; 
  std::vector<Bin> hist_2CC;
  std::vector<Bin> hist_1OC;
  std::vector<Bin> hist_2OC; 
  std::vector<Bin> hist_1OO;
  std::string Cname = "C";
  std::string Oname = "O";
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
      // Computation at step 1
      if ( step == 1 )
	{
	  // Getting index of C and O atoms
	  atom_indexesC = makeVec(0,32);
	  atom_indexesO = makeVec(32,atom_list.size());

	  // Calculating first step contact matrix
	  contact_matrix_init = makeContactMatrix( atom_list, box );
	  // Creating Histograms
 	  hist_1CC = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesC , Cname ) , hist_start , hist_end2 , nb_box );
	  hist_2CC = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesC , Cname ) , hist_start2 , hist_end2 , nb_box );
	  hist_1CO = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesC , Oname ) , hist_start , hist_end , nb_box );
	  hist_2CO = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesC , Oname ) , hist_start , hist_end2 , nb_box );
	  hist_3CO = makeRegularHistogram( getNNearest( contact_matrix_init , 3 , atom_indexesC , Oname ) , hist_start2 , hist_end2 , nb_box );
	  hist_4CO = makeRegularHistogram( getNNearest( contact_matrix_init , 4 , atom_indexesC , Oname ) , hist_start2 , hist_end2 , nb_box );
	  hist_1OC = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesO , Cname ) , hist_start2 , hist_end2 , nb_box );
	  hist_2OC = makeRegularHistogram( getNNearest( contact_matrix_init , 2 , atom_indexesO , Cname ) , hist_start2 , hist_end2 , nb_box );
	  hist_1OO = makeRegularHistogram( getNNearest( contact_matrix_init , 1 , atom_indexesO , Oname ) , hist_start2 , hist_end2 , nb_box );
	}
      if( step % comp_step == 0 && !(debug) && step >= start_step && step <=  end_step )
	{
	  //----------------
	  // Contact Matrix
	  //-------------------------------------------------------------------
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  //-------------------------------------------------------------------
	  // Nearest Neighbours
	  //--------------------------------------------------------------------
	  hist_1CO = addHistograms( hist_1CO , makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesC , Oname ) , hist_start , hist_end , nb_box ) );
	  hist_2CO = addHistograms( hist_2CO , makeRegularHistogram( getNNearest( contact_matrix , 2 , atom_indexesC , Oname ) , hist_start , hist_end2 , nb_box ) );
	  hist_3CO = addHistograms( hist_3CO , makeRegularHistogram( getNNearest( contact_matrix , 3 , atom_indexesC , Oname ) , hist_start2 , hist_end2 , nb_box ) );
	  hist_4CO = addHistograms( hist_4CO , makeRegularHistogram( getNNearest( contact_matrix , 4 , atom_indexesC , Oname  ) , hist_start2 , hist_end2 , nb_box ) );
	  hist_1CC = addHistograms( hist_1CC , makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesC , Cname ) , hist_start , hist_end2 , nb_box ) );
	  hist_2CC = addHistograms( hist_2CC , makeRegularHistogram( getNNearest( contact_matrix , 2, atom_indexesC , Cname ) , hist_start2 , hist_end2 , nb_box ) );
	  hist_1OC = addHistograms( hist_1OC, makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesO ,  Cname ) , hist_start2 , hist_end2 , nb_box ) );
	  hist_2OC = addHistograms( hist_2OC, makeRegularHistogram( getNNearest( contact_matrix , 2 , atom_indexesO ,  Cname ) , hist_start2 , hist_end2 , nb_box ) );
	  hist_1OO = addHistograms( hist_1OO, makeRegularHistogram( getNNearest( contact_matrix , 1 , atom_indexesO ,  Oname ) , hist_start2 , hist_end2 , nb_box ) );
	  //--------------------------------------------------------------------
	  // Print step
	  //-----------------------------------------
	  std::cout << "step " << step << std::endl;
	  //-----------------------------------------
	  }
      step++;
    } while( atom_list.size() != 0 );
  //-----------------------------------------------------------------------------
  writeHistogram( outputCO_1nn , normalizeHistogram( hist_1CO ) );
  writeHistogram( outputCO_2nn , normalizeHistogram( hist_2CO ) );
  writeHistogram( outputCO_3nn , normalizeHistogram( hist_3CO ) );
  writeHistogram( outputCO_4nn , normalizeHistogram( hist_4CO ) );
  writeHistogram( outputCC_1nn , normalizeHistogram( hist_1CC ) );
  writeHistogram( outputCC_2nn , normalizeHistogram( hist_2CC ) );
  writeHistogram( outputOC_1nn , normalizeHistogram( hist_1OC ) );
  writeHistogram( outputOC_2nn , normalizeHistogram( hist_2OC ) );
  writeHistogram( outputOO_1nn , normalizeHistogram( hist_1OO ) );
    
  //Closing fluxes
  //----------------------
  input.close();
  outputCO_1nn.close();
  outputCO_2nn.close();
  outputCO_3nn.close();
  outputCO_4nn.close();
  outputCC_1nn.close();
  outputCC_2nn.close();
  outputOC_1nn.close();
  outputOC_2nn.close();
  outputOO_1nn.close();
  //----------------------
  
  return 0;
}
