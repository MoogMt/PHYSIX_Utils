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
int main(void)
{
  //--------
  // Input
  //---------------------------------
  std::ifstream input("TRAJEC.xyz");
  //---------------------------------

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
  ContactMatrix cm;    // Contact Matrix
  //--------------------------------------------------
  
  //------------------------------------
  // Histograms
  //------------------------------------
  // Technical values
  double hist_start = 0.95;
  double hist_end   = 3.00;
  int nb_box = 300;
  //------------------------------------
  std::vector<Bin> hist_1CO;
  std::vector<Bin> hist_2CO;
  std::vector<Bin> hist_3CO;
  std::vector<Bin> hist_4CO;
  std::vector<Bin> hist_1CC; 
  std::vector<Bin> hist_2CC;
  std::vector<Bin> hist_1OC;
  std::vector<Bin> hist_2OC; 
  std::vector<Bin> hist_1OO;
  //---------------------------------------------------

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
  
  //--------
  // Output
  //--------------------------------------------------------------------------------
  std::ofstream outputCC_1nn ("1nearestCC.dat",  std::ios::out );
  std::ofstream outputCC_2nn ("2nearestCC.dat",  std::ios::out );
  std::ofstream outputCO_1nn ("1nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_2nn ("2nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_3nn ("3nearestCO.dat",  std::ios::out );
  std::ofstream outputCO_4nn ("4nearestCO.dat",  std::ios::out );
  std::ofstream outputOC_1nn ("1nearestOC.dat",  std::ios::out );
  std::ofstream outputOC_2nn ("2nearestOC.dat",  std::ios::out );
  std::ofstream outputOO_1nn ("1nearestOO.dat",  std::ios::out );
  //--------------------------------------------------------------------------------

  //--------------
  // Atom Indexes
  //-------------------------
  std::vector<int> atom_indexesC;
  std::vector<int> atom_indexesO;
  //-------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input , atom_list , lut_list, true, true ) )
    {
      makeContactMatrix( cm , atom_list, cell , cut_off , lut_list );
      if ( step == 1 )
	{
	  //-----------------------------------
	  // Getting atoms indexes for O and C
	  //---------------------------------------------------------------------------------
	  atom_indexesC = cm.lut_list.types[0].atom_index;
	  atom_indexesO = cm.lut_list.types[1].atom_index;
	  //---------------------------------------------------------------------------------
	  //----------
	  // CC bonds
	  //---------------------------------------------------------------------------------
	  hist_1CC = makeRegularHistogram( getNNearest( cm , 1 , atom_indexesC , "C" ) , hist_start , hist_end , nb_box );
	  hist_2CC = makeRegularHistogram( getNNearest( cm , 2 , atom_indexesC , "C" ) , hist_start , hist_end , nb_box );
	  //---------------------------------------------------------------------------------
	  // CO bonds
	  //---------------------------------------------------------------------------------
	  hist_1CO = makeRegularHistogram( getNNearest( cm , 1 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box );
	  hist_2CO = makeRegularHistogram( getNNearest( cm , 2 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box );
	  hist_3CO = makeRegularHistogram( getNNearest( cm , 3 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box );
	  hist_4CO = makeRegularHistogram( getNNearest( cm , 4 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box );
	  //---------------------------------------------------------------------------------
	  // OC bonds
	  //---------------------------------------------------------------------------------
	  hist_1OC = makeRegularHistogram( getNNearest( cm , 1 , atom_indexesO , "C" ) , hist_start , hist_end , nb_box );
	  hist_2OC = makeRegularHistogram( getNNearest( cm , 2 , atom_indexesO , "C" ) , hist_start , hist_end , nb_box );
	  //---------------------------------------------------------------------------------
	  // OO Bonds
	  //---------------------------------------------------------------------------------
	  hist_1OO = makeRegularHistogram( getNNearest( cm , 1 , atom_indexesO , "O" ) , hist_start , hist_end , nb_box );
	  //---------------------------------------------------------------------------------
	}
      else if ( step  % comp_step == 0 ) 
	{
	  //----------
	  // CC bonds
	  //---------------------------------------------------------------------------------
	  hist_1CC = addHistograms( hist_1CC , makeRegularHistogram( getNNearest( cm , 1 , atom_indexesC , "C" ) , hist_start , hist_end ,  nb_box ) );
	  hist_2CC = addHistograms( hist_2CC , makeRegularHistogram( getNNearest( cm , 2, atom_indexesC  , "C" ) , hist_start , hist_end , nb_box ) );
	  //---------------------------------------------------------------------------------
	  // CO bonds
	  //---------------------------------------------------------------------------------	  
	  hist_1CO = addHistograms( hist_1CO , makeRegularHistogram( getNNearest( cm , 1 , atom_indexesC , "O" ) , hist_start , hist_end ,  nb_box ) );
	  hist_2CO = addHistograms( hist_2CO , makeRegularHistogram( getNNearest( cm , 2 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box ) );
	  hist_3CO = addHistograms( hist_3CO , makeRegularHistogram( getNNearest( cm , 3 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box ) );
	  hist_4CO = addHistograms( hist_4CO , makeRegularHistogram( getNNearest( cm , 4 , atom_indexesC , "O" ) , hist_start , hist_end , nb_box ) );
	  //---------------------------------------------------------------------------------
	  // OC bonds
	  //---------------------------------------------------------------------------------
	  hist_1OC = addHistograms( hist_1OC , makeRegularHistogram( getNNearest( cm , 1 , atom_indexesO , "C" ) , hist_start , hist_end , nb_box ) );
	  hist_2OC = addHistograms( hist_2OC , makeRegularHistogram( getNNearest( cm , 2 , atom_indexesO , "C" ) , hist_start , hist_end , nb_box ) );
	  //---------------------------------------------------------------------------------
	  // OO Bonds
	  //---------------------------------------------------------------------------------
	  hist_1OO = addHistograms( hist_1OO , makeRegularHistogram( getNNearest( cm , 1 , atom_indexesO , "O" ) , hist_start , hist_end , nb_box ) );
	}
      step++;
     }
  //----------------------------------------------------

  //----------------------
  // Writting histograms
  //--------------------------------------------------------------------------
  writeHistogram( outputCC_1nn , normalizeHistogram( hist_1CC ) );
  writeHistogram( outputCC_2nn , normalizeHistogram( hist_2CC ) );
  writeHistogram( outputCO_1nn , normalizeHistogram( hist_1CO ) );
  writeHistogram( outputCO_2nn , normalizeHistogram( hist_2CO ) );
  writeHistogram( outputCO_3nn , normalizeHistogram( hist_3CO ) );
  writeHistogram( outputCO_4nn , normalizeHistogram( hist_4CO ) );
  writeHistogram( outputOC_1nn , normalizeHistogram( hist_1OC ) );
  writeHistogram( outputOC_2nn , normalizeHistogram( hist_2OC ) );
  writeHistogram( outputOO_1nn , normalizeHistogram( hist_1OO ) );
  //--------------------------------------------------------------------------

  //-----------------------------
  // Closing flux for histograms
  //-----------------------------------
  outputCC_1nn.close();
  outputCC_2nn.close();
  outputCO_1nn.close();
  outputCO_2nn.close();
  outputCO_3nn.close();
  outputCO_4nn.close();
  outputOC_1nn.close();
  outputOC_2nn.close();
  outputOO_1nn.close();
  //-----------------------------------
  
  //--------------
  //Closing fluxes
  //----------------------
  input.close();
   //----------------------
  
  return 0;
}
