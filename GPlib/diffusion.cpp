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
  std::ifstream input40("40GPa/TRAJEC.xyz");
  std::ifstream input45("45GPa/TRAJEC.xyz");
  std::ifstream input50("50GPa/TRAJEC.xyz");
  std::ifstream input60("60GPa/TRAJEC.xyz");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream diffusion40("diffusion_40.dat");
  std::ofstream diffusion45("diffusion_45.dat");
  std::ofstream diffusion50("diffusion_50.dat");
  std::ofstream diffusion60("diffusion_60.dat");
  //-------------------------------------
  
  //----------------------
  // Physical parameters
  //--------------------------------------
  int step       = 1;  // Step counter
  int start_step = 5000; // Start step
  int end_step   = 120000;
  int comp_step  = 1; // Frequency of computation
  //--------------------------------------

  //---------------
  // Initializers
  //--------------------------------------------------
  AtomList  atom_list;     // Atoms in cell
  std::vector<double> r; // Atoms_old in cell
  AllTypeLUT lut_list; // LUT for types
  std::vector<double> x0, y0, z0;
  std::vector<double> x, y, z;
  double d_40=0, d_45=0, d_50=0 , d_60=0;
  int count = 0;
  //--------------------------------------------------

  //-----------------
  // Reading Cut-Off
  //-------------------------------------------------------------------
  CutOffMatrix cut_off;
  if ( ! readCutOff( "40GPa/cut_off.dat" , cut_off , lut_list ) )
    {
      return 1;
    }
  //------------------------------------------------------------------

  //------
  // Cell
  //------------------------------------------------------------------
  Cell cell40;
  if ( ! readParamCell( "40GPa/cell.param" , cell40 ) )
    {
      return 1;
    }
  Cell cell45;
  if ( ! readParamCell( "45GPa/cell.param" , cell45 ) )
    {
      return 1;
    }
  Cell cell50;
  if ( ! readParamCell( "50GPa/cell.param" , cell50 ) )
    {
      return 1;
    }
  Cell cell60;
  if ( ! readParamCell( "60GPa/cell.param" , cell60 ) )
    {
      return 1;
    }
  //------------------------------------------------------------------

  //----------------------------------------------------
  computeDiff( diffusion40 , input40 , comp_step , start_step , end_step , atom_list , lut_list , cell40 );
  //computeDiff( diffusion45 , input45 , comp_step , start_step , end_step , atom_list , lut_list , cell45 );
  //computeDiff( diffusion50 , input50 , comp_step , start_step , end_step , atom_list , lut_list , cell50 );
  //computeDiff( diffusion60 , input60 , comp_step , start_step , end_step , atom_list , lut_list , cell60 );
  //----------------------------------------------------
  
  //---------------
  //Closing fluxes
  //----------------------
  diffusion40.close();
  diffusion45.close();
  diffusion50.close();
  diffusion60.close();
  //----------------------
      
  return 0;
}
