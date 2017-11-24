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
  std::ofstream diffcoef("diffcoef.dat");
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
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //------------------------------------------------------------------
  
  //-------------------
  // Reading XYZ file
  //----------------------------------------------------
  while( readStepXYZfast( input40 , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = backIn( atom_list.x , cell.a );
	  y0 = backIn( atom_list.y , cell.b );
	  z0 = backIn( atom_list.z , cell.c );
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( backIn( atom_list.x, cell.a ) , x0 );
	  std::vector<double> y = difference( backIn( atom_list.y, cell.b ) , y0 );
	  std::vector<double> z = difference( backIn( atom_list.z, cell.c ) , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  diffusion40 << step-start_step << " " << average( r ) << std::endl;
	  count++;
	  if ( step > start_step + 2000 )
	    {
	      d_40 += average( r );
	    }
	}      
      std::cout << step << std::endl;
      step++;
     }
  d_40 /= (double)(count);
  //----------------------------------------------------
  step=0;
  //----------------------------------------------------
  while( readStepXYZfast( input45 , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = backIn( atom_list.x, cell.a );
	  y0 = backIn( atom_list.y, cell.b );
	  z0 = backIn( atom_list.z, cell.c );
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( backIn( atom_list.x, cell.a ) , x0 );
	  std::vector<double> y = difference( backIn( atom_list.y, cell.b ) , y0 );
	  std::vector<double> z = difference( backIn( atom_list.z, cell.c ) , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  diffusion45 << step-start_step << " " << average( r ) << std::endl;
	  count++;
	  if ( step > start_step + 2000 )
	    {
	      d_45 += average( r );
	    }
	}      
      std::cout << step << std::endl;
      step++;
     }
  d_45 /= (double)(count);
  //----------------------------------------------------
  step=0;
  //----------------------------------------------------
  while( readStepXYZfast( input50 , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = backIn( atom_list.x, cell.a );
	  y0 = backIn( atom_list.y, cell.b );
	  z0 = backIn( atom_list.z, cell.c );
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( backIn( atom_list.x , cell.a ) , x0 );
	  std::vector<double> y = difference( backIn( atom_list.y , cell.b ) , y0 );
	  std::vector<double> z = difference( backIn( atom_list.z , cell.c ) , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  diffusion50 << step-start_step << " " << average( r ) << std::endl;
	  count++;
	  if ( step > start_step + 2000 )
	    {
	      d_50 += average( r );
	    }
	}      
      std::cout << step << std::endl;
      step++;
     }
  d_50 /= (double)(count);
  //----------------------------------------------------
  step=0;
  //----------------------------------------------------
  while( readStepXYZfast( input60 , atom_list , lut_list, true, true ) )
    {
      if ( step == start_step )
	{
	  x0 = backIn( atom_list.x , cell.a );
	  y0 = backIn( atom_list.y , cell.b );
	  z0 = backIn( atom_list.z , cell.c );
	}
      else if ( step % comp_step == 0 && step > start_step && step < end_step )
	{
	  std::vector<double> x = difference( backIn( atom_list.x , cell.a ) , x0 );
	  std::vector<double> y = difference( backIn( atom_list.y , cell.b ) , y0 );
	  std::vector<double> z = difference( backIn( atom_list.z , cell.c ) , z0 );
	  std::vector<double> r = square( squaroot( addVector( addVector( square( x ), square( y ) ), square( z ) ) ) );
	  diffusion60 << step-start_step << " " << average( r ) << std::endl;
	  count++;
	  if ( step > start_step + 2000 )
	    {
	      d_60 += average( r );
	    }
	}      
      std::cout << step << std::endl;
      step++;
     }
  d_60 /= (double)(count);
  //----------------------------------------------------

  //----------------------------------------------------
  diffcoef << 40 << " " << d_40 << std::endl;
  diffcoef << 45 << " " << d_45 << std::endl;
  diffcoef << 50 << " " << d_50 << std::endl;
  diffcoef << 60 << " " << d_60 << std::endl;
  //----------------------------------------------------
  
  //---------------
  //Closing fluxes
  //----------------------
  diffusion40.close();
  diffusion45.close();
  diffusion50.close();
  diffusion60.close();
  diffcoef.close();
  //----------------------
      
  return 0;
}
