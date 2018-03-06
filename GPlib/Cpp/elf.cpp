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
  std::ifstream input("ELF.cube");
  //--------------------------------

  //--------
  // Output
  //-------------------------------------
  std::ofstream output("ELF_out.dat");
  //-------------------------------------
  
  //---------------
  // Initializers
  //-------------------------------------
  AtomList  atom_list; // Atoms in cell
  //-------------------------------------

  //--------------------
  // Reading Cell File
  //-------------------------------------------------------------------
  Cell cell;
  if ( ! readParamCell( "cell.param" , cell ) )
    {
      return 1;
    }
  //-------------------------------------------------------------------

  std::string line;
  int nb_atoms=0;
  double ox, oy, oz;
  int nb_vox[3];
  double x,y,z;
  int atom_type;
  double dummy;
  double cell_x[3], cell_y[3], cell_z[3];
  AtomList atoms;
  
  //-------------------------------
  // Read useless two first lines
  //-------------------------------------------------------------------
  for ( int i=0; i<2; i++)
    {
  
      if ( ! std::getline(input,line) )
	{
	  std::cout << "Problem with the file ELF.cube" <<std::endl;
	  return 1;
	}
    }
  //-------------------------------------------------------------------

  //----------------------------------------------------------------------
  // Reading number of atoms and origin of the origin point of the density
  //------------------------------------------------------------------------
  if( std::getline(input,line) )
    {
      std::istringstream it_string(line);
      if ( ! ( it_string >> nb_atoms >> ox >> oy >> oz ) ) return 1;
    }
  else return 1;
  //-----------------------------------------------------------------------

  //--------------------------------------------
  // Reading number of voxels and cell vectors 
  //-----------------------------------------------------------------------
  for( int i=0 ; i<3 ; i++ )
    {
      if( std::getline(input,line) )
	{
	  if ( ! ( it_string >> nb_vox[i] >> cell_x[i] >> cell_y[i] >> cell_z[i] ) ) return 1;
	}
      else return 1;
    }
  //-----------------------------------------------------------------------

  //---------------------------
  // Reading atomic positions
  //----------------------------------------------------------------------------------------
  // Initiate the atom_list
  //-------------------------------------
  atoms.assign( nb_atoms, "");
  atoms.x.assign( nb_atoms, 0);
  atoms.y.assign( nb_atoms, 0);
  atoms.z.assign( nb_atoms, 0);
  atoms.index.assign( nb_atoms, 0);
  //-------------------------------------
  // Reading the atomic positions and all
  //----------------------------------------------------------------------------------
  for( int i=0 ; i < nb_atoms ; i++ )
    {
      if ( std::getline(input,line) )
	{
	  if ( ! ( it_string >> atom_type >> dummy >> x >> y >> z ) ) return 1;
	}
      else return 1;
    }
  //----------------------------------------------------------------------------------------
  
  return 0;
}
