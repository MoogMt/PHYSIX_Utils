// External
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

// Local Files
#include "utils.h"
#include "atom.h"
#include "cell.h"
#include "contact_matrix.h"
#include "xyz.h"


//================
// MAIN PROGRAM
//=====================================================================
int main(void)
{
  
  // Input
  //------------------------------------
  if( ! fileExists("TRAJEC.xyz",true) )
    {
      return 1;
    }
  std::ifstream input("TRAJEC.xyz");
  //------------------------------------

  // Output
  //--------------------------------------
  std::ofstream output_coord("coord.dat",std::ios::out | std::ios::app);  
  //----------------------

  //---------------------
  // Physical parameters
  //----------------------
  Cell box = {9.5,9.5,9.5,90,90,90}; // Definition of simulation box
  double cut_off_radius = 1.6;       // Cut-Off for molecules
  int step = 1;                      // Step counter
  int comp_step=1;                   // The number of step you wait to compute CM
  std::vector<Atom> atom_list;       // Atoms in cell
  //----------------------------

  // Reading XYZ file
  //---------------------------------------
  do
    {
      atom_list = readstepXYZ( input ); // Read one line
      if( step % comp_step == 0 )
	{
	  Contact_Matrix contact_matrix = makeContactMatrix(atom_list,box);
	  output_coord << step << " ";
	  writeCoordinance(output_coord,contact_matrix,"C",cut_off_radius,step,false);
	  writeCoordinance(output_coord,contact_matrix,"O",cut_off_radius,step,false);
	  output_coord << std::endl;
	  
	}
      step++;
    } while( atom_list.size() != 0 );
  //--------------------------------------

  input.close(); // Closing flux
 
  return 0;
}
