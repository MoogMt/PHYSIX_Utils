#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

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
  std::ifstream input("TRAJEC.xyz");

  //----------------------
  // Physical parameters
  //----------------------
  Cell box = {9.0,9.0,9.0,90,90,90}; // Definition of simulation box
  double cut_off_radius = 1.6;      // Cut-Off for molecules
  int step = 1;                      // Step counter
  int comp_step=5;                   // The number of step you wait to compute CM

  std::vector<Atom> atom_list;            // Atoms in cell
  //----------------------------

  // Reading XYZ file
  //---------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if( step % comp_step == 0 )
	{
	  Contact_Matrix contact_matrix = makeContactMatrix(atom_list,box);
	  std::vector<int> coordO  = getCoordinances("O",contact_matrix,cut_off_radius);
	  std::vector<int> coordC  = getCoordinances("C",contact_matrix,cut_off_radius);
	  std::cout << step << " ";
	  for ( int i=0; i < coordC.size() ; i++ )
	    {
	      std::cout << coordC[i] << " ";
	    }
	  for ( int i=0; i < coordO.size() ; i++ )
	    {
	      std::cout << coordO[i] << " ";
	    }
	  std::cout << average(coordC) << " ";
	  std::cout << average(coordO) << " ";
	  std::cout << std::endl;
	}
      step++;
    } while( atom_list.size() != 0 );
  //--------------------------------------

  input.close(); // Closing flux
 
  return 0;
}
