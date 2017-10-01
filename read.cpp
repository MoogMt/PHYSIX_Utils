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
  //------------------------------------
  std::ifstream input("TRAJEC.xyz");
  std::ofstream output_contact  ("contact.dat",  std::ios::out | std::ios::app );
  std::ofstream output_nn ("nearest.dat",  std::ios::out | std::ios::app );

  //----------------------
  // Physical parameters
  //----------------------
  Cell box = {9.0,9.0,9.0,90,90,90}; // Definition of simulation box
  double cut_off_radius = 1.6;       // Cut-Off for molecules
  int step = 1;                      // Step counter
  int comp_step=5;                   // The number of step you wait to compute CM
  std::vector<Atom> atom_list;            // Atoms in cell
  std::vector<int> atom_indexes;
  std::vector<int> nearest = makeVec(1,4+1);
  //----------------------------

  // Reading XYZ file
  //---------------------------------------
  do
    {
      atom_list=readstepXYZ( input ); // Read one line
      if ( step == 1 )
	{
	  atom_indexes = makeVec(0,atom_list.size());
	}
      if( step % comp_step == 0 )
	{
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  writeAtomContact( output_contact , contact_matrix , atom_indexes );
	  writeNearest( output_nn, contact_matrix , nearest, atom_indexes , step );
	  std::cout << "step " << step << std::endl;
	}
      step++;
    } while( atom_list.size() != 0 );
  //--------------------------------------

  //closing fluxes
  input.close(); 
  output_nn.close();
  
  return 0;
}
