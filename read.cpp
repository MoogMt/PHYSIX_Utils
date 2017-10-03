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
  std::ofstream outputC_contact  ("contactC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputO_contact  ("contactO.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  std::ofstream outputC_1nn ("1nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_2nn ("2nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_3nn ("3nearestC.dat",  std::ios::out | std::ios::app );
  std::ofstream outputC_4nn ("4nearestC.dat",  std::ios::out | std::ios::app );
  //--------------------------------------------------------------------------------
  std::ofstream outputO_1nn ("1nearestO.dat",  std::ios::out | std::ios::app );
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
	}
      if( step % comp_step == 0 && !(debug) )
	{
	  //----------------
	  // Contact Matrix
	  //-------------------------------------------------------------------
	  Contact_Matrix contact_matrix = makeContactMatrix( atom_list , box );
	  //-------------------------------------------------------------------
	  // Atom Contact
	  //-------------------------------------------------------------------
	  writeAtomContact( outputC_contact , contact_matrix , atom_indexesC );
	  writeAtomContact( outputO_contact , contact_matrix , atom_indexesO );
	  //-------------------------------------------------------------------
	  // Nearest Neighbours
	  //--------------------------------------------------------------------
	  // Carbon
	  writeNearest( outputC_1nn , contact_matrix , 1 , atom_indexesC , step );
	  writeNearest( outputC_2nn , contact_matrix , 2 , atom_indexesC , step );
	  writeNearest( outputC_3nn , contact_matrix , 3 , atom_indexesC , step );
	  writeNearest( outputC_4nn , contact_matrix , 4 , atom_indexesC , step );
	  // Oxygen
	  writeNearest( outputO_1nn , contact_matrix , 1 , atom_indexesO , step );
	  writeNearest( outputO_2nn , contact_matrix , 2 , atom_indexesO , step );
	  //--------------------------------------------------------------------
	  // Print step
	  //-----------------------------------------
	  std::cout << "step " << step << std::endl;
	  //-----------------------------------------
	  }
      step++;
    } while( atom_list.size() != 0 );
  //-----------------------------------------------------------------------------

  //Closing fluxes
  //----------------------
  input.close();
  outputC_contact.close();
  outputO_contact.close();
  outputC_1nn.close();
  outputC_2nn.close();
  outputC_3nn.close();
  outputC_4nn.close();
  outputO_1nn.close();
  outputO_2nn.close();
  //----------------------
  
  return 0;
}
