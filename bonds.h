#ifndef BONDS_H
#define BONDS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"

//======
// BOND
//===============================
struct Bond
{
  int atom1_index, atom2_index;
  double strengh;
};
//===============================

//==========
// FUNCTIONS
//========================================
Bond emptyBond();
Bond makeBond( atom_list[i] , atom_list[j] );
std::vector<Bond> makeBonds(std::vector<Atom> atom_list, double cut_off_radius ) 
std::vector<Bond> makeBonds(Contact_Matrix contact_matrix , double cut_off_radius )
//=======================================

#endif 
