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

//-------
// CELL
//---------------------------
struct Bond
{
  int atom1_index, atom2_index;
  double strengh;
};
//----------------------------

//==========
// FUNCTIONS
//========================================

std::vector<Bond> makeBonds(std::vector<Atom> atom_list, double cut_off_radius ) 
std::vector<Bond> makeBonds(Contact_Matrix contact_matrix , double cut_off_radius )
//=======================================

#endif 
