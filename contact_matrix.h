#ifndef CONTACT_MATRIX_H
#define CONTACT_MATRIX_H

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

//----------------
// CONTACT MATRIX
//---------------------------------------------
struct Contact_Matrix
{
  std::vector<std::string> types;  // Types of the atoms
  std::vector<double> matrix; // Contact matrix
};
//---------------------------------------------

//===========
// FUNCTIONS
//====================================================================
// Computes the contact matrix from an atom list and a cell
Contact_Matrix makeContactMatrix(std::vector<Atom> atom_list, Cell box);
// Gets the list of distances of all atoms with regard to an atom
std::vector<double> getAtomContact(Contact_Matrix contact_matrix, int atom_index);
// Gets the coordinance of an atom
std::vector<int> getCoordinances(std::string type, Contact_Matrix contact_matrix, double cut_off_radius);
//====================================================================

#endif 
