#ifndef CELL_H
#define CELL_H

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

//-------
// CELL
//---------------------------
struct Cell
{
  double a, b, c;             // lengths parameters
  double alpha, beta, gamma;  // angular parameters
};
//----------------------------

//==========
// FUNCTIONS
//========================================
double backIn( double x, double a );
Atom wrapPBC(Atom atom_in, Cell box );
std::vector<Atom> pbc( Atom atom, Cell box );
double distanceAtoms( std::vector<Atom> atoms , int i , int j , Cell box );
Cell readParamCellStep( std::ofstream & file );
Cell readParamCell( std::string file_name );
//=======================================

#endif 
