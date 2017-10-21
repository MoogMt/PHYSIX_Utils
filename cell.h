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
// PBC
//----------------------------------------------------------------------------------------------
double backIn( double x, double a );
Atom wrapPBC(Atom atom_in, Cell box );
std::vector<Atom> wrapPBC( std::vector<Atom> atoms , Cell cell );
void wrapPBC( AtomList & atom_list , int i , Cell cell );
void wrapPBC( AtomList & atom_list , Cell cell );
std::vector<Atom> pbcImages( Atom atom, Cell box );
//----------------------------------------------------------------------------------------------
// Calculating Distances
//-----------------------------------------------------------------------------------------------
double distAtoms1D( double x1, double x2, double a );
double distanceAtoms( std::vector<Atom> atoms, int i, int j, Cell box , bool wrap = true , bool sqrt_test = true );
double distanceAtomsSq( AtomList & atom_list ,  int i , int  j , Cell cell , bool wrap );
//----------------------------------------------------------------------------------------------
// Change Box
//----------------------------------------------------------------------------------------------
Cell compressBox( Cell cell , double frac_a , double frac_b , double frac_c );
//----------------------------------------------------------------------------------------------
// IO
//----------------------------------------------------------------------------------------------
// Writting Atom distances
void writeAtomDistances( std::ofstream & file , std::vector<Atom> atom_list , std::vector<int> atom_index, Cell box);
// Reading Cell Param File
Cell readParamCellStep( std::ofstream & file );
Cell readParamCell( std::string file_name );
//=======================================

#endif 
