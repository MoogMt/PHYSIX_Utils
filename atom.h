#ifndef ATOMS_H
#define ATOMS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

//======
// ATOM
//==============================
struct Atom
{
  std::string name;    // name
  double x, y, z; // position
  int index;      // index
};
//==============================

//==========
// MOLECULE
//=============================================
struct Molecule
{
  std::string name;        // Name of molecule
  std::vector<Atom> atoms; // Atoms
};
//=============================================

//===========
// TYPES LUT
//=============================================
struct indexes
{
  std::vector<int> list;
};
struct typeLUT
{
  std::vector<std::string> type;
  std::vector<std::list> type_lists;
};
//=============================================

//===========
// Distances
//=============================================================================================
double distanceAtoms(Atom i, Atom j);
//=============================================================================================

//=====
// LUT
//=============================================================================================
bool typeExist( std::string type, std::string type_vector[] , int size );
typeLUT makeLUT ( std::vector<Atom> atoms);
typeLUT makeLUT ( std::vector<Atom> atoms, int n_type);
int getTypeId( std::string specie , typeLUT type_LUT , bool msg = false );
//=============================================================================================

//=======
// MOVE 
//=============================================================================================
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c );
//=============================================================================================

//=====
// IO
//=============================================================================================

void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
//=============================================================================================
double distanceAtoms(Atom i, Atom j);
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c );
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
//================================================

#endif 
