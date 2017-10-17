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
template <int N>
struct typeLUT
{
  std::string types[N];
};
//=============================================

//============
// FUNCTIONS
//=============================================================================================
double distanceAtoms(Atom i, Atom j);
bool typeExist( std::string type, std::string type_vector[] , int size );
template <int N> typeLUT<N> makeLUT ( std::vector<Atom> atoms);
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c );
void printAtoms( std::vector<Atom> atoms );
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
//=============================================================================================

#endif 
