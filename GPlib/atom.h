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

//===========
// ATOM LIST
//================================
struct AtomList
{
  std::vector<std::string> names;
  std::vector<double> x, y, z;
  std::vector<int> index;
};
//================================
  
//=============================
// TYPES LUT ( Look Up Table )
//=============================================
struct TypeLUT
{
  std::string type_name;
  int type_index;
  std::vector<int> atom_index;
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
// Tests
//---------------------------------------------------------------
bool testType ( TypeLUT lut , std::string specie );
bool testType ( TypeLUT lut , int index );
bool typeExists ( std::vector<TypeLUT> list , std::string specie );
bool typeExists( std::vector<TypeLUT> list , int index );
//---------------------------------------------------------------
// Modify
//---------------------------------------------------------------
void addAtom2LUT( std::vector<TypeLUT>  & list , std::string name , int index  );
void addAtom2LUT( std::vector<TypeLUT> & list , Atom atom );
//---------------------------------------------------------------
// Create
//---------------------------------------------------------------
TypeLUT makeLUT( std::string name , int index );
TypeLUT makeLUT( std::string name , int index , std::vector<int> atom_index );
std::vector<TypeLUT> makeLUT( std::vector<Atom> atoms );
//=============================================================================================

//=======
// MOVE 
//=============================================================================================
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c );
//=============================================================================================

//=====
// IO
//=============================================================================================
// WRITE
//-----------------------------------------------------------------------------------------
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie );
//=============================================================================================

#endif 
