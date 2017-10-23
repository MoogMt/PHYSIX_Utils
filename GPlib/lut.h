#ifndef LUT_H
#define LUT_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"

//=====
// LUT
//===============================
struct AllTypeLUT
{
  std::vector<int> type_index;
  std::vector<std::string> type_name;
  std::vector<TypeLUT> types;
};
struct TypeLUT
{
  std::string type_name;
  int index;
  std::vector<int> atom_index;
};
//===============================

//=======
// Tests
//=============================================================================================
bool testType ( TypeLUT lut , std::string specie );
bool testType ( TypeLUT lut , int index );
bool typeExists ( std::vector<TypeLUT> list , std::string specie );
bool typeExists( std::vector<TypeLUT> list , int index );
//=============================================================================================

//========
// Modify
//=============================================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , std::string name , int index  );
void addAtom2LUT( std::vector<TypeLUT> & list , Atom atom );
//=============================================================================================

//========
// Create
//=============================================================================================
TypeLUT makeLUT( std::string name , int index );
TypeLUT makeLUT( std::string name , int index , std::vector<int> atom_index );
std::vector<TypeLUT> makeLUT( std::vector<Atom> atoms );
//=============================================================================================

#endif
