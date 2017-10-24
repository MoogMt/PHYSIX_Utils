1#ifndef LUT_H
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
//==============================
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
bool testType ( const AllTypeLut lut_all, const int atom_index ,  const std::string specie );
bool testType ( const AllTypeLUT lut_all, const int atom_index ,  const int type_index );
bool testType ( const TypeLUT lut , const std::string specie );
bool testType ( const TypeLUT lut , const int index );
//--------------------------------------------------------------------------
bool typeExists ( const std::vector<TypeLUT> list , const std::string specie );
bool typeExists( const std::vector<TypeLUT> list , const int index );
//=============================================================================================

//========
// Modify
//=============================================================================================
void addAtom2LUT( std::vector<TypeLUT>  & list , const std::string name , const int index  );
void addAtom2LUT( std::vector<TypeLUT> & list , const Atom atom );
void addAtom2LUT( AllTypeLUT & list , const std::string name , const int index  )
void addAtom2LUT( AllTypeLUT & list , const Atom atom );
//=============================================================================================

//========
// Create
//=============================================================================================
TypeLUT makeLUT( const std::string name , const int index );
TypeLUT makeLUT( const std::string name , const int index , const std::vector<int> atom_index );
std::vector<TypeLUT> makeLUT( const std::vector<Atom> atoms );
//=============================================================================================

#endif
