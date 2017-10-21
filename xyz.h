#ifndef XYZ_H
#define XYZ_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"

//======
// READ
//=======================================================================
bool readStepXYZ( std::ifstream & file , AtomList & atoms , std::vector<TypeLUT> & lut_list , bool same_type=true , bool verbose = true );
bool readStepXYZ( std::ifstream & file , std::vector<Atom> & atoms , std::vector<TypeLUT> & lut_list , bool same_type = true , bool verbose = true );
//----------------------------------------------------
std::vector<Atom> readstepXYZ( std::ifstream & file );
//=======================================================================

//=======
// WRITE
//=======================================================================
void writeXYZ( std::ofstream & file , std::vector<Atom> atom_list );
void writeXYZ( std::ofstream & file , std::vector<Atom> atom_list , int step );
//=======================================================================
  
#endif 
