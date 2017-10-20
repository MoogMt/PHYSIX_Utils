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
bool readStepXYZ( std::ifstream & file , std::vector<Atom> & atom , std::vector<typeLUT> & lut_list , bool verbose = false );
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
