#ifndef CUTOFF_H
#define CUTOFF_H

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

//================
// Cut-Off Matrix
//====================================================================
struct CutOffMatrix
{
  std::vector<TypeLUT> list_lut;
  std::vector<double> matrix;
};
//====================================================================

//=============
// FUNCTIONS
//=======================================================================================
// Read Cut Off Matrix
//---------------------------------------------------------------------------------------
CutOffMatrix readCutOff( std::string file , std::vector<TypeLUT> & lut_list );
//-------------
// Get Cut Off 
//---------------------------------------------------------------------------------------
double getCutOff( CutOffMatrix cut_off_matrix , int i , int j );
//=======================================================================================

#endif 
