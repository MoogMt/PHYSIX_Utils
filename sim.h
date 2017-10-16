#ifndef SIM_H
#define SIM_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"
#include "cell.h"
#include "utils.h"

//=====
// SIM
//===========================
struct Sim
{
  std::vector<Atom> atoms;
  Cell cell;
};
//===========================

//==========
// Function
//===========================
Sim emptySim();
Sim compressBox( Sim sim_set , double frac_a , double frac_b , double frac_c );
//===========================

#endif 
