#ifndef XYZ_H
#define XYZ_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>

#include "atom.h"

std::vector<Atom> readstepXYZ(std::ifstream& file);

#endif 
