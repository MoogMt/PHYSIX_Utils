#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>  
#include <stdlib.h>    
#include <math.h>


// Min 
double min(std::vector<double> vector);
int min( int int1, int int2 );

// Max
int max( int int1, int int2 );

// Averages
double average(std::vector<int> data);
double average(std::vector<double> data);

// Sums
int sumFromO(int integer);
int sumBtw(int int1, int int2);

// Compute separator
int computeSep(int atom_index, int nb_atoms);

// String related
bool fileExists(const std::string file_name, const bool message );

#endif 
