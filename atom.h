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

//-------
// ATOM
//------------------------------
struct Atom
{
  std::string name;    // name
  double x, y, z; // position
  int index;      // index
};
//------------------------------

//-------
// BOND
//--------------------------------
struct Bond
{
  int index_i, index_j;  // Index of atoms i and j
  int time;              // Nb of steps the bond lives
};
//--------------------------------

//-----------
// MOLECULE
//--------------------
struct Molecule
{
  std::string name;        // Name of molecule
  std::vector<Atom> atoms; // Atoms
};
//-------------------


//============
// FUNCTIONS
//================================================
void printAtoms( std::vector<Atom> atoms );
double distanceAtoms(Atom i, Atom j);
//================================================

#endif 
