#ifndef MOLECULES_H
#define MOLECULES_H

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
#include "cutoff.h"
#include "contact_matrix.h"

//==========
// MOLECULE
//=============================================
struct MoleculeBasic
{
  std::vector<std::string> names;
  std::vector<int> atom_index;
};
//=============================================


//======
// MAKE
//=====================================================================
MoleculeBasic startMolecule( std::string name_atom , int index_atom );
std::vector<MoleculeBasic> makeMolecules( ContactMatrix & cm );
//=====================================================================

//=======
// PRINT
//===========================================================
void printMolecule( MoleculeBasic molecule );
void printMoleculeSize( MoleculeBasic molecule, bool toline );
void printMolecules( std::vector<MoleculeBasic> molecules );
void printMoleculesSize( std::vector<MoleculeBasic> molecules );
//===========================================================

#endif
