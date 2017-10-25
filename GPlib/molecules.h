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
#include "lut.h"

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
MoleculeBasic startMolecule( const std::string name_atom , const int index_atom );
std::vector<MoleculeBasic> makeMolecules( const ContactMatrix & cm );
//=====================================================================

//=======
// PRINT
//===========================================================
void printMolecule( const MoleculeBasic molecule );
void printMoleculeSize( const MoleculeBasic molecule, const bool toline );
void printMolecules( const std::vector<MoleculeBasic> molecules );
void printMoleculesSize( const std::vector<MoleculeBasic> molecules );
//===========================================================

#endif
