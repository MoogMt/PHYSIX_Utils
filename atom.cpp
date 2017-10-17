#include "atom.h"

//===========
// DISTANCES
//======================================================================
double distanceAtoms(Atom i, Atom j)
// Returns the distance between two atoms
{
  double x2 = i.x - j.x;
  double y2 = i.y - j.y;
  double z2 = i.z - j.z;
  double dist=sqrt( pow(x2,2.0) + pow(y2,2.0) + pow(z2,2.0) );
  return dist;
}
//======================================================================

//===========
// TYPES LUT
//======================================================================
bool typeExist( std::string type, std::vector<std::string> type_vector )
// Check that a specific type is present in a vector
{
  for ( int i=0 ; i < type_vector.size() ; i++ )
    {
      if ( type_vector[i] == type )
	{
	  return true;
	}
    }
  return false;
}
//----------------------------------------------------------------
typeLUT makeLUT ( std::vector<Atom> atoms)
// Creates a LUT table for types from an atom list
{
  std::vector<std::string> types;
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      if ( ! typeExist( atoms[i].name , types ) )
	{
	  types.push_back( atoms[i].name );
	}
    }
  return { types };
}
//----------------------------------------------------------------
typeLUT makeLUT ( std::vector<Atom> atoms, int n_type)
// Creates a LUT table for types from an atom list, with extra information about the number of types that are present in the table
{
  std::vector<std::string> types;
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      if ( ! typeExist( atoms[i].name , types ) )
	{
	  types.push_back( atoms[i].name );
	}
      if ( types.size() == n_type ) break;
    }
  return { types };
}
//----------------------------------------------------------------
int getTypeId( std::string type , typeLUT type_LUT, bool msg )
// Returns the type id of a specific type
{
  for ( int i=0 ; i < type_LUT.type.size() ; i++ )
    {
      if ( type == type_LUT.type[i] ) return i;
    }
  if ( msg ) std::cout << "No such specie in the atom list!" << std::endl;
  return -1;
}
//======================================================================

//=======
// MOVE
//===========================================================================================
std::vector<Atom> compressAtoms( std::vector<Atom> atoms, double frac_a , double frac_b , double frac_c )
{
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      atoms[i].x *= frac_a;
      atoms[i].y *= frac_b;
      atoms[i].z *= frac_c;
    }
  return atoms;
}
//===========================================================================================



//=======
// IO
//======================================================================================
//-------
// PRINT
//--------------------------------------------------------------------------------------
void printAtoms( std::vector<Atom> atoms )
// Prints all informations on a given atoms in the console
{
  for( int i=0; i < atoms.size(); i++ )
    {
      std::cout << "ATOM " << i+1 << " | " << atoms[i].name << " | "  << atoms[i].x << " " << " " << atoms[i].y << " " << atoms[i].z << std::endl;
    }
  return ;
}
//--------
// WRITE
//--------------------------------------------------------------------------------------
void writePositions( std::ofstream & file , std::vector<Atom> atoms, std::string specie )
{
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      if ( atoms[i].name == specie )
	{
	  file << atoms[i].x << " ";
	  file << atoms[i].y << " ";
	  file << atoms[i].z << " ";
	  file << std::endl;
	}
    }
  file << std::endl;
  return;
}
//=======================================================================================

