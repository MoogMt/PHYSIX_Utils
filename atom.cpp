
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
bool testType( typeLUT lut , std::string specie )
{
  if ( lut.type_name == specie ) return true;
  else return false;
}
//---------------------------------------------------------------
bool testType( typeLUT lut , int index )
{
  if ( lut.type_index == index ) return true;
  else return false;
}
//---------------------------------------------------------------
bool typeExists( std::vector<typeLUT> list , std::string specie )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], specie ) ) return true;
    }
  return false;
}
//---------------------------------------------------------------
bool typeExists( std::vector<typeLUT> list , int index )
{
  for ( int i = 0 ; i < list.size() ; i++ )
    {
      if ( testType( list[i], index ) ) return true;
    }
  return false;
}
//======================================================================

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
//----------------------------------------------------------------------
