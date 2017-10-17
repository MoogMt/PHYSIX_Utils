#include "atom.h"

//===========
// DISTANCES
//======================================================================
double distanceAtoms(Atom i, Atom j)
// Distance between two atoms
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
bool typeExist( std::string type, std::string type_vector[] , int size )
{
  for ( int i=0 ; i < size ; i++ )
    {
      if ( type_vector[i] == type )
	{
	  return true;
	}
    }
  return false;
}
//----------------------------------------------------------------
template <int N> typeLUT<N> makeLUT ( std::vector<Atom> atoms)
{
  int count_type = 0;
  std::string types[N];
  for ( int i=0 ; i < atoms.size() ; i++ )
    {
      if ( ! typeExist( atoms[i].name , types , count_type ) )
	{
	  types[i] =  atoms[i].name;
	  count_type++;
	}
      if ( count_type == N ) break;
    }
  return { types };
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
// PRINT
//======================================================================
void printAtoms( std::vector<Atom> atoms )
// Prints all informations on a given atoms in the console
{
  for( int i=0; i < atoms.size(); i++ )
    {
      std::cout << "ATOM " << i+1 << " | " << atoms[i].name << " | "  << atoms[i].x << " " << " " << atoms[i].y << " " << atoms[i].z << std::endl;
    }
  return ;
}
//======================================================================

//=======
// WRITE
//=======================================================================================
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

